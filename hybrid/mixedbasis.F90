!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine generates the mixed basis set used to evaluate the  !
! exchange term in HF/hybinp functional calculations or EXX           !
! calculations. In the latter case a second mixed basis set is setup  !
! for the OEP integral equation.                                      !
! In all cases the mixed basis consists of IR plane waves             !
!                                                                     !
! IR:                                                                 !
!    M_{\vec{k},\vec{G}} = 1/\sqrt{V} \exp{i(\vec{k}+\vec{G})}        !
!                                                                     !
! which are zero in the MT spheres and radial functions times         !
! spherical harmonics in the MT spheres                               !
!                                                                     !
! MT:                                                                 !
!     a                a                                              !
!    M              = U   Y                                           !
!     PL               PL  LM                                         !
!                                                                     !
!            where     a    a  a                                      !
!                     U  = u  u                                       !
!                      PL   pl  p'l'                                  !
!                                                                     !
!               and    L \in {|l-l'|,...,l+l'}                        !
!                                                                     !
!               and    P counts the different combinations of         !
!                      pl,p'l' which contribute to L                  !
!                                                                     !
!                                               M.Betzinger (09/07)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_mixedbasis

CONTAINS

   SUBROUTINE mixedbasis(atoms, kpts, input, cell, xcpot, mpinp, mpdata, hybinp, hybdat,&
                         enpara, fmpi, v, iterHF)

      USE m_judft
      USE m_types
      USE m_constants
      USE m_loddop, ONLY: loddop
      USE m_intgrf, ONLY: intgrf_init, intgrf
      use m_rorder, only: rorderpf
      USE m_hybrid_core
      USE m_wrapper
      USE m_eig66_io

      IMPLICIT NONE

      TYPE(t_xcpot_inbuild), INTENT(IN)    :: xcpot
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      TYPE(t_mpdata), intent(inout)  :: mpdata
      TYPE(t_mpinp), intent(in)     :: mpinp
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat
      TYPE(t_enpara), INTENT(IN)    :: enpara
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_potden), INTENT(IN)    :: v

      integer, intent(in) :: iterHF

      ! local type variables
      TYPE(t_usdus)                   ::  usdus

      ! local scalars
      INTEGER                         ::  jspin, itype, l1, l2, l, n_radbasfn, full_n_radbasfn, n1, n2
      INTEGER                         ::  i_basfn, i, n_grid_pt,j
      REAL                            ::  rdum, rdum1, max_momentum, momentum

      ! - local arrays -

      REAL                            ::  bashlp(atoms%jmtd)

      REAL, ALLOCATABLE               ::  bas1(:,:,:,:,:), bas2(:,:,:,:,:)
      REAL, ALLOCATABLE               ::  basmhlp(:,:,:,:)
      REAL, ALLOCATABLE               ::  gridf(:,:), vr0(:,:,:)

      LOGICAL, ALLOCATABLE            ::  selecmat(:,:,:,:)
      LOGICAL, ALLOCATABLE            ::  seleco(:,:), selecu(:,:)

      CHARACTER, PARAMETER            :: lchar(0:38) = (/'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', &
                                                         'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', &
                                                         'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x'/)


      IF (fmpi%irank == 0) WRITE (oUnit, '(//A,I2,A)') '### subroutine: mixedbasis ###'

      IF (xcpot%is_name("exx")) CALL judft_error("EXX is not implemented in this version", calledby='mixedbasis')

      ! Deallocate arrays which might have been allocated in a previous run of this subroutine
      IF (ALLOCATED(mpdata%n_g)) deallocate(mpdata%n_g)
      IF (ALLOCATED(mpdata%num_radbasfn)) deallocate(mpdata%num_radbasfn)
      IF (ALLOCATED(mpdata%gptm_ptr)) deallocate(mpdata%gptm_ptr)
      IF (ALLOCATED(mpdata%g)) deallocate(mpdata%g)
      IF (ALLOCATED(mpdata%radbasfn_mt)) deallocate(mpdata%radbasfn_mt)

      CALL usdus%init(atoms, input%jspins)

      if(.not. allocated(mpdata%num_radfun_per_l)) &
         allocate(mpdata%num_radfun_per_l(0:atoms%lmaxd, atoms%ntype), source=-1)

      call mpdata%set_num_radfun_per_l(atoms)
      call mpdata%init(hybinp, hybdat, atoms)

      ! initialize gridf for radial integration
      CALL intgrf_init(atoms%ntype, atoms%jmtd, atoms%jri, atoms%dx, atoms%rmsh, gridf)

      allocate(vr0(atoms%jmtd, atoms%ntype, input%jspins), source=0.0)

      vr0(:,:,:) = v%mt(:,0, :,:)

      ! calculate radial basisfunctions u and u' with
      ! the spherical part of the potential vr0 and store them in
      ! bas1 = large component ,bas2 = small component

      call gen_bas_fun(atoms, enpara, gridf, input, mpdata, fmpi, vr0, usdus, bas1, bas2)

      ! - - - - - - SETUP OF THE MIXED BASIS IN THE IR - - - - - - -

      ! construct G-vectors with cutoff smaller than gcutm
      call mpdata%gen_gvec(mpinp, cell, kpts, fmpi)

      ! - - - - - - - - Set up MT product basis for the non-local exchange potential  - - - - - - - - - -

      IF (fmpi%irank == 0) THEN
         WRITE (oUnit, '(A)') 'MT product basis for non-local exchange potential:'
         WRITE (oUnit, '(A)') 'Reduction due to overlap (quality of orthonormality, should be < 1.0E-06)'
      END IF

      allocate(mpdata%num_radbasfn(0:maxval(hybinp%lcutm1), atoms%ntype), source=0)
      allocate(seleco(maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd))
      allocate(selecu(maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd))

      ! determine maximal indices of (radial) mixed-basis functions (->num_radbasfn)
      ! (will be reduced later-on due to overlap)
      DO itype = 1, atoms%ntype
         seleco = .FALSE.
         selecu = .FALSE.
         seleco(1, 0:hybinp%select1(1, itype)) = .TRUE.
         selecu(1, 0:hybinp%select1(3, itype)) = .TRUE.
         seleco(2, 0:hybinp%select1(2, itype)) = .TRUE.
         selecu(2, 0:hybinp%select1(4, itype)) = .TRUE.

         ! include local orbitals
         IF (maxval(mpdata%num_radfun_per_l) >= 3) THEN
            seleco(3:,:) = .TRUE.
            selecu(3:,:) = .TRUE.
         END IF

         DO l = 0, hybinp%lcutm1(itype)
            n_radbasfn = 0

            !
            ! valence * valence
            !
            if(.not. allocated(selecmat)) then
               allocate(selecmat(maxval(mpdata%num_radfun_per_l), &
                                 0:atoms%lmaxd, &
                                 maxval(mpdata%num_radfun_per_l), &
                                 0:atoms%lmaxd))
            endif
            selecmat = calc_selecmat(atoms, mpdata, seleco, selecu)

            DO l1 = 0, atoms%lmax(itype)
               DO l2 = 0, atoms%lmax(itype)
                  IF (l >= ABS(l1 - l2) .AND. l <= l1 + l2) THEN
                     DO n1 = 1, mpdata%num_radfun_per_l(l1, itype)
                        DO n2 = 1, mpdata%num_radfun_per_l(l2, itype)
                           IF (selecmat(n1, l1, n2, l2)) THEN
                              n_radbasfn = n_radbasfn + 1
                              selecmat(n2, l2, n1, l1) = .FALSE. ! prevent double counting of products (a*b = b*a)
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
            IF (n_radbasfn == 0 .AND. fmpi%irank == 0) &
               WRITE (oUnit, '(A)') 'mixedbasis: Warning!  No basis-function product of '//lchar(l)// &
               '-angular momentum defined.'
            mpdata%num_radbasfn(l, itype) = n_radbasfn*input%jspins
         END DO
      END DO

      allocate(mpdata%radbasfn_mt(atoms%jmtd,&
                            maxval(mpdata%num_radbasfn), &
                            0:maxval(hybinp%lcutm1), &
                            atoms%ntype), source=0.0)

      ! Define product bases and reduce them according to overlap

      DO itype = 1, atoms%ntype
         seleco = .FALSE.
         selecu = .FALSE.
         seleco(1, 0:hybinp%select1(1, itype)) = .TRUE.
         selecu(1, 0:hybinp%select1(3, itype)) = .TRUE.
         seleco(2, 0:hybinp%select1(2, itype)) = .TRUE.
         selecu(2, 0:hybinp%select1(4, itype)) = .TRUE.
         ! include lo's
         IF (maxval(mpdata%num_radfun_per_l) >= 3) THEN
            seleco(3:,:) = .TRUE.
            selecu(3:,:) = .TRUE.
         END IF

         n_grid_pt = atoms%jri(itype)
         DO l = 0, hybinp%lcutm1(itype)
            full_n_radbasfn = mpdata%num_radbasfn(l, itype)
            ! allow for zero product-basis functions for
            ! current l-quantum number
            IF (n_radbasfn == 0) THEN
               IF (fmpi%irank == 0) WRITE (oUnit, '(6X,A,'':   0 ->   0'')') lchar(l)
               CYCLE
            END IF

            ! set up the overlap matrix
            i_basfn = 0

            ! valence*valence
            selecmat =  calc_selecmat(atoms, mpdata,seleco, selecu)

            DO l1 = 0, atoms%lmax(itype)
               DO l2 = 0, atoms%lmax(itype)
                  IF (l >= ABS(l1 - l2) .AND. l <= l1 + l2) THEN
                     DO n1 = 1, mpdata%num_radfun_per_l(l1, itype)
                        DO n2 = 1, mpdata%num_radfun_per_l(l2, itype)

                           IF (selecmat(n1, l1, n2, l2)) THEN
                              DO jspin = 1, input%jspins
                                 i_basfn = i_basfn + 1
                                 IF (i_basfn > full_n_radbasfn) call judft_error('got too many product functions', hint='This is a BUG, please report', calledby='mixedbasis')

                                 mpdata%radbasfn_mt(:n_grid_pt, i_basfn, l, itype) &
                                    = (   bas1(:n_grid_pt, n1, l1, itype, jspin) &
                                        * bas1(:n_grid_pt, n2, l2, itype, jspin) &
                                        + bas2(:n_grid_pt, n1, l1, itype, jspin) &
                                        * bas2(:n_grid_pt, n2, l2, itype, jspin) &
                                      ) / atoms%rmsh(:n_grid_pt, itype)


                              END DO !jspin
                              ! prevent double counting of products (a*b = b*a)
                              selecmat(n2, l2, n1, l1) = .FALSE.
                           END IF
                        END DO !n2
                     END DO !n1
                  ENDIF
               END DO !l2
            END DO  !l1


            IF (i_basfn /= full_n_radbasfn) call judft_error('counting error for product functions', hint='This is a BUG, please report', calledby='mixedbasis')


         END DO !l
         IF (fmpi%irank == 0) WRITE (oUnit, '(6X,A,I7)') 'Total:', SUM(mpdata%num_radbasfn(0:hybinp%lcutm1(itype), itype))
      END DO ! itype

      !normalize radbasfn_mt
      call mpdata%normalize(atoms, hybinp, gridf)
      call mpdata%reduce_linear_dep(mpinp,atoms, fmpi, hybinp, gridf, iterHF)

      allocate(basmhlp(atoms%jmtd, maxval(mpdata%num_radbasfn), 0:maxval(hybinp%lcutm1), atoms%ntype))
      basmhlp(1:atoms%jmtd, 1:maxval(mpdata%num_radbasfn), 0:maxval(hybinp%lcutm1), 1:atoms%ntype) &
         = mpdata%radbasfn_mt(1:atoms%jmtd, 1:maxval(mpdata%num_radbasfn), 0:maxval(hybinp%lcutm1), 1:atoms%ntype)
      deallocate(mpdata%radbasfn_mt)
      allocate(mpdata%radbasfn_mt(atoms%jmtd, maxval(mpdata%num_radbasfn), 0:maxval(hybinp%lcutm1), atoms%ntype))
      mpdata%radbasfn_mt = basmhlp

      deallocate(basmhlp, seleco, selecu, selecmat)

      !
      ! now we build linear combinations of the radial functions
      ! such that they possess no moment except one radial function in each l-channel
      !
      IF (fmpi%irank == 0) THEN
         WRITE (oUnit, '(/,A,/,A)') 'Build linear combinations of radial '// &
            'functions in each l-channel,', &
            'such that they possess no multipolmoment'// &
            ' except the last function:'

         WRITE (oUnit, '(/,17x,A)') 'moment  (quality of orthonormality)'
      END IF

      DO itype = 1, atoms%ntype
         n_grid_pt = atoms%jri(itype)

         IF (atoms%ntype > 1 .AND. fmpi%irank == 0) WRITE (oUnit, '(6X,A,I3)') 'Atom type', itype

         DO l = 0, hybinp%lcutm1(itype)
            ! determine radial function with the largest moment
            ! this function is used to build the linear combinations
            max_momentum = 0
            DO i = 1, mpdata%num_radbasfn(l, itype)
               momentum = intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*mpdata%radbasfn_mt(:n_grid_pt, i, l, itype), &
                              atoms, itype, gridf)
               IF (ABS(momentum) > max_momentum) THEN
                  n_radbasfn = i
                  max_momentum = momentum
               END IF
            END DO

            j = 0
            bashlp(:atoms%jri(itype)) = mpdata%radbasfn_mt(:atoms%jri(itype), n_radbasfn, l, itype)
            DO i = 1, mpdata%num_radbasfn(l, itype)
               IF (i == n_radbasfn) CYCLE
               j = j + 1
               mpdata%radbasfn_mt(:atoms%jri(itype), j, l, itype) = mpdata%radbasfn_mt(:atoms%jri(itype), i, l, itype)
            END DO
            mpdata%radbasfn_mt(:atoms%jri(itype), mpdata%num_radbasfn(l, itype), l, itype) = bashlp(:atoms%jri(itype))
         END DO


         DO l = 0, hybinp%lcutm1(itype)
            IF (fmpi%irank == 0) WRITE (oUnit, '(6X,A)') lchar(l)//':'

            IF (mpdata%num_radbasfn(l, itype) == 0) THEN
               IF (fmpi%irank == 0) WRITE (oUnit, '(6X,A,'':   0 ->    '')') lchar(l)
               CYCLE
            END IF

            n_radbasfn = mpdata%num_radbasfn(l, itype)
            DO i = 1, n_radbasfn-1
               ! calculate moment of radial function i
               rdum1 = intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*mpdata%radbasfn_mt(:n_grid_pt, i, l, itype), &
                              atoms, itype, gridf)

               rdum = intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*mpdata%radbasfn_mt(:n_grid_pt, n_radbasfn, l, itype), &
                             atoms, itype, gridf)

               bashlp(:n_grid_pt) = mpdata%radbasfn_mt(:n_grid_pt, n_radbasfn, l, itype)

               IF (SQRT(rdum**2 + rdum1**2) <= 1E-06 .AND. fmpi%irank == 0) &
                  WRITE (oUnit, *) 'Warning: Norm is smaller than 1E-06!'

               ! change function n_radbasfn such that n_radbasfn is orthogonal to i
               ! since the functions radbasfn_mt have been orthogonal on input
               ! the linear combination does not destroy the orthogonality to the residual functions
               mpdata%radbasfn_mt(:n_grid_pt, n_radbasfn, l, itype) = rdum/SQRT(rdum**2 + rdum1**2)*bashlp(:n_grid_pt) &
                                                + rdum1/SQRT(rdum**2 + rdum1**2)*mpdata%radbasfn_mt(:n_grid_pt, i, l, itype)

               ! combine basis function i and n_radbasfn so that they possess no momemt
               mpdata%radbasfn_mt(:n_grid_pt, i, l, itype) = rdum1/SQRT(rdum**2 + rdum1**2)*bashlp(:n_grid_pt) &
                                                - rdum/SQRT(rdum**2 + rdum1**2)*mpdata%radbasfn_mt(:n_grid_pt, i, l, itype)

               rdum1 = intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*mpdata%radbasfn_mt(:n_grid_pt, i, l, itype), &
                              atoms, itype, gridf)

               IF (rdum1 > 1E-10) call judft_error('moment of radial function does not vanish', calledby='mixedbasis')

               IF (fmpi%irank == 0) WRITE (oUnit, '(6x,I4,'' ->  '',ES8.1)') i, rdum1
            END DO
            call mpdata%check_orthonormality(atoms, fmpi, l, itype, gridf)
         ENDDO


      END DO

      call mpdata%check_radbasfn(atoms, hybinp)

      !count basis functions
      hybdat%nbasp = 0
      DO itype = 1, atoms%ntype
         DO i = 1, atoms%neq(itype)
            DO l = 0, hybinp%lcutm1(itype)
               hybdat%nbasp = hybdat%nbasp + (2*l+1) * mpdata%num_radbasfn(l, itype)
            END DO
         END DO
      END DO
      hybdat%nbasm = hybdat%nbasp + mpdata%n_g

      hybdat%maxlmindx = 0
      do itype = 1,atoms%ntype
         hybdat%maxlmindx = max(hybdat%maxlmindx,&
                                SUM([(mpdata%num_radfun_per_l(l, itype)*(2*l + 1), l=0, atoms%lmax(itype))])&
                                )
      enddo
   END SUBROUTINE mixedbasis

   subroutine gen_bas_fun(atoms, enpara, gridf, input, mpdata, fmpi, vr0, usdus, bas1, bas2)
      use m_judft
      use m_types
      USE m_radfun, ONLY: radfun
      USE m_radflo, ONLY: radflo
      USE m_intgrf,   ONLY: intgrf
      implicit NONE
      type(t_atoms), intent(in)        :: atoms
      type(t_enpara), intent(in)       :: enpara
      type(t_input), intent(in)        :: input
      TYPE(t_mpdata), intent(in)      :: mpdata
      type(t_mpi), intent(in)          :: fmpi
      type(t_usdus), intent(inout)     :: usdus

      REAL, ALLOCATABLE, INTENT(INOUT) :: bas1(:,:,:,:,:), bas2(:,:,:,:,:)
      REAL, intent(in)                 :: vr0(:,:,:), gridf(:,:)

      REAL    ::   u(atoms%jmtd, 2, 0:atoms%lmaxd)
      REAL    ::  du(atoms%jmtd, 2, 0:atoms%lmaxd)
      REAL    :: flo(atoms%jmtd, 2, atoms%nlod)

      REAL    :: uuilon(atoms%nlod, atoms%ntype)
      REAL    :: duilon(atoms%nlod, atoms%ntype)
      REAL    :: ulouilopn(atoms%nlod, atoms%nlod, atoms%ntype)
      REAL    :: wronk, norm

      INTEGER :: itype, jspin, i, l, ilo, ok
      INTEGER :: n_grid_pt, noded, nodem
      INTEGER :: l_idx(0:atoms%lmaxd)
      u  = 0.0
      du = 0.0

      ! this is 5-D array. it could cause Problems in bigger systems
      allocate(bas1(atoms%jmtd,    &
                    maxval(mpdata%num_radfun_per_l), &
                    0:atoms%lmaxd, &
                    atoms%ntype,   &
                    input%jspins),   source=0.0, stat=ok)
      if(ok /= 0) call judft_error("Can't allocate bas1 array. Stat= " // int2str(ok))

      allocate(bas2, source=bas1, stat=ok)
      if(ok /= 0) call judft_error("Can't allocate bas1 array. Stat= " // int2str(ok))

      DO itype = 1, atoms%ntype
         n_grid_pt = atoms%jri(itype) ! number of radial gridpoints
         DO jspin = 1, input%jspins
            DO l = 0, atoms%lmax(itype)
               CALL radfun(l, itype, jspin, enpara%el0(l, itype, jspin), vr0(:,itype, jspin), atoms, &
                           u(:,:,l), du(:,:,l), usdus, nodem, noded, wronk)
            END DO
            bas1(1:n_grid_pt, 1, 0:atoms%lmaxd, itype, jspin)  = u(1:n_grid_pt, 1, 0:atoms%lmaxd)
            bas2(1:n_grid_pt, 1, 0:atoms%lmaxd, itype, jspin)  = u(1:n_grid_pt, 2, 0:atoms%lmaxd)
            bas1(1:n_grid_pt, 2, 0:atoms%lmaxd, itype, jspin) = du(1:n_grid_pt, 1, 0:atoms%lmaxd)
            bas2(1:n_grid_pt, 2, 0:atoms%lmaxd, itype, jspin) = du(1:n_grid_pt, 2, 0:atoms%lmaxd)

            ! generate radial functions for local orbitals
            IF (atoms%nlo(itype) >= 1) THEN
               CALL radflo(atoms, itype, jspin, enpara%ello0(1, 1, jspin), vr0(:,itype, jspin), &
                           u, du, fmpi, usdus, uuilon, duilon, ulouilopn, flo)

               l_idx = 2
               DO ilo = 1, atoms%nlo(itype)
                  l = atoms%llo(ilo, itype)
                  l_idx(l) = l_idx(l) + 1
                  bas1(1:n_grid_pt, l_idx(l), atoms%llo(ilo, itype), itype, jspin) = flo(1:n_grid_pt, 1, ilo)
                  bas2(1:n_grid_pt, l_idx(l), atoms%llo(ilo, itype), itype, jspin) = flo(1:n_grid_pt, 2, ilo)
               END DO
            END IF
         END DO
      END DO

      ! the radial functions are normalized
      DO jspin = 1, input%jspins
         DO itype = 1, atoms%ntype
            DO l = 0, atoms%lmax(itype)
               DO i = 1, mpdata%num_radfun_per_l(l, itype)
                  norm = sqrt(intgrf(bas1(:,i, l, itype, jspin)**2 + bas2(:,i, l, itype, jspin)**2, &
                                atoms, itype, gridf))
                  bas1(:atoms%jri(itype), i, l, itype, jspin) = bas1(:atoms%jri(itype), i, l, itype, jspin)/norm
                  bas2(:atoms%jri(itype), i, l, itype, jspin) = bas2(:atoms%jri(itype), i, l, itype, jspin)/norm
               END DO
            END DO
         END DO
      END DO
   end subroutine gen_bas_fun

   function calc_selecmat(atoms,mpdata,seleco, selecu) result(selecmat)
      ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
      use m_types
      use m_judft
      implicit NONE

      type(t_atoms),  intent(in) :: atoms
      TYPE(t_mpdata), intent(in) :: mpdata
      LOGICAL, intent(in) :: seleco(maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd)
      LOGICAL, intent(in) :: selecu(maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd)
      LOGICAL  ::  selecmat(maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, &
                            maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd)
      integer                       :: n1, l1, n2, l2

      ! column-major means left-most index varies the fastest
      do l2=0,atoms%lmaxd
         do n2=1,maxval(mpdata%num_radfun_per_l)
            do l1=0,atoms%lmaxd
               do n1=1,maxval(mpdata%num_radfun_per_l)
                  selecmat(n1,l1,n2,l2) = seleco(n1, l1) .AND. selecu(n2, l2)
               enddo
            enddo
         enddo
      enddo
   end function calc_selecmat
END MODULE m_mixedbasis
