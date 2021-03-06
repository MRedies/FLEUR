!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   This module generates the little group of k and the extended irr. !
!   BZ. Furthermore it calculates the irr. representation             !
!                                                                     !
!   P(R,T)\phi_n,k = \sum_{n'} rep_v(n',n) *\phi_n',k        !
!   where                                                             !
!         P  is an element of the little group of k                   !
!         n' runs over group of degenerat states belonging to n.      !
!                                                                     !
!                                             M.Betzinger (09/07)     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE m_symm_hf

   use m_judft
   USE m_types
   USE m_types_hybdat
   USE m_constants
   USE m_util
   USE m_intgrf
   USE m_io_hybinp

CONTAINS

   SUBROUTINE symm_hf_init(fi, nk, nsymop, rrot, psym)
      use m_juDFT
      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      INTEGER, INTENT(IN)    :: nk
      INTEGER, INTENT(INOUT) :: nsymop
      INTEGER, INTENT(INOUT) :: rrot(:, :, :) ! 3,3,fi%sym%nsym
      INTEGER, INTENT(INOUT) :: psym(:) ! Note: psym is only filled up to index nsymop

      INTEGER :: i
      REAL    :: rotkpt(3)

      CALL timestart("symm_hf_init")
      nsymop = 0
      ! calculate rotations in reciprocal space
      DO i = 1, fi%sym%nsym
         IF(i <= fi%sym%nop) THEN
            rrot(:, :, i) = transpose(fi%sym%mrot(:, :, fi%sym%invtab(i)))
         ELSE
            rrot(:, :, i) = -rrot(:, :, i - fi%sym%nop)
         END IF
      END DO

      ! determine little group of k., i.e. those symmetry operations
      ! which keep bk(:,nk) invariant
      ! nsymop :: number of such symmetry-operations
      ! psym   :: points to the symmetry-operation

      psym = 0
      nsymop = 0
      DO i = 1, fi%sym%nsym
         rotkpt = matmul(rrot(:, :, i), fi%kpts%bkf(:, nk))

         !transfer rotkpt into BZ
         rotkpt = fi%kpts%to_first_bz(rotkpt)

         !check if rotkpt is identical to bk(:,nk)
         IF(maxval(abs(rotkpt - fi%kpts%bkf(:, nk))) <= 1E-07) THEN
            nsymop = nsymop + 1
            psym(nsymop) = i
         END IF
      END DO

      WRITE(oUnit, '(A,i3)') 'current nk: ', nk
      WRITE(oUnit, '(A,3f10.5)') ' fi%kpts%bkf(:,nk):', fi%kpts%bkf(:, nk)
      WRITE(oUnit, '(A,i3)') ' Number of elements in the little group:', nsymop

      CALL timestop("symm_hf_init")
   END SUBROUTINE symm_hf_init

   SUBROUTINE symm_hf(fi, nk, hybdat, eig_irr, mpdata, lapw, nococonv, zmat, c_phase, jsp, &
                      rrot, nsymop, psym, n_q, parent, nsest, indx_sest)

      USE m_olap
      USE m_trafo
      use m_calc_cmt
      use m_juDFT

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_hybdat), INTENT(IN) :: hybdat

      TYPE(t_mpdata), intent(in) :: mpdata
      TYPE(t_lapw), INTENT(IN)   :: lapw
      type(t_nococonv), intent(in):: nococonv
      type(t_mat), intent(in)    :: zmat

!     - scalars -
      INTEGER, INTENT(IN)              :: nk
      INTEGER, INTENT(IN)              :: jsp
      INTEGER, INTENT(IN)              :: nsymop

!     - arrays -
      complex, intent(in)              :: c_phase(hybdat%nbands(nk))
      INTEGER, INTENT(IN)              :: rrot(:, :, :)
      INTEGER, INTENT(IN)              :: psym(:)
      INTEGER, INTENT(INOUT)           :: parent(fi%kpts%nkptf)
      INTEGER, INTENT(INOUT)           :: nsest(hybdat%nbands(nk))
      INTEGER, INTENT(INOUT)           :: indx_sest(hybdat%nbands(nk), hybdat%nbands(nk))
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: n_q(:)

      REAL, INTENT(IN)                 :: eig_irr(:, :)

!     - local scalars -
      INTEGER                         :: ikpt, ikpt1, iop, isym, iisym, m
      INTEGER                         :: itype, ieq, iatom, ratom
      INTEGER                         :: iband, iband1, iband2, iatom0
      INTEGER                         :: i, j, ic, ic1, ic2
      INTEGER                         :: ok
      INTEGER                         :: l, lm
      INTEGER                         :: n1, n2, nn
      INTEGER                         :: ndb, ndb1, ndb2
      INTEGER                         :: nrkpt
      INTEGER                         :: maxndb, nddb

      REAL                            :: tolerance, pi

      COMPLEX                         :: cdum
      COMPLEX, PARAMETER             :: img = (0.0, 1.0)

!     - local arrays -
      INTEGER                         :: neqvkpt(fi%kpts%nkptf)
      INTEGER                         :: list(fi%kpts%nkptf)
      INTEGER                         :: degenerat(hybdat%ne_eig(nk))

      REAL                            :: rotkpt(3), g(3)
      REAL, ALLOCATABLE               :: olapmt(:, :, :, :)

      COMPLEX                         :: cmt(hybdat%nbands(nk), hybdat%maxlmindx, fi%atoms%nat)
      COMPLEX                         :: carr1(hybdat%maxlmindx, fi%atoms%nat)
      COMPLEX, ALLOCATABLE             :: carr(:), wavefolap(:, :)
      COMPLEX, ALLOCATABLE             :: cmthlp(:, :, :)
      COMPLEX, ALLOCATABLE             :: cpwhlp(:, :)
      COMPLEX, ALLOCATABLE             :: trace(:, :)

      TYPE(t_mat)                      :: olappw, z
      COMPLEX, ALLOCATABLE             :: rep_d(:, :, :)
      LOGICAL, ALLOCATABLE             :: symequivalent(:, :)

      CALL timestart("symm_hf")
      parent = 0; nsest = 0; indx_sest = 0;
      WRITE(oUnit, '(A)') new_line('n')//new_line('n')//'### subroutine: symm ###'

      ! determine extented irreducible BZ of k ( EIBZ(k) ), i.e.
      ! those k-points, which can generate the whole BZ by
      ! applying the symmetry operations of the little group of k
      call timestart("calc EIBZ")
      neqvkpt = 0

      DO i = 1, fi%kpts%nkptf
         list(i) = i - 1
      END DO

      DO ikpt = 2, fi%kpts%nkptf
         DO iop = 1, nsymop

            rotkpt = matmul(rrot(:, :, psym(iop)), fi%kpts%bkf(:, ikpt))

            !determine number of rotkpt
            nrkpt = fi%kpts%get_nk(rotkpt)
            IF(nrkpt == 0) call judft_error('symm: Difference vector not found !')

            IF(list(nrkpt) /= 0) THEN
               list(nrkpt) = 0
               neqvkpt(ikpt) = neqvkpt(ikpt) + 1
               parent(nrkpt) = ikpt
            END IF
            IF(all(list == 0)) EXIT

         END DO
      END DO

      ! for the Gamma-point holds:
      parent(1) = 1
      neqvkpt(1) = 1

      WRITE(oUnit, '(A,i5)') ' Number of k-points in the EIBZ', fi%kpts%EIBZ(nk)%nkpt
      call timestop("calc EIBZ")

      ! determine the factor n_q, that means the number of symmetrie operations of the little group of bk(:,nk)
      ! which keep q (in EIBZ) invariant
      call timestart("calc n_q")
      IF(ALLOCATED(n_q)) DEALLOCATE(n_q)
      allocate(n_q(fi%kpts%EIBZ(nk)%nkpt), source=0)

      ic = 0
      n_q = 0
      DO ikpt = 1, fi%kpts%nkptf
         IF(parent(ikpt) == ikpt) THEN
            ic = ic + 1
            DO iop = 1, nsymop
               isym = psym(iop)
               rotkpt = matmul(rrot(:, :, isym), fi%kpts%bkf(:, ikpt))

               !transfer rotkpt into BZ
               rotkpt = fi%kpts%to_first_bz(rotkpt)

               !check if rotkpt is identical to bk(:,ikpt)
               IF(maxval(abs(rotkpt - fi%kpts%bkf(:, ikpt))) <= 1E-06) THEN
                  n_q(ic) = n_q(ic) + 1
               END IF
            END DO
         END IF
      END DO
      IF(ic /= fi%kpts%EIBZ(nk)%nkpt) call judft_error('symm: failure EIBZ')
      call timestop("calc n_q")

      ! calculate degeneracy:
      ! degenerat(i) = 1 state i  is not degenerat,
      ! degenerat(i) = j state i has j-1 degenerat states at {i+1,...,i+j-1}
      ! degenerat(i) = 0 state i is degenerat

      tolerance = 1E-07 !0.00001

      degenerat = 1

      WRITE(oUnit, '(A,f10.8)') ' Tolerance for determining degenerate states=', tolerance

      DO i = 1, hybdat%nbands(nk)
         DO j = i + 1, hybdat%nbands(nk)
            IF(abs(eig_irr(i, nk) - eig_irr(j, nk)) <= tolerance) THEN
               degenerat(i) = degenerat(i) + 1
            END IF
         END DO
      END DO

      DO i = 1, hybdat%ne_eig(nk)
         IF(degenerat(i) /= 1 .or. degenerat(i) /= 0) THEN
            degenerat(i + 1:i + degenerat(i) - 1) = 0
         END IF
      END DO

      ! maximal number of degenerate bands -> maxndb
      maxndb = maxval(degenerat)

      ! number of different degenerate bands/states
      nddb = count(degenerat >= 1)
#ifdef CPP_EXPLICIT_HYB
      WRITE(oUnit, *) ' Degenerate states:'
      DO iband = 1, hybdat%nbands(nk)/5 + 1
         WRITE(oUnit, '(5i5)') degenerat(iband*5 - 4:min(iband*5, hybdat%nbands(nk)))
      END DO
#endif
      ! read in cmt and z at current k-point (nk)
      call calc_cmt(fi%atoms, fi%cell, fi%input, fi%noco,nococonv, fi%hybinp, hybdat, mpdata, fi%kpts, &
                          fi%sym, fi%oned, zmat, jsp, nk, c_phase, cmt)

      IF(allocated(olapmt)) deallocate(olapmt)
      allocate(olapmt(maxval(mpdata%num_radfun_per_l), maxval(mpdata%num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), stat=ok)
      IF(ok /= 0) call judft_error('symm: failure allocation olapmt')
      olapmt = 0

      call timestart("calc olapmt")
      DO itype = 1, fi%atoms%ntype
         DO l = 0, fi%atoms%lmax(itype)
            nn = mpdata%num_radfun_per_l(l, itype)
            DO n2 = 1, nn
               DO n1 = 1, nn
                  olapmt(n1, n2, l, itype) = intgrf( &
                                             hybdat%bas1(:, n1, l, itype)*hybdat%bas1(:, n2, l, itype) &
                                             + hybdat%bas2(:, n1, l, itype)*hybdat%bas2(:, n2, l, itype), &
                                             fi%atoms, itype, hybdat%gridf)
               END DO
            END DO
         END DO
      END DO
      call timestop("calc olapmt")

      allocate(wavefolap(hybdat%nbands(nk), hybdat%nbands(nk)), carr(maxval(mpdata%num_radfun_per_l)), stat=ok)
      IF(ok /= 0) call judft_error('symm: failure allocation wfolap/maxindx')
      wavefolap = 0

      call timestart("calc wavefolap")
      iatom = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1
            lm = 0
            DO l = 0, fi%atoms%lmax(itype)
               DO M = -l, l
                  nn = mpdata%num_radfun_per_l(l, itype)
                  DO iband1 = 1, hybdat%nbands(nk)
                     carr(:nn) = matmul(olapmt(:nn, :nn, l, itype), &
                                        cmt(iband1, lm + 1:lm + nn, iatom))
                     DO iband2 = 1, iband1
                        wavefolap(iband2, iband1) &
                           = wavefolap(iband2, iband1) &
                             + dot_product(cmt(iband2, lm + 1:lm + nn, iatom), carr(:nn))
                     END DO
                  END DO
                  lm = lm + nn
               END DO
            END DO
         END DO
      END DO
      call timestop("calc wavefolap")

      DO iband1 = 1, hybdat%nbands(nk)
         DO iband2 = 1, iband1
            wavefolap(iband1, iband2) = conjg(wavefolap(iband2, iband1))
         END DO
      END DO

      allocate(symequivalent(nddb, nddb), stat=ok)
      IF(ok /= 0) call judft_error('symm: failure allocation symequivalent')
      symequivalent = .false.
      ic1 = 0
      DO iband1 = 1, hybdat%nbands(nk)
         ndb1 = degenerat(iband1)
         IF(ndb1 == 0) CYCLE
         ic1 = ic1 + 1
         ic2 = 0
         DO iband2 = 1, hybdat%nbands(nk)
            ndb2 = degenerat(iband2)
            IF(ndb2 == 0) CYCLE
            ic2 = ic2 + 1
            IF(any(abs(wavefolap(iband1:iband1 + ndb1 - 1, &
                                 iband2:iband2 + ndb2 - 1)) > 1E-9)) THEN
!                .and. ndb1 .eq. ndb2 ) THEN
               symequivalent(ic2, ic1) = .true.
            END IF
         END DO
      END DO

      !
      ! generate index field which contain the band combinations (n1,n2),
      ! which are non zero
      !
      call timestart("calc bandcombos")
      ic1 = 0
      indx_sest = 0
      nsest = 0
      DO iband1 = 1, hybdat%nbands(nk)
         ndb1 = degenerat(iband1)
         IF(ndb1 >= 1) ic1 = ic1 + 1
         i = 0
         DO WHILE(degenerat(iband1 - i) == 0)
            i = i + 1
         END DO
         ndb1 = degenerat(iband1 - i)
         ic2 = 0
         DO iband2 = 1, hybdat%nbands(nk)
            ndb2 = degenerat(iband2)
            IF(ndb2 >= 1) ic2 = ic2 + 1
            i = 0
            DO WHILE(degenerat(iband2 - i) == 0)
               i = i + 1
            END DO
            ndb2 = degenerat(iband2 - i)
            ! only upper triangular part
            IF(symequivalent(ic2, ic1) .and. iband2 <= iband1) THEN
!            IF( ndb1 .ne. ndb2 ) call judft_error('symm_hf: failure symequivalent')
               nsest(iband1) = nsest(iband1) + 1
               indx_sest(nsest(iband1), iband1) = iband2
            END IF
         END DO
      END DO
      call timestop("calc bandcombos")

      !
      ! calculate representations for core states
      ! (note that for a core state, these are proportional to the Wigner D matrices)
      !
      ! Definition of the Wigner rotation matrices
      !
      !                     -1                l
      ! P(R) Y  (r)  = Y  (R  r)  =  sum     D    (R)  Y   (r)
      !       lm        lm              m'    m m'      lm'
      !

      pi = pimach()

      call timestart("calc core repr")
      IF(hybdat%lmaxcd > fi%atoms%lmaxd) then
         call judft_error('symm_hf: The very impropable case that hybdat%lmaxcd > fi%atoms%lmaxd occurs')
      endif
      iatom = 0
      iatom0 = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1
            DO iop = 1, nsymop
               isym = psym(iop)
               IF(isym <= fi%sym%nop) THEN
                  iisym = isym
               ELSE
                  iisym = isym - fi%sym%nop
               END IF

               ratom = fi%hybinp%map(iatom, isym)
               rotkpt = matmul(rrot(:, :, isym), fi%kpts%bkf(:, nk))
               g = nint(rotkpt - fi%kpts%bkf(:, nk))

               cdum = exp(-2*pi*img*dot_product(rotkpt, fi%sym%tau(:, iisym)))* &
                      exp(2*pi*img*dot_product(g, fi%atoms%taual(:, ratom)))
            END DO
         END DO
         iatom0 = iatom0 + fi%atoms%neq(itype)
      END DO
      call timestop("calc core repr")

      CALL timestop("symm_hf")
   END SUBROUTINE symm_hf
END MODULE m_symm_hf
