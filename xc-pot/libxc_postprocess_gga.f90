!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_libxc_postprocess_gga
CONTAINS

   SUBROUTINE libxc_postprocess_gga_mt(xcpot,atoms,sphhar,n,grad,rho,v_xc)
      USE m_mt_tofrom_grid
      USE m_types
      use m_judft_string

      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN)   :: xcpot
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      INTEGER,INTENT(IN)          :: n
      TYPE(t_gradients),INTENT(IN):: grad
      REAL,INTENT(IN)             :: rho(:,:)
      REAL,INTENT(INOUT)          :: v_xc(:,:)

      INTEGER :: nsp,n_sigma,i
      REAL,ALLOCATABLE:: vsigma(:,:),vsigma_mt(:,:,:)
      TYPE(t_gradients)::grad_vsigma

      n_sigma=MERGE(1,3,SIZE(v_xc,2)==1) !Number of contracted gradients in libxc 1 for non-spin-polarized, 3 otherwise
      nsp=SIZE(v_xc,1) !no of points
      ALLOCATE(vsigma(nsp,n_sigma),vsigma_mt(atoms%jri(n),0:sphhar%nlhd,n_sigma))
      vsigma_mt=0.0
      vsigma=TRANSPOSE(grad%vsigma) !create a (nsp,n_sigma) matrix
      CALL mt_from_grid(atoms,sphhar,n,n_sigma,vsigma,vsigma_mt)
      DO i=1,atoms%jri(n)
         vsigma_mt(i,:,:)=vsigma_mt(i,:,:)*atoms%rmsh(i,n)**2
      ENDDO
      ALLOCATE(grad_vsigma%gr(3,nsp,n_sigma))
      CALL mt_to_grid(xcpot%needs_grad(),n_sigma,atoms,sphhar,vsigma_mt,n,grad=grad_vsigma)

      CALL libxc_postprocess_gga(transpose(grad%vsigma),grad,grad_vsigma,rho,v_xc)
   END SUBROUTINE libxc_postprocess_gga_mt

   SUBROUTINE libxc_postprocess_gga_pw(xcpot,stars,cell,grad,rho,v_xc)
      USE m_pw_tofrom_grid
      USE m_types

      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN)   :: xcpot
      TYPE(t_stars),INTENT(IN)    :: stars
      TYPE(t_cell),INTENT(IN)     :: cell
      TYPE(t_gradients),INTENT(IN):: grad
      REAL, INTENT(IN)            :: rho(:,:)
      REAL,INTENT(INOUT)          :: v_xc(:,:)

      COMPLEX,ALLOCATABLE:: vsigma_g(:,:)
      REAL,ALLOCATABLE:: vsigma(:,:)
      TYPE(t_gradients)::grad_vsigma
      INTEGER :: nsp,n_sigma, mini_loc
      INTEGER, allocatable :: mini_array(:)
      nsp=SIZE(v_xc,1) !no of points
      n_sigma=MERGE(1,3,SIZE(v_xc,2)==1) !See in _mt routine
      ALLOCATE(vsigma_g(stars%ng3,n_sigma),vsigma(nsp,n_sigma)); vsigma_g=0.0
      vsigma=TRANSPOSE(grad%vsigma) !create a (nsp,n_sigma) matrix
      CALL pw_from_grid(xcpot%needs_grad(),stars,.FALSE.,vsigma,vsigma_g)
      !vsigma_g(:,1)=vsigma_g(:,1)*stars%nstr(:)
      ALLOCATE(grad_vsigma%gr(3,nsp,n_sigma))
      CALL pw_to_grid(xcpot%needs_grad(),n_sigma,.false.,stars,cell,vsigma_g,grad_vsigma,xcpot)

      mini_array = minloc(rho)
      mini_loc = mini_array(1)
      write (*,*) "minloc rho: ", mini_loc
      write (*,*) "rho = ", rho(mini_loc,1)
      write (*,*) "grad%gr", grad%gr(:,mini_loc,1)
      write (*,*) "min |grad| = ", minval(abs(grad%gr))
      write (*,*) "max |grad| = ", maxval(abs(grad%gr)), new_line("a")


      write (*,*) "lapl", grad%laplace(mini_loc,1)
      write (*,*) "min |lapl| = ", minval(abs(grad%laplace))
      write (*,*) "max |lapl| = ", maxval(abs(grad%laplace)), new_line("a")

      write (*,*) "vsigma = ", vsigma(mini_loc,1)
      write (*,*) "min |vsigma| = ", minval(abs(vsigma))
      write (*,*) "max |vsigma| = ", maxval(abs(vsigma)), new_line("a")

      write (*,*) "grad_vsigma = ", grad_vsigma%gr(:,mini_loc,1)
      write (*,*) "min |grad_vsigma| = ", minval(abs(grad_vsigma%gr))
      write (*,*) "max |grad_vsigma| = ", maxval(abs(grad_vsigma%gr)), new_line("a")

      write (*,*) "v_xc", v_xc(mini_loc,1)
      write (*,*) "grad vsig . grad den", &
                dot_PRODUCT(grad_vsigma%gr(:,mini_loc,1), grad%gr(:,mini_loc,1))
      write (*,*) "vsig * lapl", -2*vsigma(mini_loc,1)*grad%laplace(mini_loc,1), new_line("a")

      CALL libxc_postprocess_gga(transpose(grad%vsigma),grad,grad_vsigma,rho,v_xc)
   END SUBROUTINE libxc_postprocess_gga_pw

   SUBROUTINE libxc_postprocess_gga(vsigma,grad,grad_vsigma,rho,v_xc)
      USE m_types
      IMPLICIT NONE
      REAL,INTENT(IN)             :: vsigma(:,:)
      TYPE(t_gradients),INTENT(IN):: grad,grad_vsigma
      REAL,INTENT(IN)             :: rho(:,:)
      REAL,INTENT(INOUT)          :: v_xc(:,:)
      INTEGER:: i

      REAL, PARAMETER             :: den_cutoff = 1e-3


      IF (SIZE(v_xc,2)==1) THEN !Single spin
         DO i=1,SIZE(v_xc,1) !loop over points
            if(rho(i,1) > den_cutoff) then
               v_xc(i,1)=v_xc(i,1)-2*dot_PRODUCT(grad_vsigma%gr(:,i,1),&
                                                 grad%gr(:,i,1))-2*vsigma(i,1)*grad%laplace(i,1)
            endif
         ENDDO
      ELSE  !two spins
         DO i=1,SIZE(v_xc,1) !loop over points
            if(rho(i,1) > den_cutoff) then
               v_xc(i,1)=v_xc(i,1)-2*dot_PRODUCT(grad_vsigma%gr(:,i,1),grad%gr(:,i,1))-2*vsigma(i,1)*grad%laplace(i,1)-&
                          dot_PRODUCT(grad_vsigma%gr(:,i,2),grad%gr(:,i,2))-vsigma(i,2)*grad%laplace(i,2)
            ENDIF
            if(rho(i,2) > den_cutoff) then
               v_xc(i,2)=v_xc(i,2)-2*dot_PRODUCT(grad_vsigma%gr(:,i,3),grad%gr(:,i,2))-2*vsigma(i,3)*grad%laplace(i,2)-&
                          dot_PRODUCT(grad_vsigma%gr(:,i,2),grad%gr(:,i,1))-vsigma(i,2)*grad%laplace(i,1)
            ENDIF
         ENDDO
      END IF

   END SUBROUTINE libxc_postprocess_gga

END MODULE m_libxc_postprocess_gga
