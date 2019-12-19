!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!> This module defines two data-types used in the calculation of xc-potentials
!! a) the abstract t_xcpot which should be overwritten for actual implementations
!! b) the t_gradients that collects the gradients needed in GGAs
!!
!! Currently t_xcpot_inbuild implements the XC-pots directly build into FLEUR
!! and t_xcpot_libxc provides and interface to libxc.
!! In addition to overloading the t_xcpot datatype also mpi_bc_xcpot must be adjusted
!! for additional implementations.
MODULE m_types_xcpot
   USE m_types_potden
   IMPLICIT NONE
   PRIVATE
<<<<<<< variant A
   PUBLIC :: t_xcpot,t_gradients


   TYPE,ABSTRACT,EXTENDS(t_fleurinput_base) :: t_xcpot
>>>>>>> variant B
   PUBLIC           :: t_xcpot,t_gradients

   TYPE t_kinED
      logical             :: set
      real, allocatable   :: is(:,:)   ! (nsp*jmtd, jspins)
      real, allocatable   :: mt(:,:,:) ! (nsp*jmtd, jspins, local num of types)
   contains
      procedure           :: alloc_mt => kED_alloc_mt
   END TYPE t_kinED

   TYPE,ABSTRACT,EXTENDS(t_fleurinput_base) :: t_xcpot
======= end
      REAL :: gmaxxc
      TYPE(t_kinED)    :: kinED
   CONTAINS
      PROCEDURE        :: vxc_is_LDA => xcpot_vxc_is_LDA
      PROCEDURE        :: vxc_is_GGA => xcpot_vxc_is_GGA

      PROCEDURE        :: vx_is_LDA     => xcpot_vx_is_LDA
      PROCEDURE        :: vx_is_GGA     => xcpot_vx_is_GGA
      PROCEDURE        :: vx_is_MetaGGA => xcpot_vx_is_MetaGGA

      PROCEDURE        :: vc_is_LDA => xcpot_vc_is_LDA
      PROCEDURE        :: vc_is_GGA => xcpot_vc_is_GGA

      PROCEDURE        :: exc_is_LDA     => xcpot_exc_is_LDA
      PROCEDURE        :: exc_is_GGA     => xcpot_exc_is_GGA
      PROCEDURE        :: exc_is_MetaGGA => xcpot_exc_is_MetaGGA

      PROCEDURE        :: needs_grad => xcpot_needs_grad
      PROCEDURE        :: is_hybrid  => xcpot_is_hybrid

      PROCEDURE        :: get_exchange_weight => xcpot_get_exchange_weight
      PROCEDURE        :: get_vxc             => xcpot_get_vxc
      PROCEDURE        :: get_exc             => xcpot_get_exc
      PROCEDURE,NOPASS :: alloc_gradients     => xcpot_alloc_gradients
   END TYPE t_xcpot

   TYPE t_gradients
      !Naming convention:
      !t,u,d as last letter for total,up,down
      !agr for absolute value of gradient
      !g2r for laplacien of gradient
      !+??
      REAL,ALLOCATABLE :: agrt(:),agru(:),agrd(:)
      REAL,ALLOCATABLE :: g2ru(:),g2rd(:),gggrt(:)
      REAL,ALLOCATABLE :: gggru(:),gzgr(:),g2rt(:)
      REAL,ALLOCATABLE :: gggrd(:),grgru(:),grgrd(:)
      !These are the contracted Gradients used in libxc
      REAL,ALLOCATABLE :: sigma(:,:)
      REAL,ALLOCATABLE :: vsigma(:,:)
      REAL,ALLOCATABLE :: gr(:,:,:)
      REAL,ALLOCATABLE :: laplace(:,:)
      !For sourcefree calculations, the cartesian components of the gradient
      !on the spherical grid are needed (no spin dependency needed though).
      REAL,ALLOCATABLE :: grxyz(:)
   END TYPE t_gradients
CONTAINS
   SUBROUTINE kED_alloc_mt(kED,nsp_x_jmtd, jspins, n_start, n_types, n_stride)
      IMPLICIT NONE
      class(t_kinED), intent(inout)   :: kED
      integer, intent(in)            :: nsp_x_jmtd, jspins, n_start, n_types, n_stride
      integer                        :: cnt, n

      if(.not. allocated(kED%mt)) then
         cnt = 0
         do n = n_start,n_types,n_stride
            cnt = cnt + 1
         enddo
         allocate(kED%mt(nsp_x_jmtd, jspins, cnt), source=0.0)
      endif
   end subroutine kED_alloc_mt

   ! LDA
   LOGICAL FUNCTION xcpot_vc_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_vc_is_LDA=.false.
   END FUNCTION xcpot_vc_is_LDA

   LOGICAL FUNCTION xcpot_vx_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_vx_is_LDA=.false.
   END FUNCTION xcpot_vx_is_LDA

   LOGICAL FUNCTION xcpot_vxc_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_vxc_is_LDA = xcpot%vx_is_LDA() .and. xcpot%vc_is_LDA()
   END FUNCTION xcpot_vxc_is_LDA

   LOGICAL FUNCTION xcpot_exc_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_exc_is_LDA=.false.
   END FUNCTION xcpot_exc_is_LDA

   ! GGA
   LOGICAL FUNCTION xcpot_vc_is_GGA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_vc_is_GGA=.false.
   END FUNCTION xcpot_vc_is_GGA

   LOGICAL FUNCTION xcpot_vx_is_GGA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_vx_is_GGA=.false.
   END FUNCTION xcpot_vx_is_GGA

   LOGICAL FUNCTION xcpot_vxc_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_vxc_is_gga= xcpot%vx_is_GGA() .and. xcpot%vc_is_GGA()
   END FUNCTION xcpot_vxc_is_gga

   LOGICAL FUNCTION xcpot_exc_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_exc_is_gga=.false.
   END FUNCTION xcpot_exc_is_gga

   LOGICAL FUNCTION xcpot_vx_is_MetaGGA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_vx_is_MetaGGA=.false.
   END FUNCTION xcpot_vx_is_MetaGGA

   LOGICAL FUNCTION xcpot_exc_is_MetaGGA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_exc_is_MetaGGA=.false.
   END FUNCTION xcpot_exc_is_MetaGGA

   LOGICAL FUNCTION xcpot_needs_grad(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot

      xcpot_needs_grad= xcpot%vc_is_gga() .or. xcpot%vx_is_MetaGGA()
   END FUNCTION xcpot_needs_grad

   LOGICAL FUNCTION xcpot_is_hybrid(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      xcpot_is_hybrid=.FALSE.
   END FUNCTION xcpot_is_hybrid

   FUNCTION xcpot_get_exchange_weight(xcpot) RESULT(a_ex)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN):: xcpot
      REAL:: a_ex
      a_ex=-1
   END FUNCTION xcpot_get_exchange_weight

   SUBROUTINE xcpot_get_vxc(xcpot,jspins,rh,vxc,vx,grad, kinED_KS)
      USE m_judft
      IMPLICIT NONE

      CLASS(t_xcpot),INTENT(IN) :: xcpot
      INTEGER, INTENT (IN)     :: jspins
      !--> charge density
      REAL,INTENT (IN)         :: rh(:,:)
      !---> xc potential
      REAL, INTENT (OUT)       :: vxc (:,:),vx(:,:)
      TYPE(t_gradients),OPTIONAL,INTENT(INOUT)::grad
      REAL, INTENT(IN),OPTIONAL:: kinED_KS(:,:)

      vxc = 0.0
      vx  = 0.0
      call juDFT_error("Can't use XC-parrent class")
   END SUBROUTINE xcpot_get_vxc

   SUBROUTINE xcpot_get_exc(xcpot,jspins,rh,exc,grad,kinED_KS, mt_call)
      USE m_types_misc
      USE m_judft
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE

      CLASS(t_xcpot),INTENT(IN)             :: xcpot
      INTEGER, INTENT (IN)                  :: jspins
      !--> charge density
      REAL,INTENT (IN)                      :: rh(:,:)
      !--> kinetic energy density
      !---> xc energy density
      REAL, INTENT (OUT)                    :: exc (:)
      TYPE(t_gradients),OPTIONAL,INTENT(IN) :: grad
      LOGICAL, OPTIONAL, INTENT(IN)         :: mt_call
      REAL, INTENT(IN), OPTIONAL            :: kinED_KS(:,:)

      exc = 0.0
      call juDFT_error("Can't use XC-parrent class")
   END SUBROUTINE xcpot_get_exc

   SUBROUTINE xcpot_alloc_gradients(ngrid,jspins,grad)
      IMPLICIT NONE

      INTEGER, INTENT (IN)         :: jspins,ngrid
      TYPE(t_gradients),INTENT(INOUT):: grad

      IF (allocated(grad%agrt)) THEN
         DEALLOCATE(grad%agrt,grad%agru,grad%agrd)
         DEALLOCATE(grad%g2ru,grad%g2rd,grad%gggrt)
         DEALLOCATE(grad%gggru,grad%gzgr,grad%g2rt)
         DEALLOCATE(grad%gggrd,grad%grgru,grad%grgrd)
      ENDIF
      !For the in-build xc-pots
      ALLOCATE(grad%agrt(ngrid),grad%agru(ngrid),grad%agrd(ngrid))
      ALLOCATE(grad%g2ru(ngrid),grad%g2rd(ngrid),grad%gggrt(ngrid))
      ALLOCATE(grad%gggru(ngrid),grad%gzgr(ngrid),grad%g2rt(ngrid))
      ALLOCATE(grad%gggrd(ngrid),grad%grgru(ngrid),grad%grgrd(ngrid))
  END SUBROUTINE xcpot_alloc_gradients

END MODULE m_types_xcpot