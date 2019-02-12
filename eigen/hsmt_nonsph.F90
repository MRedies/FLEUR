!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_nonsph
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC hsmt_nonsph
  INTERFACE hsmt_nonsph_noMPI
    module procedure hsmt_nonsph_noMPI_cpu
#ifdef CPP_GPU
    module procedure hsmt_nonsph_noMPI_gpu
#endif
  END INTERFACE
CONTAINS
  SUBROUTINE hsmt_nonsph(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
    USE m_types
    USE m_hsmt_ab
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_tlmplm),INTENT(IN)     :: td
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)          :: n,isp,iintsp,jintsp
    COMPLEX,INTENT(IN)            :: chi
    !     .. Array Arguments ..
#if defined CPP_GPU
    REAL,MANAGED,INTENT(IN)    :: fj(:,:,:),gj(:,:,:)
#else
    REAL,INTENT(IN)            :: fj(:,0:,:),gj(:,0:,:)
#endif
    CLASS(t_mat),INTENT(INOUT)     ::hmat

    !     .. Local Arrays
#if defined CPP_GPU
    COMPLEX,ALLOCATABLE,DEVICE  :: h_loc_dev(:,:)
    COMPLEX,ALLOCATABLE,MANAGED :: ab(:,:,:)
#else
    COMPLEX,ALLOCATABLE:: ab(:,:,:)
#endif

    !     .. Local Scalars
    INTEGER  nn, na, i, ab_size

    CALL timestart("non-spherical setup")

    IF (hmat%l_real) THEN
       IF (ANY(SHAPE(hmat%data_c)/=SHAPE(hmat%data_r))) THEN
          DEALLOCATE(hmat%data_c)
          ALLOCATE(hmat%data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       ENDIF
       hmat%data_c=0.0
    ENDIF

    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2,MIN(jintsp,iintsp):MAX(jintsp,iintsp)))

    DO nn = 1,atoms%neq(n)
       na = SUM(atoms%neq(:n-1))+nn
       IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN

          CALL timestart("hsmt_abNSPH")
          DO i = MIN(jintsp,iintsp),MAX(jintsp,iintsp)
              CALL hsmt_ab(mpi,sym,atoms,noco,isp,i,n,na,cell,lapw,fj,gj,ab(:,:,i),ab_size,.TRUE.) 
          ENDDO
          CALL timestop("hsmt_abNSPH")

          IF (mpi%n_size==1) THEN
#if defined CPP_GPU
             ALLOCATE(h_loc_dev(size(td%h_loc,1),size(td%h_loc,2)))
             h_loc_dev(1:,1:) = CONJG(td%h_loc(0:,0:,n,isp)) 

             CALL hsmt_nonsph_noMPI(n,na,mpi,atoms,isp,iintsp,jintsp,chi,lapw,h_loc_dev,ab,ab_size,hmat)
#else
             CALL hsmt_nonsph_noMPI(n,na,mpi,atoms,isp,iintsp,jintsp,chi,lapw,td,ab,ab_size,hmat)
#endif
          ELSE
             CALL hsmt_nonsph_MPI(n,na,mpi,atoms,isp,iintsp,jintsp,chi,lapw,td,ab,ab_size,hmat)
          ENDIF

       ENDIF
    ENDDO

    IF (hmat%l_real) THEN
       hmat%data_r=hmat%data_r+REAL(hmat%data_c)
    ENDIF

    CALL timestop("non-spherical setup")
  END SUBROUTINE hsmt_nonsph

  SUBROUTINE hsmt_nonsph_MPI(n,na,mpi,atoms,isp,iintsp,jintsp,chi,lapw,td,ab,ab_size,hmat)
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),   INTENT(IN)     :: mpi
    TYPE(t_atoms), INTENT(IN)     :: atoms
    TYPE(t_lapw),  INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)     :: td
    INTEGER,       INTENT(IN)     :: n,na,isp,iintsp,jintsp,ab_size
    COMPLEX,       INTENT(IN)     :: chi,ab(:,:,:)
    CLASS(t_mat),  INTENT(INOUT)  :: hmat
    !     ..
    !     .. Local Variables
    COMPLEX,ALLOCATABLE :: ab1(:,:),ab_select(:,:)
    INTEGER  l,ll,m
    REAL     rchi

    ALLOCATE(ab1(lapw%nv(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))
    ALLOCATE(ab_select(lapw%num_local_cols(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))
    rchi=MERGE(REAL(chi),REAL(chi)*2,(atoms%invsat(na)==0))
    IF (iintsp==jintsp) THEN
       !Even if iintsp=jintsp=2 in the calling subroutine (hsmt_nonsph),
       !the local count starts with 1
       CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab(1,1,1),SIZE(ab,1),&
                   td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
       !Cut out of ab1 only the needed elements here
       ab_select=ab1(mpi%n_rank+1:lapw%nv(jintsp):mpi%n_size,:)
       CALL zgemm("N","T",lapw%nv(iintsp),lapw%num_local_cols(iintsp),ab_size,CMPLX(rchi,0.0),CONJG(ab1),&
                   SIZE(ab1,1),ab_select,lapw%num_local_cols(iintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))
    ELSE
       CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab(1,1,jintsp),SIZE(ab,1),&
                   td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
       !Cut out of ab1 only the needed elements here
       ab_select=ab1(mpi%n_rank+1:lapw%nv(jintsp):mpi%n_size,:)
       CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab(1,1,iintsp),SIZE(ab,1),td%h_loc(:,:,n,isp),&
                   SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
       CALL zgemm("N","t",lapw%nv(iintsp),lapw%num_local_cols(jintsp),ab_size,chi,conjg(ab1),SIZE(ab1,1),&
                   ab_select,lapw%num_local_cols(jintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))   
    ENDIF

  END SUBROUTINE hsmt_nonsph_MPI

  SUBROUTINE hsmt_nonsph_noMPI_cpu(n,na,mpi,atoms,isp,iintsp,jintsp,chi,lapw,td,ab,ab_size,hmat)
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),   INTENT(IN)    :: mpi
    TYPE(t_atoms), INTENT(IN)    :: atoms
    TYPE(t_lapw),  INTENT(IN)    :: lapw
    TYPE(t_tlmplm),INTENT(IN)    :: td
    INTEGER,       INTENT(IN)    :: n,na,isp,iintsp,jintsp,ab_size
    COMPLEX,       INTENT(in)    :: chi,ab(:,:,:)
    CLASS(t_mat),  INTENT(INOUT) :: hmat
    !     ..
    !     .. Local Variables
    COMPLEX, ALLOCATABLE :: ab1(:,:),ab2(:,:)
    REAL     rchi
    INTEGER  l,ll,m

    ALLOCATE(ab1(MAXVAL(lapw%nv),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))
    IF (iintsp.NE.jintsp) ALLOCATE(ab2(MAXVAL(lapw%nv),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))
    rchi=MERGE(REAL(chi),REAL(chi)*2,(atoms%invsat(na)==0))
    IF (iintsp==jintsp) THEN
       !Even if iintsp=jintsp=2 in the calling subroutine (hsmt_nonsph),
       !the local count starts with 1
       CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab(1,1,1),SIZE(ab,1),&
                   td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
       IF (isp<3) THEN
          CALL zherk("U","N",lapw%nv(iintsp),ab_size,Rchi,CONJG(ab1),SIZE(ab1,1),1.0,hmat%data_c,SIZE(hmat%data_c,1))
       ELSE !This is the case of a local off-diagonal contribution.
            !It is not Hermitian, so we need to USE zgemm CALL
          CALL zgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,CONJG(ab(1,1,jintsp)),SIZE(ab,1),&
                      ab1,SIZE(ab1,1),CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
       ENDIF
    ELSE  !here the l_ss off-diagonal part starts
       CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab(1,1,jintsp),SIZE(ab,1),&
                   td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
       CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab(1,1,iintsp),SIZE(ab,1),&
                   td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab2,SIZE(ab2,1))
       CALL zgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,conjg(ab2),SIZE(ab2,1),&
                   ab1,SIZE(ab1,1),CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
    ENDIF
    
 END SUBROUTINE hsmt_nonsph_noMPI_cpu

#if defined CPP_GPU
  SUBROUTINE hsmt_nonsph_noMPI_gpu(n,na,mpi,atoms,isp,iintsp,jintsp,chi,lapw,h_loc_dev,ab_dev,ab_size,hmat)
  !note that basically all matrices in the GPU version are conjugates of their cpu counterparts
    USE m_types
  !   cublas: required to use generic BLAS interface
  !   cudafor: required to use CUDA runtime API routines 
  !   nvtx: profiling
    USE cublas   
    USE cudafor
    USE nvtx

    IMPLICIT NONE
    TYPE(t_mpi),   INTENT(IN)         :: mpi
    TYPE(t_atoms), INTENT(IN)         :: atoms
    TYPE(t_lapw),  INTENT(IN)         :: lapw
    INTEGER,       INTENT(IN)         :: n,na,isp,iintsp,jintsp,ab_size
    COMPLEX,       INTENT(IN)         :: chi
    COMPLEX,       INTENT(IN), DEVICE :: h_loc_dev(:,:)
    COMPLEX,       INTENT(IN), DEVICE :: ab_dev(:,:,:)
    CLASS(t_mat),  INTENT(INOUT)      :: hmat
    !     ..
    !     .. Local Variables
    COMPLEX, ALLOCATABLE, DEVICE :: ab1_dev(:,:), ab2_dev(:,:)
    REAL rchi
    INTEGER  l,ll,m,i,j,istat

    CALL nvtxStartRange("hsmt_nonsph",1)    

    ALLOCATE(ab1_dev(lapw%nv(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))
    IF (iintsp.NE.jintsp) ALLOCATE(ab2_dev(lapw%nv(iintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))
    rchi=MERGE(REAL(chi),REAL(chi)*2,(atoms%invsat(na)==0))
    IF (iintsp==jintsp) THEN
       !Even if iintsp=jintsp=2 in the calling subroutine (hsmt_nonsph),
       !the local count starts with 1
       CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab_dev(1,1,1),SIZE(ab_dev,1),&
                   h_loc_dev,SIZE(h_loc_dev,1),CMPLX(0.,0.),ab1_dev,SIZE(ab1_dev,1))
       CALL zherk("U","N",lapw%nv(iintsp),ab_size,Rchi,ab1_dev,SIZE(ab1_dev,1),1.0,hmat%data_c,SIZE(hmat%data_c,1))
       istat = cudaDeviceSynchronize() 
    ELSE  !here the l_ss off-diagonal part starts
       CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab_dev(1,1,jintsp),SIZE(ab_dev,1),&
                   h_loc_dev,SIZE(h_loc_dev,1),CMPLX(0.,0.),ab1_dev,SIZE(ab1_dev,1))
       CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab_dev(1,1,iintsp),SIZE(ab_dev,1),&
                   h_loc_dev,SIZE(h_loc_dev,1),CMPLX(0.,0.),ab2_dev,SIZE(ab2_dev,1))
       !$cuf kernel do<<<*,256>>>
       DO i = 1,size(ab1_dev,2)
         DO j = 1,size(ab1_dev,1)
            ab1_dev(j,i) = conjg(ab1_dev(j,i))
         ENDDO
       ENDDO
       CALL zgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,ab2_dev,SIZE(ab2_dev,1),&
                   ab1_dev,SIZE(ab1_dev,1),CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
    ENDIF

    CALL nvtxEndRange
 END SUBROUTINE hsmt_nonsph_noMPI_gpu
#endif

  
END MODULE m_hsmt_nonsph
