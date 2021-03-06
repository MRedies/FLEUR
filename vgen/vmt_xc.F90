   !--------------------------------------------------------------------------------
   ! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
   ! This file is part of FLEUR and available as free software under the conditions
   ! of the MIT license as expressed in the LICENSE file in more detail.
   !--------------------------------------------------------------------------------
MODULE m_vmt_xc
#ifdef CPP_MPI 
   use mpi 
#endif
   USE m_judft
      !.....------------------------------------------------------------------
      !     Calculate the GGA xc-potential in the MT-spheres
      !.....------------------------------------------------------------------
      !     instead of vmtxcor.f: the different exchange-correlation
      !     potentials defined through the key icorr are called through
      !     the driver subroutine vxcallg.f, subroutines vectorized
      !     ** r.pentcheva 22.01.96
      !     *********************************************************
      !     angular mesh calculated on speacial gauss-legendre points
      !     in order to use orthogonality of lattice harmonics and
      !     avoid a least square fit
      !     ** r.pentcheva 04.03.96
      !     *********************************************************
      !     MPI and OpenMP parallelization
      !             U.Alekseeva, February 2017
      !     *********************************************************

   CONTAINS
      SUBROUTINE vmt_xc(fmpi,sphhar,atoms,&
                        den,xcpot,input,sym,EnergyDen,kinED,noco,vTot,vx,exc)
#include"cpp_double.h"
         use m_libxc_postprocess_gga
         USE m_mt_tofrom_grid
         USE m_types_xcpot_inbuild
         USE m_types
         USE m_metagga
         USE m_juDFT_string
         IMPLICIT NONE

         CLASS(t_xcpot),INTENT(IN)      :: xcpot
         TYPE(t_mpi),INTENT(IN)         :: fmpi
         TYPE(t_input),INTENT(IN)       :: input
         TYPE(t_sym),INTENT(IN)         :: sym
         TYPE(t_sphhar),INTENT(IN)      :: sphhar
         TYPE(t_atoms),INTENT(IN)       :: atoms
         TYPE(t_potden),INTENT(IN)      :: den,EnergyDen
         TYPE(t_noco), INTENT(IN)       :: noco
         TYPE(t_potden),INTENT(INOUT)   :: vTot,vx,exc
         TYPE(t_kinED),INTENT(IN)       :: kinED
         !     ..
         !     .. Local Scalars ..
         TYPE(t_gradients)     :: grad, tmp_grad
         TYPE(t_xcpot_inbuild) :: xcpot_tmp
         TYPE(t_potden)        :: vTot_tmp
         TYPE(t_sphhar)        :: tmp_sphhar
         REAL, ALLOCATABLE     :: ch(:,:),v_x(:,:),v_xc(:,:),e_xc(:,:)
         INTEGER               :: n,nsp,nt,jr, loc_n
         INTEGER               :: i, j, idx, cnt
         REAL                  :: divi

         !     ..

         !locals for fmpi
         integer :: ierr
         integer:: n_start,n_stride
         REAL,ALLOCATABLE:: xcl(:,:)
         LOGICAL :: lda_atom(atoms%ntype),l_libxc, perform_MetaGGA
         !.....------------------------------------------------------------------
         perform_MetaGGA = ALLOCATED(EnergyDen%mt) &
                         .AND. (xcpot%exc_is_MetaGGA() .or. xcpot%vx_is_MetaGGA())
         lda_atom=.FALSE.; l_libxc=.FALSE.
         SELECT TYPE(xcpot)
         TYPE IS(t_xcpot_inbuild)
            lda_atom=xcpot%lda_atom
            IF (ANY(lda_atom)) THEN
               IF((.NOT.xcpot%is_name("pw91"))) &
                  CALL judft_warn("Using locally LDA only possible with pw91 functional")
               !TODO: check this code and the functionality
               !xcpot_tmp%inbuild_name="l91"
               !xcpot_tmp%l_relativistic=.FALSE.
               CALL xcpot_tmp%init(atoms%ntype)
               ALLOCATE(xcl(SIZE(v_xc,1),SIZE(v_xc,2)))
            ENDIF
         CLASS DEFAULT
            l_libxc=.true. !libxc!!
         END SELECT

         nsp=atoms%nsp()
         !ALLOCATE(ch(nsp*atoms%jmtd,input%jspins),v_x(nsp*atoms%jmtd,input%jspins),v_xc(nsp*atoms%jmtd,input%jspins),e_xc(nsp*atoms%jmtd,input%jspins))
         !IF (xcpot%needs_grad()) CALL xcpot%alloc_gradients(SIZE(ch,1),input%jspins,grad)

         CALL init_mt_grid(input%jspins,atoms,sphhar,xcpot%needs_grad(),sym)

#ifdef CPP_MPI
         n_start=fmpi%irank+1
         n_stride=fmpi%isize
         IF (fmpi%irank>0) THEN
            vTot%mt=0.0
            vx%mt=0.0
            exc%mt=0.0
         ENDIF
#else
         n_start=1
         n_stride=1
#endif
         loc_n = 0
         !TODO: MetaGGA
         DO n = n_start,atoms%ntype,n_stride
            ALLOCATE(ch(nsp*atoms%jri(n),input%jspins),v_x(nsp*atoms%jri(n),input%jspins),&
                     v_xc(nsp*atoms%jri(n),input%jspins),e_xc(nsp*atoms%jri(n),input%jspins))
            IF (xcpot%needs_grad()) CALL xcpot%alloc_gradients(SIZE(ch,1),input%jspins,grad)
            loc_n = loc_n + 1

            CALL mt_to_grid(xcpot%needs_grad(), input%jspins, atoms,sym,sphhar,.True.,den%mt(:,0:,n,:),n,noco,grad,ch)

            !
            !         calculate the ex.-cor. potential
#ifdef CPP_LIBXC
            if(perform_MetaGGA .and. kinED%set) then
              CALL xcpot%get_vxc(input%jspins,ch,v_xc&
                   , v_x,grad, kinEnergyDen_KS=kinED%mt(:,:,loc_n))
            else
               CALL xcpot%get_vxc(input%jspins,ch,v_xc&
                  , v_x,grad)
            endif
#else
               CALL xcpot%get_vxc(input%jspins,ch,v_xc&
                  , v_x,grad)
#endif
            IF (lda_atom(n)) THEN
               ! Use local part of pw91 for this atom
               CALL xcpot_tmp%get_vxc(input%jspins,ch(:nsp*atoms%jri(n),:),xcl(:nsp*atoms%jri(n),:),v_x(:nsp*atoms%jri(n),:),grad)
               !Mix the potentials
               divi = 1.0 / (atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(1,n))
               nt=0
               DO jr=1,atoms%jri(n)
                  v_xc(nt+1:nt+nsp,:) = ( xcl(nt+1:nt+nsp,:) * ( atoms%rmsh(atoms%jri(n),n) &
                                          - atoms%rmsh(jr,n) ) &
                                          + v_xc(nt+1:nt+nsp,:) * ( atoms%rmsh(jr,n) &
                                          - atoms%rmsh(1,n) ) &
                                         ) * divi
                  nt=nt+nsp
               ENDDO
            ENDIF

            !Add postprocessing for libxc
            IF (l_libxc.AND.xcpot%needs_grad()) CALL libxc_postprocess_gga_mt(xcpot,atoms,sym,sphhar,noco,n,v_xc,grad, atom_num=n)

            CALL mt_from_grid(atoms,sym,sphhar,n,input%jspins,v_xc,vTot%mt(:,0:,n,:))
            CALL mt_from_grid(atoms,sym,sphhar,n,input%jspins,v_x,vx%mt(:,0:,n,:))

            IF (ALLOCATED(exc%mt)) THEN
               !
               !           calculate the ex.-cor energy density
               !
#ifdef CPP_LIBXC
               IF(perform_MetaGGA .and. kinED%set) THEN
                  CALL xcpot%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),&
                     e_xc(:nsp*atoms%jri(n),1),grad, &
                     kinEnergyDen_KS=kinED%mt(:,:,loc_n), mt_call=.True.)
               ELSE
                  CALL xcpot%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),&
                     e_xc(:nsp*atoms%jri(n),1),grad, mt_call=.True.)
               ENDIF
#else
               CALL xcpot%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),&
                       e_xc(:nsp*atoms%jri(n),1),grad, mt_call=.True.)
#endif
               !write (*,*) "cut first ", cut_ratio, " number of points"
               !where(cut_mask) e_xc(:,1) = 0.0

               IF (lda_atom(n)) THEN
                  ! Use local part of pw91 for this atom
                  CALL xcpot_tmp%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),xcl(:nsp*atoms%jri(n),1),grad)
                  !Mix the potentials
                  nt=0
                  DO jr=1,atoms%jri(n)
                     e_xc(nt+1:nt+nsp,1) = ( xcl(nt+1:nt+nsp,1) * ( atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(jr,n) ) +&
                                            e_xc(nt+1:nt+nsp,1) * ( atoms%rmsh(jr,n) - atoms%rmsh(1,n) ) ) * divi
                     nt=nt+nsp
                  END DO
               ENDIF
               CALL mt_from_grid(atoms,sym,sphhar,n,1,e_xc,exc%mt(:,0:,n,:))
            ENDIF
            DEALLOCATE (ch,v_x,v_xc,e_xc)
         ENDDO

         CALL finish_mt_grid()
#ifdef CPP_MPI
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,vx%mt,SIZE(vx%mt),CPP_MPI_REAL,MPI_SUM,fmpi%mpi_comm,ierr)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,vTot%mt,SIZE(vTot%mt),CPP_MPI_REAL,MPI_SUM,fmpi%mpi_comm,ierr)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,exc%mt,SIZE(exc%mt),CPP_MPI_REAL,MPI_SUM,fmpi%mpi_comm,ierr)
#endif
         !
         RETURN
      END SUBROUTINE vmt_xc
   END MODULE m_vmt_xc
