MODULE m_hubbard1_setup

   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_uj2f
   USE m_doubleCounting
   USE m_hubbard1Distance
   USE m_occmtx
   USE m_hubbard1_io
   USE m_types_selfen
   USE m_add_selfen
   USE m_mpi_bc_tool
   USE m_greensf_io
#ifdef CPP_EDSOLVER
   USE EDsolver, only: EDsolver_from_cfg
#endif
#ifdef CPP_MPI
   use mpi
#endif
   IMPLICIT NONE

#include"cpp_double.h"

   CHARACTER(len=30), PARAMETER :: hubbard1CalcFolder = "Hubbard1"
   CHARACTER(len=30), PARAMETER :: hubbard1Outfile    = "out"

   CONTAINS

   SUBROUTINE hubbard1_setup(atoms,gfinp,hub1inp,input,fmpi,noco,pot,gdft,hub1data,results,den)

      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_gfinp),    INTENT(IN)     :: gfinp
      TYPE(t_hub1inp),  INTENT(IN)     :: hub1inp
      TYPE(t_input),    INTENT(IN)     :: input
      TYPE(t_mpi),      INTENT(IN)     :: fmpi
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_potden),   INTENT(IN)     :: pot
      TYPE(t_greensf),  INTENT(IN)     :: gdft(:) !green's function calculated from the Kohn-Sham system
      TYPE(t_hub1data), INTENT(INOUT)  :: hub1data
      TYPE(t_results),  INTENT(INOUT)  :: results
      TYPE(t_potden),   INTENT(INOUT)  :: den

      INTEGER :: i_hia,nType,l,occDFT_INT,ispin,m,i_exc,n
      INTEGER :: io_error,ierr
      INTEGER :: indStart,indEnd
      INTEGER :: hubbardioUnit
      INTEGER :: n_hia_task,extra,i_hia_start,i_hia_end
      REAL    :: U,J
      LOGICAL :: l_firstIT_HIA,l_ccfexist,l_bathexist,l_amf

      CHARACTER(len=300) :: folder
      TYPE(t_greensf),ALLOCATABLE :: gu(:)
      TYPE(t_selfen), ALLOCATABLE :: selfen(:)

#ifdef CPP_HDF
      INTEGER(HID_T)     :: greensf_fileID
#endif

      REAL    :: mu_dc(input%jspins)
      REAL    :: f0(atoms%n_hia),f2(atoms%n_hia)
      REAL    :: f4(atoms%n_hia),f6(atoms%n_hia)
      REAL    :: occDFT(atoms%n_hia,input%jspins)
      COMPLEX :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,3)
      COMPLEX, ALLOCATABLE :: e(:)
      COMPLEX, ALLOCATABLE :: ctmp(:)

      !Check if the EDsolver library is linked
#ifndef CPP_EDSOLVER
      CALL juDFT_error("No solver linked for Hubbard 1", hint="Link the edsolver library",calledby="hubbard1_setup")
#endif

      IF(fmpi%irank.EQ.0) THEN
         !-------------------------------------------
         ! Create the Input for the Hubbard 1 Solver
         !-------------------------------------------
         CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(hubbard1CalcFolder)))

         !Positions of the DFT+HIA elements in all DFT+U related arrays
         indStart = atoms%n_u+1
         indEnd   = atoms%n_u+atoms%n_hia

         ! calculate slater integrals from u and j
         CALL uj2f(input%jspins,atoms%lda_u(indStart:indEnd),atoms%n_hia,f0,f2,f4,f6)

         DO i_hia = 1, atoms%n_hia

            l = atoms%lda_u(atoms%n_u+i_hia)%l
            U = atoms%lda_u(atoms%n_u+i_hia)%U
            J = atoms%lda_u(atoms%n_u+i_hia)%J
            l_amf = atoms%lda_u(atoms%n_u+i_hia)%l_amf
            nType = atoms%lda_u(atoms%n_u+i_hia)%atomType

            IF(ALL(ABS(gdft(i_hia)%gmmpMat).LT.1e-12)) THEN
               CALL juDFT_error("Hubbard-1 has no DFT greensf available",calledby="hubbard1_setup")
            ENDIF

            !Create Subfolder (if there are multiple Hubbard 1 procedures)
            CALL hubbard1_path(atoms,i_hia,folder)
            CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(folder)))

            !-------------------------------------------------------
            ! Calculate the DFT occupation of the correlated shell
            !-------------------------------------------------------
            CALL occmtx(gdft(i_hia),gfinp,input,atoms,mmpMat(:,:,i_hia,:))

            !For the first iteration we can fix the occupation and magnetic moments in the inp.xml file
            l_firstIT_HIA = hub1data%iter.EQ.1 .AND.ALL(ABS(den%mmpMat(:,:,indStart:indEnd,:)).LT.1e-12)
            IF(l_firstIT_HIA) THEN
               IF(hub1inp%init_occ(i_hia) > -9e98) THEN
                  occDFT(i_hia,1) = MIN(2.0*l+1.0,hub1inp%init_occ(i_hia))
                  IF(input%jspins.EQ.2) occDFT(i_hia,2) = MAX(0.0,hub1inp%init_occ(i_hia)-2.0*l-1.0)
               ELSE
                  occDFT(i_hia,:) = 0.0
                  DO ispin = 1, input%jspins
                     DO m = -l, l
                        occDFT(i_hia,ispin) = occDFT(i_hia,ispin) + REAL(mmpMat(m,m,i_hia,ispin))
                     ENDDO
                  ENDDO
               ENDIF

               DO i_exc = 1, hub1inp%n_exc(i_hia)
                  IF(hub1inp%init_mom(i_hia,i_exc) > -9e98) THEN
                     hub1data%mag_mom(i_hia,i_exc) = hub1inp%init_mom(i_hia,i_exc)
                  ENDIF
               ENDDO
            ELSE
               occDFT(i_hia,:) = 0.0
               DO ispin = 1, input%jspins
                  DO m = -l, l
                     occDFT(i_hia,ispin) = occDFT(i_hia,ispin) + REAL(mmpMat(m,m,i_hia,ispin))
                  ENDDO
               ENDDO
            ENDIF
            !Nearest Integer occupation
            occDFT_INT = ANINT(SUM(occDFT(i_hia,:)))

            !Initial Information (We are already on irank 0)
            WRITE(oUnit,9010) nType
9010        FORMAT(/,"Setup for Hubbard 1 solver for atom ", I3, ": ")
            WRITE(oUnit,"(A)") "Everything related to the solver (e.g. mu_dc) is given in eV"
            WRITE(oUnit,"(/,A)") "Occupation from DFT-Green's function:"
            WRITE(oUnit,9020) 'spin-up','spin-dn'
9020        FORMAT(TR8,A7,TR3,A7)
            WRITE(oUnit,9030) occDFT(i_hia,:)
9030        FORMAT(TR7,f8.4,TR2,f8.4)

            !--------------------------------------------------------------------------
            ! Calculate the chemical potential for the solver
            ! This is equal to the double-counting correction used in DFT+U
            !--------------------------------------------------------------------------
            ! V_FLL = U (n - 1/2) - J (n - 1) / 2
            ! V_AMF = U n/2 + 2l/[2(2l+1)] (U-J) n
            !--------------------------------------------------------------------------
            mu_dc = doubleCountingPot(U,J,l,l_amf,.NOT.hub1inp%l_dftspinpol,occDFT(i_hia,:),&
                                      l_write=fmpi%irank==0)

            !-------------------------------------------------------
            ! Check for additional input files
            !-------------------------------------------------------
            !Is a crystal field matrix present in the work directory (overwrites the calculated matrix)
            INQUIRE(file=TRIM(ADJUSTL(cfg_file_ccf)),exist=l_ccfexist)
            IF(l_ccfexist) CALL read_ccfmat(hub1data%ccfmat(i_hia,-l:l,-l:l),l)
            !Is a bath parameter file present
            INQUIRE(file=TRIM(ADJUSTL(cfg_file_bath)),exist=l_bathexist)
            !Copy the bath file to the Hubbard 1 solver if its present
            IF(l_bathexist) CALL SYSTEM('cp ' // TRIM(ADJUSTL(cfg_file_bath)) // ' ' // TRIM(ADJUSTL(folder)))

            !-------------------------------------------------------
            ! Write the main config files
            !-------------------------------------------------------
            CALL write_hubbard1_input(folder,i_hia,l,f0(i_hia),f2(i_hia),f4(i_hia),f6(i_hia),&
                                      hub1inp,hub1data,mu_dc(1),occDFT_INT,l_bathexist,l_firstIT_HIA)
         ENDDO
      ENDIF !fmpi%irank == 0

      IF(fmpi%irank.EQ.0) THEN
         WRITE(*,*) "Calculating new density matrix ..."
      ENDIF

      !Argument order different because occDFT is not allocatable
      CALL mpi_bc(0,fmpi%mpi_comm,occDFT)

      !Initializations
      ALLOCATE(gu(atoms%n_hia))
      ALLOCATE(selfen(atoms%n_hia))
      DO i_hia = 1, atoms%n_hia
         CALL gu(i_hia)%init(gdft(i_hia)%elem,gfinp,atoms,input,contour_in=gdft(i_hia)%contour)
         CALL selfen(i_hia)%init(lmaxU_const,gdft(i_hia)%contour%nz,input%jspins,&
                                 noco%l_noco.AND.(noco%l_soc.OR.gfinp%l_mperp).OR.hub1inp%l_fullmatch)
      ENDDO

#ifdef CPP_MPI
      !distribute the individual hubbard1 elements over the ranks
      n_hia_task = FLOOR(REAL(atoms%n_hia)/(fmpi%isize))
      extra = atoms%n_hia - n_hia_task*fmpi%isize
      i_hia_start = fmpi%irank*n_hia_task + 1 + extra
      i_hia_end   =(fmpi%irank+1)*n_hia_task   + extra
      IF(fmpi%irank < extra) THEN
         i_hia_start = i_hia_start - (extra - fmpi%irank)
         i_hia_end   = i_hia_end   - (extra - fmpi%irank - 1)
      ENDIF
#else
      i_hia_start = 1
      i_hia_end   = atoms%n_hia
#endif

#ifdef CPP_MPI
      !Make sure that the ranks are synchronized
      CALL MPI_BARRIER(fmpi%mpi_comm,ierr)
#endif

      mmpMat = cmplx_0
      !------------------------------------------------------------
      ! This loop runs the solver
      !------------------------------------------------------------
      DO i_hia = i_hia_start, i_hia_end

         IF(i_hia > atoms%n_hia .OR. i_hia < 1) CYCLE

         nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
         l = atoms%lda_u(atoms%n_u+i_hia)%l

         CALL hubbard1_path(atoms,i_hia,folder)

         ALLOCATE(e(gdft(i_hia)%contour%nz),source=cmplx_0)

         CALL timestart("Hubbard 1: EDsolver")
         !We have to change into the Hubbard1 directory so that the solver routines can read the config
         CALL CHDIR(TRIM(ADJUSTL(folder)))
#ifdef CPP_EDSOLVER
         !Open the output file for the solver
         hubbardioUnit = 4000+i_hia
         OPEN(unit=hubbardioUnit, file=TRIM(ADJUSTL(hubbard1Outfile)),&
              status="replace", action="write", iostat=io_error)
         IF(io_error/=0) CALL juDFT_error("Error in opening EDsolver out file",calledby="hubbard1_setup")
         e = gdft(i_hia)%contour%e*hartree_to_ev_const
         CALL EDsolver_from_cfg(2*(2*l+1),gdft(i_hia)%contour%nz,e,selfen(i_hia)%data(:,:,:,1),1,hubbardioUnit)
         !---------------------------------------------------
         ! Calculate selfenergy on lower contour explicitly
         ! Mainly out of paranoia :D
         ! No rediagonalization (last argument switches this)
         !---------------------------------------------------
         e = conjg(gdft(i_hia)%contour%e)*hartree_to_ev_const
         CALL EDsolver_from_cfg(2*(2*l+1),gdft(i_hia)%contour%nz,e,selfen(i_hia)%data(:,:,:,2),0,hubbardioUnit)
         CLOSE(hubbardioUnit, iostat=io_error)
         IF(io_error/=0) CALL juDFT_error("Error in closing EDsolver out file",calledby="hubbard1_setup")
#endif
         IF(atoms%n_hia>1) THEN
            CALL CHDIR("../../")
         ELSE
            CALL CHDIR("../")
         ENDIF
         CALL timestop("Hubbard 1: EDsolver")

         DEALLOCATE(e)

         !-------------------------------------------
         ! Postprocess selfenergy
         !-------------------------------------------
         CALL selfen(i_hia)%postProcess(input%jspins,pot%mmpMat(:,:,atoms%n_u+i_hia,:))

         !----------------------------------------------------------------------
         ! Solution of the Dyson Equation
         !----------------------------------------------------------------------
         ! G(z)^(-1) = G_0(z)^(-1) - mu - Sigma(z)
         !----------------------------------------------------------------------
         ! Sigma(z) is the self-energy from the impurity solver
         ! We introduce an additional chemical potential mu, which is determined
         ! so that the occupation of the correlated orbital does not change
         !----------------------------------------------------------------------
#ifdef CPP_DEBUG
         OPEN(unit=1337, file=TRIM(ADJUSTL(folder)) // 'mu',&
              status="replace", action="write", iostat=io_error)
#endif

         CALL timestart("Hubbard 1: Add Selfenergy")
         CALL add_selfen(gdft(i_hia),selfen(i_hia),gfinp,input,atoms,&
                         occDFT(i_hia,:),gu(i_hia),mmpMat(:,:,i_hia,:))
         CALL timestop("Hubbard 1: Add Selfenergy")

#ifdef CPP_DEBUG
         CLOSE(unit=1337)
#endif

      ENDDO

      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*)
         WRITE(oUnit,'(A)') "Calculated mu to match Self-energy to DFT-GF"
      ENDIF
      !Collect the impurity Green's Function
      DO i_hia = 1, atoms%n_hia
         CALL gu(i_hia)%collect(fmpi%mpi_comm)
         CALL selfen(i_hia)%collect(fmpi%mpi_comm)
         IF(fmpi%irank.EQ.0) THEN
            !We found the chemical potential to within the desired accuracy
            WRITE(oUnit,*) 'i_hia: ',i_hia, "    muMatch = ", selfen(i_hia)%muMatch(:)
         ENDIF
      ENDDO


#ifdef CPP_HDF
      IF(fmpi%irank.EQ.0) THEN
         !------------------------------
         !Write out DFT Green's Function
         !------------------------------
         CALL timestart("Hubbard 1: IO/Write")
         CALL openGreensFFile(greensf_fileID, input, gfinp, atoms, inFilename="greensf_DFT.hdf")
         CALL writeGreensFData(greensf_fileID, input, gfinp, atoms, &
                               GREENSF_HUBBARD_CONST, gdft, mmpmat)
         CALL closeGreensFFile(greensf_fileID)

         !-------------------------------------
         !Write out correlated Green's Function
         !-------------------------------------
         CALL openGreensFFile(greensf_fileID, input, gfinp, atoms, inFilename="greensf_IMP.hdf")
         CALL writeGreensFData(greensf_fileID, input, gfinp, atoms, &
                              GREENSF_HUBBARD_CONST, gu, mmpmat,selfen=selfen)
         CALL closeGreensFFile(greensf_fileID)
         CALL timestop("Hubbard 1: IO/Write")
      ENDIF
#endif


#ifdef CPP_MPI
      !Collect the density matrix to rank 0
      n = SIZE(mmpMat)
      ALLOCATE(ctmp(n))
      CALL MPI_REDUCE(mmpMat,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF(fmpi%irank.EQ.0) CALL CPP_BLAS_ccopy(n,ctmp,1,mmpMat,1)
      DEALLOCATE(ctmp)
#endif

      !--------------------------------------------------------------------
      ! Calculate Distances from last density matrix and update den%mmpmat
      !--------------------------------------------------------------------
      results%last_mmpMatdistance = 0.0
      results%last_occdistance = 0.0

      IF(fmpi%irank.EQ.0) THEN
         DO i_hia = 1, atoms%n_hia
            CALL hubbard1Distance(den%mmpMat(:,:,atoms%n_u+i_hia,:),mmpMat(:,:,i_hia,:),results)
            DO ispin = 1, MERGE(3,input%jspins,gfinp%l_mperp)
               den%mmpMat(-lmaxU_const:,-lmaxU_const:,atoms%n_u+i_hia,ispin) = mmpMat(-lmaxU_const:,-lmaxU_const:,i_hia,ispin)
            ENDDO
         ENDDO
      ENDIF

      !Broadcast the density matrix
      CALL mpi_bc(den%mmpMat,0,fmpi%mpi_comm)
      CALL mpi_bc(results%last_occdistance,0,fmpi%mpi_comm)
      CALL mpi_bc(results%last_mmpMatdistance,0,fmpi%mpi_comm)

      IF(fmpi%irank.EQ.0) THEN
         WRITE(*,*) "Hubbard 1 Iteration: ", hub1data%iter
         WRITE(*,*) "Distances: "
         WRITE(*,*) "-----------------------------------------------------"
         WRITE(*,*) "Occupation Distance: " , results%last_occdistance
         WRITE(*,*) "Element Distance:    " , results%last_mmpMatdistance
         WRITE(*,*) "-----------------------------------------------------"
         WRITE(oUnit,*) "nmmp occupation distance: ", results%last_occdistance
         WRITE(oUnit,*) "nmmp element distance:    ", results%last_mmpMatdistance
         WRITE(oUnit,FMT=8140) hub1data%iter
8140     FORMAT (/,5x,'******* Hubbard 1 it=',i3,'  is completed********',/,/)
      ENDIF

   END SUBROUTINE hubbard1_setup

   SUBROUTINE hubbard1_path(atoms,i_hia,xPath)

      !Defines the folder structure
      ! The Solver is run in the subdirectories
      ! Hubbard1/ if only one Hubbard1 procedure is run
      ! Hubbard1/atom_label_l if there are more

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      INTEGER,             INTENT(IN)  :: i_hia
      CHARACTER(len=300),  INTENT(OUT) :: xPath

      CHARACTER(len=300) :: folder,fmt
      CHARACTER(len=1),PARAMETER :: spdfg(0:4) = ['s','p','d','f','g']
      INTEGER nType,l

      nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
      l = atoms%lda_u(atoms%n_u+i_hia)%l
      xPath = TRIM(ADJUSTL(hubbard1CalcFolder))
      IF(atoms%n_hia>1) THEN
         WRITE(fmt,'("(A",I2.2,",A1,A1,A1)")') LEN(TRIM(ADJUSTL(atoms%label(nType))))
         WRITE(folder,fmt) TRIM(ADJUSTL(atoms%label(nType))),"_",spdfg(l),"/"
      ELSE
         folder=""
      ENDIF
      xPath = TRIM(ADJUSTL(xPath)) // "/" // folder

   END SUBROUTINE hubbard1_path

END MODULE m_hubbard1_setup
