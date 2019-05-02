MODULE m_hubbard1_setup

   USE m_juDFT

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE hubbard1_setup(iterHIA,atoms,hub1,sym,mpi,noco,input,inDen,usdus,pot,gdft,l_runinfleur,l_runhia,results)

      USE m_types
      USE m_hubbard1_io
      USE m_uj2f
      USE m_constants
      USE m_gfcalc
      USE m_umtx
      USE m_vmmp

      INTEGER,          INTENT(INOUT)  :: iterHIA !number of iteration 
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_hub1ham),  INTENT(IN)     :: hub1
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_mpi),      INTENT(IN)     :: mpi
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_input),    INTENT(IN)     :: input
      TYPE(t_usdus),    INTENT(IN)     :: usdus
      TYPE(t_potden),   INTENT(IN)     :: inDen
      TYPE(t_potden),   INTENT(INOUT)  :: pot
      TYPE(t_results),  INTENT(INOUT)  :: results
      TYPE(t_greensf),  INTENT(IN)     :: gdft !green's function in the mt-sphere including the potential correction form dft+hubbard1 
      LOGICAL,          INTENT(IN)     :: l_runinfleur !Determines wether we call the the solver here or run separately
      LOGICAL,          INTENT(IN)     :: l_runhia
      
      INTEGER i_hia,i_gf,n,l,n_occ,ispin,m,matsize,i,iz,N_mu,i_mu,j,io_error
      REAL mu_dc,beta,e_lda_hia
      CHARACTER(len=300) :: cwd
      CHARACTER(len=300) :: path
      CHARACTER(len=300) :: folder
      CHARACTER(len=300) :: message
      TYPE(t_greensf)    :: gu

      REAL     f0(hub1%n_hia,input%jspins),f2(hub1%n_hia,input%jspins)
      REAL     f4(hub1%n_hia,input%jspins),f6(hub1%n_hia,input%jspins)
      REAL     u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
               -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia)

      COMPLEX  mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,input%jspins)
      COMPLEX  mmpMat_in(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,input%jspins)
      COMPLEX  selfen(hub1%n_hia,gdft%nz,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1))
      REAL     e(gdft%nz),dist(input%jspins)
      REAL     n_l(hub1%n_hia,input%jspins)
      LOGICAL  l_exist




      IF(ANY(inDen%mmpMat(:,:,atoms%n_u+1:hub1%n_hia+atoms%n_u,:).NE.0.0).OR.l_runhia) THEN
         !Get the slater integrals from the U and J parameters
         CALL uj2f(input%jspins,hub1%lda_u(:),hub1%n_hia,f0,f2,f4,f6)
         f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
         f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
         f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
         f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2
         CALL umtx(hub1%lda_u(:),hub1%n_hia,f0(:,1),f2(:,1),f4(:,1),f6(:,1),u)

         mmpMat = inDen%mmpMat(:,:,atoms%n_u+1:hub1%n_hia,:)
         CALL v_mmp(sym,atoms,hub1%lda_u(:),hub1%n_hia,input%jspins,.true.,mmpMat,&
         u,f0,f2,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+hub1%n_hia,:),e_lda_hia)

         IF(l_runhia.AND.ANY(gdft%gmmpMat(:,:,:,:,:,:,:).NE.0.0)) THEN 
            IF(gdft%mode.EQ.2) CALL juDFT_error("This energy contour is not supported at the moment for DFT+Hubbard1",calledby="hubbard1_setup")
            !The onsite green's function was calculated but the solver 
            !was not yet run
            !--> write out the configuration for the hubbard 1 solver 

            CALL gu%init(input,lmaxU_const,atoms,.true.,noco,nz_in=gdft%nz, e_in=gdft%e,de_in=gdft%de,matsub_in=gdft%nmatsub)

            !Get the working directory
            CALL get_environment_variable('PWD',cwd)
            !Create a folder where to store the files from the solver
            !Is this applicable on all systems where fleur can be run?
            iterHIA = iterHIA + 1
            WRITE(folder,"(A5,I3.3)") "hub1_" , iterHIA
            path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
            CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)))

            DO i_hia = 1, hub1%n_hia
               n = hub1%lda_u(i_hia)%atomType
               l = hub1%lda_u(i_hia)%l
               CALL indexgf(atoms,l,n,i_gf)
               matsize = 2*(2*l+1)
               !Create a subdirectory for the atomType and shell
               WRITE(folder,"(A4,I2.2,A2,I1.1)") "atom",n,"_l",l
               CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)))
            
               !calculate the occupation of the correlated shell
               CALL occmtx(gdft,i_gf,atoms,sym,input%jspins,input%onsite_beta*hartree_to_ev_const,results%ef,mmpMat(:,:,i_hia,:))
               n_l(i_hia,:) = 0.0
               DO ispin = 1, input%jspins
                  DO m = -l, l
                     n_l(i_hia,ispin) = n_l(i_hia,ispin) + mmpMat(m,m,i_hia,ispin)
                  ENDDO
               ENDDO
               WRITE(*,*) "OCCUPATION: ", SUM(n_l(i_hia,:))

               !calculate the chemical potential
               CALL mudc(hub1%lda_u(i_hia)%U,hub1%lda_u(i_hia)%J,l,n_l(i_hia,:),hub1%lda_u(i_hia)%l_amf,mu_dc,input%jspins)
               
               n_occ = ANINT(SUM(n_l(i_hia,:)))
               n_l(i_hia,1) = SUM(n_l(i_hia,:))
               !Check wether the hubbard 1 solver was run:
               INQUIRE(file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/" // "se.atom",exist=l_exist)

               IF(l_exist) THEN
                  CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,1:gdft%nz-gdft%nmatsub,1:2*(2*l+1),1:2*(2*l+1)),gdft%nz-gdft%nmatsub,2*(2*l+1),e(:),.false.)
                  IF(gdft%nmatsub.GT.0) CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,gdft%nz-gdft%nmatsub:gdft%nz,1:2*(2*l+1),1:2*(2*l+1)),gdft%nmatsub,2*(2*l+1),e(:),.true.)
               ELSE
                  CALL write_hubbard1_input(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),&
                                             hub1%xi(i_hia),-hub1%bz(i_hia),MAX(1,n_occ-hub1%n_exc),MIN(2*(2*l+1),n_occ+hub1%n_exc),hub1%beta,mu_dc,hub1%l_ccf(i_hia),&
                                             gdft%nz-gdft%nmatsub,gdft%nmatsub,REAL(gdft%e(1)),REAL(gdft%e(gdft%nz-gdft%nmatsub)),AIMAG(gdft%e(1)))
                  IF(hub1%l_ccf(i_hia)) THEN
                     CALL write_ccfmat(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",hub1%ccfmat(i_hia,-l:l,-l:l),l)
                  ENDIF
                  OPEN(unit=1337,file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/" // "contour.dat",status="replace",action="write")

                  DO iz = 1, gdft%nz
                     WRITE(1337,"(2f14.8)") REAL(gdft%e(iz))*hartree_to_ev_const, AIMAG(gdft%e(iz))*hartree_to_ev_const
                  ENDDO

                  CLOSE(unit= 1337)
                  
                  !There has to be a better solution
                  !Maybe use CALL System() to start the solver from here
                  !EXPERIMENTAL:
                  IF(l_runinfleur) THEN
                     CALL timestart("Hubbard 1 solver")
                     CALL CHDIR(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/")
                     WRITE(*,"(A)") "Running Hubbard 1 solver"
                     CALL SYSTEM("/home/henning/GIT/hub2new4sp/eigen > out")
                     CALL SYSTEM("/home/henning/GIT/hub2new4sp/selfen > outselfen") 
                     CALL SYSTEM("/home/henning/GIT/hub2new4sp/angmom > outangmom") 
                     WRITE(*,"(A)") "Hubbard 1 solver finished" 
                     CALL timestop("Hubbard 1 solver")
                     CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,1:gdft%nz-gdft%nmatsub,1:2*(2*l+1),1:2*(2*l+1)),gdft%nz-gdft%nmatsub,2*(2*l+1),e(:),.false.)
                     IF(gdft%nmatsub.GT.0) CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,gdft%nz-gdft%nmatsub:gdft%nz,1:2*(2*l+1),1:2*(2*l+1)),gdft%nmatsub,2*(2*l+1),e(:),.true.)
                  ENDIF
               ENDIF
            ENDDO
            CALL CHDIR(TRIM(ADJUSTL(cwd)))

            IF(l_exist.OR.l_runinfleur) THEN
               CALL add_selfen(gdft,gu,selfen,atoms,hub1,sym,input%onsite_beta,results%ef,input%jspins,n_l(:,1),mu_dc/hartree_to_ev_const,mmpMat)
               ! calculate potential matrix and total energy correction
               CALL v_mmp(sym,atoms,hub1%lda_u,hub1%n_hia,input%jspins,.true.,mmpMat,&
                     u,f0,f2,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+hub1%n_hia,:),results%e_ldau)
               !
               ! Output the density of states from the two green's functions
               !
               DO i_hia = 1, hub1%n_hia
                  n = hub1%lda_u(i_hia)%atomType
                  l = hub1%lda_u(i_hia)%l
                  WRITE(folder,"(A5,I3.3)") "hub1_" , iterHIA
                  path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
                  WRITE(folder,"(A4,I2.2,A2,I1.1)") "atom",n,"_l",l
                  CALL CHDIR(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/")
                  CALL indexgf(atoms,l,n,i_gf)
                  CALL ldosmtx("gdft",gdft,i_gf,atoms,sym,input%jspins)
                  CALL ldosmtx("g",gu,i_gf,atoms,sym,input%jspins)
               ENDDO
               CALL CHDIR(TRIM(ADJUSTL(cwd)))
               INQUIRE(file="n_mmpmat_hubbard1",exist = l_exist)
               IF(l_exist) THEN
                  OPEN(unit=1337,file="n_mmpmat_hubbard1",status="old",action="read",iostat=io_error)
                  IF(io_error.NE.0) CALL juDFT_error("IO-error in density matrix",calledby="hubbard1_setup")
                  READ(1337,"(I3.3,2f14.8)")
                  READ(1337,"(7f14.8)") mmpMat_in(:,:,:,:)
                  CLOSE(unit=1337)
                  CALL n_mmp_dist(mmpMat_in,mmpMat,hub1,results,input%jspins)
               ENDIF 
   
               !Write out the density matrix (because i don't wanna change inDen to intent(inout)
               OPEN(unit=1337,file="n_mmpmat_hubbard1",status="replace",action="write",iostat=io_error)
               IF(io_error.NE.0) CALL juDFT_error("IO-error in density matrix",calledby="hubbard1_setup")
               WRITE(1337,"(I3.3,2f14.8)") iterHIA,results%last_occdistance,results%last_mmpMatdistance
               WRITE(1337,"(7f14.8)") mmpMat(:,:,:,:)
               CLOSE(unit=1337)
               WRITE(*,"(7f14.8)") REAL(pot%mmpMat(:,:,:,:))
            ELSE
               WRITE(message,9010) iterHIA
               CALL juDFT_END(message)
9010           FORMAT("Hubbard 1 input written into hub1_",I3.3)
            ENDIF
         ENDIF
      ELSE
         !occupation matrix is zero and LDA+Hubbard 1 shouldn't be run yet
         !There is nothing to be done yet just set the potential correction to 0
         pot%mmpMat(:,:,atoms%n_u+1:hub1%n_hia+atoms%n_u,:) = CMPLX(0.0,0.0)
      ENDIF
      
     
   END SUBROUTINE hubbard1_setup

   SUBROUTINE mudc(U,J,l,n,l_amf,mu,jspins)

      REAL,    INTENT(IN)  :: U
      REAL,    INTENT(IN)  :: J 
      INTEGER, INTENT(IN)  :: l
      REAL,    INTENT(IN)  :: n(jspins)
      LOGICAL, INTENT(IN)  :: l_amf
      REAL,    INTENT(OUT) :: mu
      INTEGER, INTENT(IN)  :: jspins

      REAL vdcfll1,vdcfll2
      REAL vdcamf1,vdcamf2
      REAL nup,ndn

      IF(jspins.EQ.2) THEN
         nup = n(1)
         ndn = n(2)
      ELSE
         nup = 0.5 * n(1)
         ndn = nup
      ENDIF

      vdcfll1= u*(nup+ndn -0.5) - j*(nup-0.5)
      vdcfll2= u*(nup+ndn -0.5) - j*(ndn-0.5)
      vdcamf1= u*ndn+2.0*l/(2.0*l+1)*(u-j)*nup
      vdcamf2= u*nup+2.0*l/(2.0*l+1)*(u-j)*ndn
      WRITE(*,"(A)") 'Double counting chemical potential:'
      WRITE(*,9010) 'FLL: ','spin-up','spin-dn','(up+dn)/2','up-dn'
      WRITE(*,9020) vdcfll1,vdcfll2,(vdcfll1+vdcfll2)/2,vdcfll1-vdcfll2
      WRITE(*,9010) 'AMF: ','spin-up','spin-dn','(up+dn)/2','up-dn'
      WRITE(*,9020)  vdcamf1,vdcamf2,(vdcamf1+vdcamf2)/2,vdcamf1-vdcamf2

      IF(l_amf) THEN
         WRITE(*,"(A)") "Using the around-mean-field limit"
         mu = (vdcamf1+vdcamf2)/2
      ELSE
      WRITE(*,"(A)") "Using the fully-localized limit"
         mu = (vdcfll1+vdcfll2)/2
      ENDIF 
      WRITE(*,9030) mu

9010  FORMAT(TR3,A4,TR1,A7,TR3,A7,TR3,A9,TR3,A5)
9020  FORMAT(TR7,f8.4,TR2,f8.4,TR2,f8.4,TR4,f8.4)
9030  FORMAT(TR3,"mu = ",f7.4)
   END SUBROUTINE

   SUBROUTINE add_selfen(g,gp,selfen,atoms,hub1,sym,beta,ef,jspins,n_occ,mu_dc,mmpMat)

      !Calculates the interacting Green's function for the mt-sphere with
      !
      ! (G)^-1 = (G_0)^-1 - mu 1 - selfen
      !
      !The term mu * unity is there to ensure that the number of particles 
      !doesn't change and is determined by a two-step process
      !The occupation as a function of mu has a peak in the region where 
      !something is inside the energy interval between e_bot adn e_fermi
      !To determine where we have the same number of particles we first 
      !search for the maximum occupation
      !Then the desired chemical potential is found with the bisection method 
      !to the right of the maximum
      
      USE m_types
      USE m_constants
      USE m_gfcalc

      TYPE(t_greensf),  INTENT(IN)     :: g
      TYPE(t_greensf),  INTENT(INOUT)  :: gp
      TYPE(t_hub1ham),  INTENT(IN)     :: hub1
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym
      COMPLEX,          INTENT(IN)     :: selfen(hub1%n_hia,g%nz,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1))
      INTEGER,          INTENT(IN)     :: jspins
      REAL,             INTENT(IN)     :: beta
      REAL,             INTENT(IN)     :: ef
      REAL,             INTENT(IN)     :: n_occ(hub1%n_hia)
      REAL,             INTENT(IN)     :: mu_dc
      COMPLEX,          INTENT(OUT)    :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,jspins)

      INTEGER i_hia,l,nType,i_gf,ns,ispin,m,iz,ipm
      CHARACTER(len=6) app

      REAL mu_a,mu_b,mu_step,mu_max,n_max
      REAL mu,n
      LOGICAL l_mperp

      TYPE(t_mat) :: gmat,vmat

      !replace with noco%l_mperp
      l_mperp = .true.

      !Interval where we expect the correct mu
      mu_a = -2.0
      mu_b = 2.0
      mu_step = 0.05
      mu_max = 0.0
      n_max = 0.0

      DO i_hia = 1, hub1%n_hia
         l = hub1%lda_u(i_hia)%l
         nType = hub1%lda_u(i_hia)%atomType
         ns = 2*l+1
         CALL indexgf(atoms,l,nType,i_gf)
         !intialize the matrices
         CALL gmat%init(.false.,2*ns,2*ns)
         CALL vmat%init(.false.,2*ns,2*ns)
         !Search for the maximum of occupation
         OPEN(unit=1337,file="mu",status="replace",action="write")
         mu = mu_a
         DO WHILE(mu.LE.mu_b)
            mu = mu + mu_step
            DO iz = 1, g%nz
               vmat%data_c = selfen(i_hia,iz,1:2*ns,1:2*ns)
               IF(.NOT.l_mperp) THEN
                  !Dismiss spin-off-diagonal elements
                  vmat%data_c(1:ns,ns+1:2*ns) = 0.0
                  vmat%data_c(ns+1:2*ns,1:ns) = 0.0
               ENDIF
               DO ipm = 1, 2
                  CALL to_tmat(gmat,g%gmmpMat(1,iz,i_gf,:,:,:,ipm),jspins,2,l)
                  CALL add_pot(gmat,vmat,mu-mu_dc,(ipm.EQ.1))
                  CALL to_g(gmat,gp%gmmpMat(1,iz,i_gf,:,:,:,ipm),2,jspins,l)
               ENDDO
            ENDDO
            CALL occmtx(gp,i_gf,atoms,sym,jspins,beta*hartree_to_ev_const,ef,mmpMat(:,:,i_hia,:))
            !Calculate the trace
            n = 0.0
            DO ispin = 1, jspins
               DO m = -l, l
                  n = n + mmpMat(m,m,i_hia,ispin)
               ENDDO
            ENDDO
            WRITE(1337,*) mu,n
            IF(n.GT.n_max) THEN
               mu_max = mu
               n_max  = n
            ENDIF
         ENDDO
         CLOSE(1337)
         !Set up the interval for the bisection method (mu_max,mu_b)
         mu_a = mu_max
         DO 
            mu = (mu_a + mu_b)/2.0
            DO iz = 1, g%nz
               vmat%data_c = selfen(i_hia,iz,1:2*ns,1:2*ns)
               IF(.NOT.l_mperp) THEN
                  !Dismiss spin-off-diagonal elements
                  vmat%data_c(1:ns,ns+1:2*ns) = 0.0
                  vmat%data_c(ns+1:2*ns,1:ns) = 0.0
               ENDIF
               DO ipm = 1, 2
                  CALL to_tmat(gmat,g%gmmpMat(1,iz,i_gf,:,:,:,ipm),jspins,2,l)
                  CALL add_pot(gmat,vmat,mu-mu_dc,(ipm.EQ.1))
                  CALL to_g(gmat,gp%gmmpMat(1,iz,i_gf,:,:,:,ipm),2,jspins,l)
               ENDDO
            ENDDO
            CALL occmtx(gp,i_gf,atoms,sym,jspins,beta,ef,mmpMat(:,:,i_hia,:))
            !Calculate the trace
            n = 0.0
            DO ispin = 1, jspins
               DO m = -l, l
                  n = n + mmpMat(m,m,i_hia,ispin)
               ENDDO
            ENDDO
            IF(ABS(n-n_occ(i_hia)).LT.0.001.OR.ABS((mu_b - mu_a)/2.0).LT.0.00001) THEN
               !We found the chemical potential to within the desired accuracy
               !TODO: Write to output file
               WRITE(*,*) "Calculated mu: ", mu
               WRITE(*,*) "OCCUPATION: ", n
               EXIT
            ELSE IF((n - n_occ(i_hia)).GT.0) THEN
               !The occupation is to small --> choose the left interval
               mu_a = mu
            ELSE IF((n - n_occ(i_hia)).LT.0) THEN
               !The occupation is to big --> choose the right interval
               mu_b = mu
            ENDIF
         ENDDO
         CALL gmat%free()
         CALL vmat%free()
      ENDDO

   END SUBROUTINE add_selfen

   SUBROUTINE add_pot(gmat,vmat,mu,l_upper)

      USE m_types

      TYPE(t_mat),      INTENT(INOUT)  :: gmat
      TYPE(t_mat),      INTENT(IN)     :: vmat
      REAL,             INTENT(IN)     :: mu
      LOGICAL,          INTENT(IN)     :: l_upper !Are we in the upper half of the complex plane

      INTEGER i,j

      !WRITE(*,*) "GDFT"
      !WRITE(*,"(14f10.5)") REAL(gmat%data_c)
      !WRITE(*,"(14f10.5)") AIMAG(gmat%data_c)
      CALL gmat%inverse()
      DO i = 1, gmat%matsize1
         DO j = 1, gmat%matsize1
            IF(l_upper) THEN
               gmat%data_c(i,j) = gmat%data_c(i,j) - vmat%data_c(i,j)
            ELSE
               gmat%data_c(i,j) = gmat%data_c(i,j) - conjg(vmat%data_c(i,j))
            ENDIF
            IF(i.EQ.j) gmat%data_c(i,i) = gmat%data_c(i,i) - mu
         ENDDO
      ENDDO
      CALL gmat%inverse()
    
   END SUBROUTINE add_pot

   SUBROUTINE n_mmp_dist(n_mmp_in,n_mmp_out,hub1,results,jspins)

      USE m_types
      USE m_constants

      TYPE(t_hub1ham),     INTENT(IN)     :: hub1
      COMPLEX,             INTENT(IN)     :: n_mmp_in(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,jspins)
      COMPLEX,             INTENT(IN)     :: n_mmp_out(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,jspins)
      TYPE(t_results),     INTENT(INOUT)  :: results
      INTEGER,             INTENT(IN)     :: jspins

      INTEGER ispin,i_hia,j,k
      REAL n_in,n_out
      
      !Calculates the distance for two density matrices (maximum distance between two elements)
      n_out = 0.0
      n_in = 0.0
      results%last_mmpMatdistance = 0.0
      DO ispin = 1, jspins
         DO i_hia = 1, hub1%n_hia
            DO j = -3,3
               DO k = -3,3
                  IF((ABS(n_mmp_out(k,j,i_hia,ispin) - n_mmp_in(k,j,i_hia,ispin))).GT.results%last_mmpMatdistance) THEN
                     results%last_mmpMatdistance = ABS(n_mmp_out(k,j,i_hia,ispin) - n_mmp_in(k,j,i_hia,ispin))
                  ENDIF
                  IF(j.EQ.k) THEN
                     n_out = n_out + REAL(n_mmp_out(k,j,i_hia,ispin))
                     n_in = n_in + REAL(n_mmp_in(k,j,i_hia,ispin))
                  ENDIF
               END DO
            END DO
         END DO
      ENDDO
      results%last_occdistance = ABS(n_out-n_in)
      WRITE(*,*) "Occupation distance: ", results%last_occdistance
      WRITE(*,*) "Density matrix distance: ", results%last_mmpMatdistance
   
   END SUBROUTINE n_mmp_dist

END MODULE m_hubbard1_setup