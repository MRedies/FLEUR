MODULE m_add_selfen

   LOGICAL, PARAMETER :: l_debug = .TRUE.

   CONTAINS

   SUBROUTINE add_selfen(g,gp,selfen,atoms,noco,hub1,sym,input,ef,n_occ,mu_dc,vmmp,mmpMat)

      !Calculates the interacting Green's function for the mt-sphere with
      !
      ! (G)^-1 = (G_0)^-1 - mu 1 + V_dc - selfen
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

      IMPLICIT NONE

      TYPE(t_greensf),  INTENT(IN)     :: g
      TYPE(t_greensf),  INTENT(INOUT)  :: gp
      TYPE(t_hub1ham),  INTENT(IN)     :: hub1
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_input),    INTENT(IN)     :: input
      COMPLEX,          INTENT(IN)     :: selfen(atoms%n_hia,g%nz,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1))
      REAL,             INTENT(IN)     :: ef
      REAL,             INTENT(IN)     :: n_occ(atoms%n_hia,input%jspins)
      REAL,             INTENT(IN)     :: mu_dc
      COMPLEX,          INTENT(IN)     :: vmmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      COMPLEX,          INTENT(OUT)    :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)

      INTEGER i_hia,l,nType,ns,ispin,m,iz,ipm,spin_match,matsize,start,end,i_match
      CHARACTER(len=6) app,filename

      REAL mu_a,mu_b,mu_step,mu_max,n_max
      REAL mu,n,n_target
      TYPE(t_mat) :: gmat,vmat

      !Interval where we expect the correct mu
      mu_a = -4.0
      mu_b = 4.0
      mu_step = 0.05
      mu_max = 0.0
      n_max = 0.0

      spin_match = MERGE(1,2,noco%l_soc.OR.noco%l_noco)
      spin_match=1
      

      DO i_hia = 1, atoms%n_hia
         l = atoms%lda_hia(i_hia)%l
         nType = atoms%lda_hia(i_hia)%atomType
         ns = 2*l+1
         matsize = ns*MERGE(2,1,spin_match.EQ.1)
         CALL vmat%init(.false.,matsize,matsize)
         !Search for the maximum of occupation
         DO i_match = 1, spin_match
            n_target = MERGE(SUM(n_occ(i_hia,:)),n_occ(i_hia,i_match),spin_match.EQ.1)
            WRITE(filename,9000) i_match
            IF(l_debug) OPEN(unit=1337,file=TRIM(ADJUSTL(filename)),status="replace",action="write")
            mu = mu_a
            DO WHILE(mu.LE.mu_b)
               mu = mu + mu_step
               DO iz = 1, g%nz
                  start = MERGE(1,1+(i_match-1)*ns,spin_match.EQ.1)
                  end   = MERGE(2*ns,i_match*ns,spin_match.EQ.1)
                  vmat%data_c = selfen(i_hia,iz,start:end,start:end)
                  IF(.NOT.input%l_gfmperp.AND.spin_match.EQ.1) THEN
                     !Dismiss spin-off-diagonal elements
                     vmat%data_c(1:ns,ns+1:2*ns) = 0.0
                     vmat%data_c(ns+1:2*ns,1:ns) = 0.0
                  ENDIF
                  DO ipm = 1, 2
                     IF(spin_match.EQ.2) THEN
                        CALL g%get_gf(gmat,atoms,input,iz,l,nType,ipm.EQ.2,spin=i_match)
                     ELSE
                        CALL g%get_gf(gmat,atoms,input,iz,l,nType,ipm.EQ.2)
                     ENDIF
                     CALL add_pot(gmat,vmat,mu-mu_dc,(ipm.EQ.1))
                     CALL gp%set_gf(gmat,atoms,input,iz,l,nType,ipm.EQ.2)
                     CALL gmat%free()
                  ENDDO
               ENDDO
               CALL occmtx(gp,l,nType,atoms,sym,input,mmpMat(:,:,i_hia,:))
               !Calculate the trace
               n = 0.0
               DO ispin = 1, input%jspins
                  DO m = -l, l
                     n = n + mmpMat(m,m,i_hia,ispin)
                  ENDDO
               ENDDO
               IF(l_debug) WRITE(1337,"(2f15.8)") mu,n
               IF(n.GT.n_max) THEN
                  mu_max = mu
                  n_max  = n
               ENDIF
            ENDDO
            IF(l_debug) CLOSE(1337)

            !Sanity check for the maximum occupation
            IF(n_max.GT.2*ns+0.5) THEN
               !These oscillations seem to emerge when the lorentzian smoothing is done inadequately
               CALL juDFT_error("Something went wrong with the addition of the selfenergy",calledby="add_selfen")
            ENDIF

            !Set up the interval for the bisection method (mu_max,mu_b)
            mu_a = mu_max
            DO 
               mu = (mu_a + mu_b)/2.0
               DO iz = 1, g%nz
                  vmat%data_c = selfen(i_hia,iz,1:2*ns,1:2*ns)
                  IF(.NOT.input%l_gfmperp) THEN
                     !Dismiss spin-off-diagonal elements
                     vmat%data_c(1:ns,ns+1:2*ns) = 0.0
                     vmat%data_c(ns+1:2*ns,1:ns) = 0.0
                  ENDIF
                  DO ipm = 1, 2
                     CALL g%get_gf(gmat,atoms,input,iz,l,nType,ipm.EQ.2)
                     CALL add_pot(gmat,vmat,mu-mu_dc,(ipm.EQ.1))
                     CALL gp%set_gf(gmat,atoms,input,iz,l,nType,ipm.EQ.2)
                     CALL gmat%free()
                  ENDDO
               ENDDO
               CALL occmtx(gp,l,nType,atoms,sym,input,mmpMat(:,:,i_hia,:))
               !Calculate the trace
               n = 0.0
               DO ispin = 1, input%jspins
                  DO m = -l, l
                     n = n + mmpMat(m,m,i_hia,ispin)
                  ENDDO
               ENDDO
               IF(ABS(n-n_target).LT.0.001.OR.ABS((mu_b - mu_a)/2.0).LT.0.00001) THEN
                  !We found the chemical potential to within the desired accuracy
                  WRITE(6,"(A)") "Calculated mu to match Self-energy to DFT-GF"
                  WRITE(6,"(TR3,A4,f8.4)") "mu = ", mu
                  EXIT
               ELSE IF((n - n_target).GT.0) THEN
                  !The occupation is to small --> choose the left interval
                  mu_a = mu
               ELSE IF((n - n_target).LT.0) THEN
                  !The occupation is to big --> choose the right interval
                  mu_b = mu
               ENDIF
            ENDDO
         ENDDO
         CALL vmat%free()
      ENDDO

9000  FORMAT("mu_",I1)

   END SUBROUTINE add_selfen

   SUBROUTINE add_pot(gmat,vmat,mu,l_upper)

      USE m_types

      TYPE(t_mat),      INTENT(INOUT)  :: gmat
      TYPE(t_mat),      INTENT(IN)     :: vmat
      REAL,             INTENT(IN)     :: mu
      LOGICAL,          INTENT(IN)     :: l_upper !Are we in the upper half of the complex plane

      INTEGER i,j
      
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

END MODULE m_add_selfen