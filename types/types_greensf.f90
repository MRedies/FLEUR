MODULE m_types_greensf

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_types_greensf
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  Contains a type for onsite and intersite green's functions in the mt-sphere
   !>  It stores the energy contour in the complex plane and the corresponding   
   !>  matrix elements of the green's function
   !>  We have the following cases
   !>    -onsite
   !>       -we look at l=l' but m\=m'
   !>       -we treat non-magnetic/collinear and noco (not tested)
   !>       -we look at r=r' and spherically averaged gf
   !>    -intersite
   !>       -l\=l' and m\=m'
   !>       -r\=r' (not stored we calculate the gf by calling calc_intersite in m_intersite for specific r and r')
   !------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

      TYPE t_greensf

         !Energy contour parameters
         INTEGER  :: mode  !Determines the shape of the contour (more information in kkintgr.f90)
         INTEGER  :: nz    !number of points in the contour
         INTEGER  :: nmatsub

         !array for energy contour
         COMPLEX, ALLOCATABLE  :: e(:)  !energy points
         COMPLEX, ALLOCATABLE  :: de(:) !weights for integration

         !Arrays for Green's function
         COMPLEX, ALLOCATABLE :: gmmpMat(:,:,:,:,:,:) 

         !for radial dependence
         COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: du(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:,:)

         

         CONTAINS
            PROCEDURE, PASS :: init => greensf_init
            PROCEDURE       :: e_contour
            PROCEDURE       :: get_gf
            PROCEDURE       :: set_gf
      END TYPE t_greensf


   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE greensf_init(thisGREENSF,input,lmax,atoms,noco,nz_in,e_in,de_in,matsub_in)

         USE m_juDFT
         USE m_types_setup

         CLASS(t_greensf),    INTENT(INOUT)  :: thisGREENSF
         TYPE(t_atoms),       INTENT(IN)     :: atoms
         TYPE(t_input),       INTENT(IN)     :: input
         INTEGER,             INTENT(IN)     :: lmax
         TYPE(t_noco),        INTENT(IN)     :: noco
         !Pass a already calculated energy contour to the type 
         INTEGER, OPTIONAL,   INTENT(IN)     :: nz_in
         INTEGER, OPTIONAL,   INTENT(IN)     :: matsub_in
         COMPLEX, OPTIONAL,   INTENT(IN)     :: e_in(:)
         COMPLEX, OPTIONAL,   INTENT(IN)     :: de_in(:)

         INTEGER i,j,r_dim,l_dim,spin_dim
         REAL    tol,n
         LOGICAL l_new


         !IF(noco%l_mperp) CALL juDFT_error("NOCO + gf not implemented",calledby="greensf_init")

         !
         !Set up general parameters for the Green's function
         !
         !
         !Setting up parameters for the energy contour
         !
         thisGREENSF%mode     = input%gf_mode
         IF(PRESENT(nz_in)) THEN
            thisGREENSF%nz = nz_in 
            thisGREENSF%nmatsub = matsub_in
         ELSE
            !Parameters for the energy contour in the complex plane

            IF(thisGREENSF%mode.EQ.1) THEN
                  thisGREENSF%nz = input%gf_n1+input%gf_n2+input%gf_n3+input%gf_nmatsub
                  thisGREENSF%nmatsub = input%gf_nmatsub
            ELSE IF(thisGREENSF%mode.GE.2) THEN
               thisGREENSF%nz = input%gf_n
               thisGREENSF%nmatsub = 0
            ENDIF
         END IF

         ALLOCATE (thisGREENSF%e(thisGREENSF%nz))
         ALLOCATE (thisGREENSF%de(thisGREENSF%nz))

         IF(PRESENT(e_in)) THEN
            thisGREENSF%e(:) = e_in(:)
            thisGREENSF%de(:)= de_in(:)
         ELSE
            !If no energy contour is given it is set up to zero
            thisGREENSF%e(:) = CMPLX(0.0,0.0)
            thisGREENSF%de(:)= CMPLX(0.0,0.0)
         END IF


         IF(atoms%n_gf.GT.0) THEN !Are there Green's functions to be calculated?
            IF(input%l_gfsphavg) THEN
               spin_dim = MERGE(4,input%jspins,input%l_gfmperp)
               ALLOCATE (thisGREENSF%gmmpMat(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim,2))
               thisGREENSF%gmmpMat = 0.0
            ELSE
               ALLOCATE (thisGREENSF%uu(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim,2))
               ALLOCATE (thisGREENSF%dd(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim,2))
               ALLOCATE (thisGREENSF%du(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim,2))
               ALLOCATE (thisGREENSF%ud(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim,2))
               thisGREENSF%uu = 0.0
               thisGREENSF%dd = 0.0
               thisGREENSF%du = 0.0
               thisGREENSF%ud = 0.0
            ENDIF
         ENDIF

      END SUBROUTINE greensf_init

      SUBROUTINE e_contour(this,input,eb,et,ef)


         USE m_types_setup
         USE m_constants
         USE m_juDFT
         USE m_grule
         USE m_ExpSave

         IMPLICIT NONE

         CLASS(t_greensf),  INTENT(INOUT)  :: this
         TYPE(t_input),     INTENT(IN)     :: input
         REAL,              INTENT(IN)     :: eb  
         REAL,              INTENT(IN)     :: et
         REAL,              INTENT(IN)     :: ef
         
         INTEGER i, j, iz,nz

         REAL e1, e2, del, sigma
         COMPLEX de
         REAL r, xr, xm, c, s, a, b
         REAL x(this%nz), w(this%nz)



         IF(this%mode.EQ.1) THEN
            
            sigma = input%gf_sigma * pi_const

            IF(this%nmatsub > 0) THEN

               nz = 0
               !Left Vertical part (eb,0) -> (eb,sigma)
               de = this%nmatsub * CMPLX(0.0,sigma)
               CALL grule(input%gf_n1,x(1:(input%gf_n1)/2),w(1:(input%gf_n1)/2))
               x = -x
               DO i = 1, (input%gf_n1+3)/2-1
                  x(input%gf_n1-i+1) = -x(i)
                  w(input%gf_n1-i+1) =  w(i)
               ENDDO
               DO iz = 1, input%gf_n1
                  nz = nz + 1
                  IF(nz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
                  this%e(nz) = eb + de + de * x(iz) 
                  this%de(nz) = w(iz)*de
               ENDDO

               !Horizontal Part (eb,sigma) -> (et,sigma)
               de = (ef-30*input%gf_sigma-eb)/2.0
               CALL grule(input%gf_n2,x(1:(input%gf_n2)/2),w(1:(input%gf_n2)/2))
               x = -x
               DO i = 1, (input%gf_n2+3)/2-1
                  x(input%gf_n2-i+1) = -x(i)
                  w(input%gf_n2-i+1) =  w(i)
               ENDDO
               DO iz = 1, input%gf_n2
                  nz = nz + 1
                  IF(nz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
                  this%e(nz) = de*x(iz) + de + eb + 2 * this%nmatsub * ImagUnit * sigma
                  this%de(nz) = de*w(iz)
               ENDDO

               !Right Vertical part (et,sigma) -> infty
               CALL grule(input%gf_n3,x(1:(input%gf_n3)/2),w(1:(input%gf_n3)/2))
               x = -x
               DO i = 1, (input%gf_n3+3)/2-1
                  x(input%gf_n3-i+1) = -x(i)
                  w(input%gf_n3-i+1) =  w(i)
               ENDDO
               de = 30*input%tkb
               DO iz = 1, input%gf_n3
                  nz = nz + 1
                  IF(nz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
                  this%e(nz) = de*x(iz)+ef +  2 * this%nmatsub * ImagUnit * sigma
                  this%de(nz) = w(iz)*de/(1.0+exp_save((REAL(this%e(nz))-ef)/input%gf_sigma))
               ENDDO

               !Matsubara frequencies
               DO iz = this%nmatsub , 1, -1 
                  nz = nz + 1 
                  IF(nz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
                  this%e(nz)  = ef + (2*iz-1) * ImagUnit *sigma 
                  this%de(nz) =  -2 *ImagUnit * sigma
               ENDDO 
               WRITE(*,1000) this%nz, this%nmatsub,input%gf_n1,input%gf_n2,input%gf_n3
            ENDIF
1000     FORMAT("Energy Contour for Green's Function Integration with nz: ", I5.1,"; nmatsub: ", I5.1,"; n1: ", I5.1,"; n2: ", I5.1,"; n3: ", I5.1)
         ELSE IF(this%mode.EQ.2) THEN

            !Semicircle
            e1 = eb
            e2 = ef

            this%nmatsub = 0
            !Radius
            r  = (e2-e1)*0.5
            !midpoint
            xr = (e2+e1)*0.5

            CALL grule(this%nz,x(1:(this%nz)/2),w(1:(this%nz)/2))

            DO i = 1, this%nz/2
               x(this%nz-i+1) = -x(i)
               w(this%nz-i+1) =  w(i)
            ENDDO
            DO i = 1, this%nz
               this%e(i) = xr + ImagUnit * r * EXP(-ImagUnit*pi_const/2.0 * x(this%nz-i+1))
               this%de(i) = pi_const/2.0 * r * EXP(-ImagUnit*pi_const/2.0 * x(this%nz-i+1)) * w(this%nz-i+1)
            ENDDO

         ELSE IF(this%mode.EQ.3) THEN
            !Equidistant contour (without vertical edges)
   
            de = (et-eb)/REAL(this%nz-1)
            DO iz = 1, this%nz
               this%e(iz) = (iz-1) * de + eb + ImagUnit * input%gf_sigma
               this%de(iz) = de
            ENDDO

         ELSE

            CALL juDFT_error("Invalid mode for energy contour in Green's function calculation", calledby="init_e_contour")

         END IF

      END SUBROUTINE e_contour

      SUBROUTINE get_gf(this,gmat,atoms,input,iz,l,nType,l_conjg,spin,lp,nTypep,u,udot)

         USE m_types_mat
         USE m_types_setup
         USE m_juDFT
         USE m_ind_greensf

         !Returns the matrix belonging to energy point iz with l,lp,nType,nTypep
         !when jr (and jrp) are given return for that radial point

         IMPLICIT NONE

         CLASS(t_greensf),    INTENT(IN)  :: this
         TYPE(t_atoms),       INTENT(IN)  :: atoms 
         TYPE(t_input),       INTENT(IN)  :: input
         TYPE(t_mat),         INTENT(OUT) :: gmat !Return matrix 

         INTEGER,             INTENT(IN)  :: iz
         INTEGER,             INTENT(IN)  :: nType
         INTEGER,             INTENT(IN)  :: l 
         LOGICAL,             INTENT(IN)  :: l_conjg
         INTEGER, OPTIONAL,   INTENT(IN)  :: spin
         INTEGER, OPTIONAL,   INTENT(IN)  :: nTypep
         INTEGER, OPTIONAL,   INTENT(IN)  :: lp 
         REAL   , OPTIONAL,   INTENT(IN)  :: u(2,input%jspins)       !Radial functions at the point where you want to evaluate the greens function 
         REAL   , OPTIONAL,   INTENT(IN)  :: udot(2,input%jspins)

         INTEGER matsize1,matsize2,i_gf,i,j,ind1,ind2,ind1_start,ind2_start
         INTEGER m,mp,spin1,spin2,ipm,ispin,ispin_end
         INTEGER lp_loop
         LOGICAL l_radial

         IF(PRESENT(u).OR.PRESENT(udot).AND.input%l_gfsphavg) THEN
            CALL juDFT_error("Greens function not calculated for radial dependence", calledby="get_gf")
         ENDIF

         IF((PRESENT(u).AND..NOT.PRESENT(udot)).OR.&
            (PRESENT(udot).AND..NOT.PRESENT(u))) THEN
            CALL juDFT_error("Not a valid input: Either provide both u and udot or neither of them", calledby="get_gf")
         ENDIF

         l_radial = PRESENT(u).AND.PRESENT(udot)

         IF(PRESENT(spin)) THEN
            IF(spin.GT.4.OR.spin.LT.1) THEN
               CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
            ENDIF
         ENDIF
         
         !Determine matsize for the result gmat (if spin is given only return this digonal element)
         matsize1 = (2*l+1) * MERGE(1,input%jspins,PRESENT(spin))
         IF(PRESENT(lp)) THEN
            matsize2 = (2*lp+1) * MERGE(1,input%jspins,PRESENT(spin))
         ELSE
            matsize2 = matsize1
         ENDIF 
         CALL gmat%init(.FALSE.,matsize1,matsize2)

         IF(.NOT.PRESENT(lp)) THEN
            lp_loop = l 
         ELSE 
            lp_loop = lp
         ENDIF

         !Find the index i_gf corresponding to l,lp,nType,nTypep
         i_gf = ind_greensf(atoms,l,nType,lp,nTypep)
         ipm = MERGE(2,1,l_conjg)

         gmat%data_c = 0.0
         ispin_end = MERGE(4,input%jspins,input%l_gfmperp)

         DO ispin = MERGE(spin,1,PRESENT(spin)), MERGE(spin,ispin_end,PRESENT(spin))
            !Find the right quadrant in gmat according to the spin index
            spin1 = MERGE(ispin,1,ispin.NE.3)
            spin1 = MERGE(spin1,2,ispin.NE.4)
            spin2 = MERGE(ispin,2,ispin.NE.3)
            spin2 = MERGE(spin2,1,ispin.NE.4)
            
            IF(.NOT.PRESENT(spin)) THEN
               ind1_start = (spin1-1)*(2*l+1) 
               ind2_start = (spin2-1)*(2*lp_loop+1) 
            ELSE 
               ind1_start = 0
               ind2_start = 0
            ENDIF
            ind1 = ind1_start
            DO m = -l,l 
               ind1 = ind1 + 1 
               ind2 = ind2_start
               DO mp = -lp_loop,lp_loop
                  ind2 = ind2 + 1
                  IF(l_radial) THEN
                     gmat%data_c(ind1,ind2) = this%uu(iz,i_gf,m,mp,ispin,ipm) * u(1,spin1)*u(2,spin2) + &
                                              this%dd(iz,i_gf,m,mp,ispin,ipm) * udot(1,spin1)*udot(2,spin2) + &
                                              this%du(iz,i_gf,m,mp,ispin,ipm) * udot(1,spin1)*u(2,spin2) + &
                                              this%ud(iz,i_gf,m,mp,ispin,ipm) * u(1,spin1)*udot(2,spin2) 
                  ELSE        
                     gmat%data_c(ind1,ind2) = this%gmmpMat(iz,i_gf,m,mp,ispin,ipm)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

      END SUBROUTINE get_gf

      SUBROUTINE set_gf(this,gmat,atoms,input,iz,l,nType,l_conjg,spin,lp,nTypep) 

         USE m_types_mat
         USE m_types_setup
         USE m_juDFT
         USE m_ind_greensf

         !Sets the spherically averaged greens function matrix belonging to energy point iz with l,lp,nType,nTypep
         !equal to gmat 

         IMPLICIT NONE

         CLASS(t_greensf),    INTENT(INOUT)  :: this
         TYPE(t_atoms),       INTENT(IN)     :: atoms 
         TYPE(t_input),       INTENT(IN)     :: input
         TYPE(t_mat),         INTENT(IN)     :: gmat  

         INTEGER,             INTENT(IN)     :: iz
         INTEGER,             INTENT(IN)     :: nType
         INTEGER,             INTENT(IN)     :: l 
         LOGICAL,             INTENT(IN)     :: l_conjg
         INTEGER, OPTIONAL,   INTENT(IN)     :: spin
         INTEGER, OPTIONAL,   INTENT(IN)     :: nTypep
         INTEGER, OPTIONAL,   INTENT(IN)     :: lp 

         INTEGER matsize1,matsize2,i_gf,i,j,ind1,ind2,ind1_start,ind2_start
         INTEGER m,mp,spin1,spin2,ipm,ispin,ispin_end
         INTEGER lp_loop

         IF(PRESENT(spin)) THEN
            IF(spin.GT.4.OR.spin.LT.1) THEN
               CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
            ENDIF
         ENDIF
         
         !Determine matsize for the result gmat (if spin is given only return this digonal element)
         matsize1 = (2*l+1) * MERGE(1,input%jspins,PRESENT(spin))
         IF(PRESENT(lp)) THEN
            matsize2 = (2*lp+1) * MERGE(1,input%jspins,PRESENT(spin))
         ELSE
            matsize2 = matsize1
         ENDIF 
         
         !Check the expected matsizes against the actual
         IF(matsize1.NE.gmat%matsize1.OR.matsize2.NE.gmat%matsize2) THEN
            CALL juDFT_error("Mismatch in matsizes", calledby="set_gf")
         ENDIF

         IF(.NOT.PRESENT(lp)) THEN
            lp_loop = l 
         ELSE 
            lp_loop = lp
         ENDIF

         !Find the index i_gf corresponding to l,lp,nType,nTypep
         i_gf = ind_greensf(atoms,l,nType,lp,nTypep)
         ipm = MERGE(2,1,l_conjg)

         ispin_end = MERGE(4,input%jspins,input%l_gfmperp)

         DO ispin = MERGE(spin,1,PRESENT(spin)), MERGE(spin,ispin_end,PRESENT(spin))
            !Find the right quadrant in gmat according to the spin index
            IF(.NOT.PRESENT(spin)) THEN
               spin1 = MERGE(ispin,1,ispin.NE.3)
               spin1 = MERGE(spin1,2,ispin.NE.4)
               spin2 = MERGE(ispin,2,ispin.NE.3)
               spin2 = MERGE(spin2,1,ispin.NE.4)
               ind1_start = (spin1-1)*(2*l+1) 
               ind2_start = (spin2-1)*(2*lp_loop+1) 
            ELSE 
               ind1_start = 0
               ind2_start = 0
            ENDIF

            ind1 = ind1_start
            DO m = -l,l 
               ind1 = ind1 + 1 
               ind2 = ind2_start
               DO mp = -lp_loop,lp_loop
                  ind2 = ind2 + 1
                  this%gmmpMat(iz,i_gf,m,mp,ispin,ipm) = gmat%data_c(ind1,ind2) 
               ENDDO
            ENDDO
         ENDDO

      END SUBROUTINE set_gf

END MODULE m_types_greensf