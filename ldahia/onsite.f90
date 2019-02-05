MODULE m_gOnsite

   !These are the subroutines for the spherically averaged greens function

CONTAINS

SUBROUTINE im_gmmpMat(atoms,sym,ispin,jspins,noccbd,tetweights,wtkpt,eig,usdus,eigVecCoeffs,gOnsite)

   !This Subroutine calculates the imaginary part of the Matrix elements G^[n \sigma]_{Lm Lm'}(E+i*sigma)
   !at the current k-Point (it is called in cdnval) inside the MT-sphere (averaged over the radial component)
   !and sums over the Brillouin-Zone using the histogram method

   !It is essentially the f-density of states in a (m,mp) matrix with an additional factor - pi

   USE m_types
   USE m_constants
   USE m_differentiate

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_sym),            INTENT(IN)     :: sym
   TYPE(t_eigVecCoeffs),   INTENT(IN)     :: eigVecCoeffs
   TYPE(t_usdus),          INTENT(IN)     :: usdus
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite

   INTEGER,                INTENT(IN)     :: ispin
   INTEGER,                INTENT(IN)     :: jspins
   INTEGER,                INTENT(IN)     :: noccbd

   REAL,                   INTENT(IN)     :: wtkpt
   REAL,                   INTENT(IN)     :: tetweights(:,:)
   REAL,                   INTENT(IN)     :: eig(noccbd)
   
   LOGICAL l_zero
   INTEGER i_hia, i, j, n, nn, natom, l, m, mp, lm, lmp, it,is, isi, ilo, ilop
   REAL fac, wk
   REAL, ALLOCATABLE :: dos_weights(:)

   COMPLEX n_tmp(3,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),nr_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
   COMPLEX n1_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const), d_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

   wk = wtkpt/gOnsite%del

   ALLOCATE(dos_weights(gOnsite%ne))

   !One HIA per atom at the moment
   DO i_hia = 1, atoms%n_hia
      l = atoms%lda_hia(i_hia)%l
      n = atoms%lda_hia(i_hia)%atomType  

      !finding the right starting index
      natom = SUM(atoms%neq(:n-1))

      DO nn = 1, atoms%neq(n)
         natom = natom +1
         fac = 1.0  /  ( sym%invarind(natom) * atoms%neq(n) )
         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(natom,l,n,ispin,wk,noccbd,i_hia,fac) &
         !$OMP SHARED(atoms,sym,eigVecCoeffs,usdus,gOnsite,eig,tetweights) &
         !$OMP PRIVATE(j,m,mp,lm,lmp,ilo,ilop,it,is,isi,l_zero) &
         !$OMP PRIVATE(n_tmp,n1_tmp,nr_tmp,d_tmp,dos_weights)

         !$OMP DO
         DO i = 1, noccbd
            l_zero = .true.
            IF(gOnsite%l_tetra) THEN
               !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
               IF(ANY(tetweights(:,i).NE.0.0)) l_zero = .false.
            ELSE
               !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
               j = NINT((eig(i)-gOnsite%e_bot)/gOnsite%del)+1
               IF( (j.LE.gOnsite%ne).AND.(j.GE.1) ) l_zero = .false.
            END IF

            IF(l_zero) CYCLE

            n_tmp(:,:,:) = cmplx(0.0,0.0)
            !
            ! contribution from states
            !
            DO m = -l, l
               lm = l*(l+1)+m
               DO mp = -l,l
                  lmp = l*(l+1)+mp
                  IF(gOnsite%nr(i_hia).EQ.1) THEN
                     n_tmp(1,m,mp) = n_tmp(1,m,mp) -  pi_const*&
                                  (conjg(eigVecCoeffs%acof(i,lmp,natom,ispin))*eigVecCoeffs%acof(i,lm,natom,ispin) +&
                                   conjg(eigVecCoeffs%bcof(i,lmp,natom,ispin))*eigVecCoeffs%bcof(i,lm,natom,ispin) *&
                                   usdus%ddn(l,n,ispin))
                  ELSE
                     n_tmp(1,m,mp) = n_tmp(1,m,mp) - pi_const * conjg(eigVecCoeffs%acof(i,lm,natom,ispin))*eigVecCoeffs%acof(i,lmp,natom,ispin)
                     n_tmp(2,m,mp) = n_tmp(2,m,mp) - pi_const * conjg(eigVecCoeffs%bcof(i,lm,natom,ispin))*eigVecCoeffs%bcof(i,lmp,natom,ispin)
                     n_tmp(3,m,mp) = n_tmp(3,m,mp) - pi_const * (conjg(eigVecCoeffs%acof(i,lm,natom,ispin))*eigVecCoeffs%bcof(i,lmp,natom,ispin)+&
                                                                conjg(eigVecCoeffs%bcof(i,lm,natom,ispin))*eigVecCoeffs%acof(i,lmp,natom,ispin)) 
                  END IF
               ENDDO
            ENDDO
            !
            ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
            !
            IF(gOnsite%nr(i_hia).EQ.1) THEN
               DO ilo = 1, atoms%nlo(n)
                  IF(atoms%llo(ilo,n).EQ.l) THEN
                     DO m = -l, l
                        lm = l*(l+1)+m
                        DO mp = -l, l
                           lmp = l*(l+1)+mp

                           n_tmp(1,m,mp) = n_tmp(1,m,mp) - pi_const *(  usdus%uulon(ilo,n,ispin) * (&
                                    conjg(eigVecCoeffs%acof(i,lmp,natom,ispin))*eigVecCoeffs%ccof(m,i,ilo,natom,ispin) +&
                                    conjg(eigVecCoeffs%ccof(mp,i,ilo,natom,ispin))*eigVecCoeffs%acof(i,lm,natom,ispin) )&
                                    + usdus%dulon(ilo,n,ispin) * (&
                                    conjg(eigVecCoeffs%bcof(i,lmp,natom,ispin))*eigVecCoeffs%ccof(m,i,ilo,natom,ispin) +&
                                    conjg(eigVecCoeffs%ccof(mp,i,ilo,natom,ispin))*eigVecCoeffs%bcof(i,lm,natom,ispin)))

                           DO ilop = 1, atoms%nlo(n)
                              IF (atoms%llo(ilop,n).EQ.l) THEN

                               n_tmp(1,m,mp) = n_tmp(1,m,mp) - pi_const * usdus%uloulopn(ilo,ilop,n,ispin) *&
                                    conjg(eigVecCoeffs%ccof(mp,i,ilop,natom,ispin)) *eigVecCoeffs%ccof(m,i,ilo,natom,ispin)

                              ENDIF
                           ENDDO

                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               !
               !  n_mmp should be rotated by D_mm' ; compare force_a21; taken from n_mat.f90
               !
               fac = 1./atoms%neq(n)
                  DO m = -l,l
                     DO mp = -l,l
                        IF(gOnsite%l_tetra) THEN
                           !We need to differentiate the weights with respect to energy (can maybe be done analytically)
                           CALL diff3(tetweights(:,i),gOnsite%del,dos_weights(:))
                           DO j = 1, gOnsite%ne
                              gOnsite%im_gmmpMat(1,j,i_hia,m,mp,ispin) = gOnsite%im_gmmpMat(1,j,i_hia,m,mp,ispin) + n_tmp(1,m,mp) * fac * dos_weights(j)
                           ENDDO
                        ELSE    
                           gOnsite%im_gmmpMat(1,j,i_hia,m,mp,ispin) = gOnsite%im_gmmpMat(1,j,i_hia,m,mp,ispin) + n_tmp(1,m,mp) * fac * wk
                        END IF
                     ENDDO
                  ENDDO
            ELSE
               fac = 1./atoms%neq(n)
               DO m = -l,l
                  DO mp = -l,l
                     IF(gOnsite%l_tetra) THEN
                        !We need to differentiate the weights with respect to energy (can maybe be done analytically)
                        CALL diff3(tetweights(:,i),gOnsite%del,dos_weights(:))
                        DO j = 1, gOnsite%ne
                           gOnsite%uu(j,i_hia,m,mp,ispin) = gOnsite%uu(j,i_hia,m,mp,ispin) + n_tmp(1,m,mp) * fac * dos_weights(j)      
                           gOnsite%dd(j,i_hia,m,mp,ispin) = gOnsite%dd(j,i_hia,m,mp,ispin) + n_tmp(2,m,mp) * fac * dos_weights(j)
                           gOnsite%du(j,i_hia,m,mp,ispin) = gOnsite%du(j,i_hia,m,mp,ispin) + n_tmp(3,m,mp) * fac * dos_weights(j)
                        ENDDO
                     ELSE
                        gOnsite%uu(j,i_hia,m,mp,ispin) = gOnsite%uu(j,i_hia,m,mp,ispin) + n_tmp(1,m,mp) * fac * wk
                        gOnsite%dd(j,i_hia,m,mp,ispin) = gOnsite%dd(j,i_hia,m,mp,ispin) + n_tmp(2,m,mp) * fac * wk
                        gOnsite%du(j,i_hia,m,mp,ispin) = gOnsite%du(j,i_hia,m,mp,ispin) + n_tmp(3,m,mp) * fac * wk
                     END IF
                  ENDDO
               ENDDO
            ENDIF       
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO
   ENDDO

END SUBROUTINE im_gmmpMat

SUBROUTINE calc_onsite(atoms,enpara,vr,jspin,jspins,gOnsite,ef,sym)

   USE m_types
   USE m_constants
   USE m_juDFT
   USE m_smooth
   USE m_kkintgr
   USE m_radfun

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_enpara),         INTENT(IN)     :: enpara
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite
   TYPE(t_sym),            INTENT(IN)     :: sym

   INTEGER,                INTENT(IN)     :: jspin,jspins

   REAL,                   INTENT(IN)     :: ef
   REAL,                   INTENT(IN)     :: vr(atoms%jmtd,atoms%ntype,jspins)

   TYPE(t_usdus) usdus
   INTEGER i, i_hia, l, m, mp, jr, noded, nodeu, n, j
   REAL wronk
   CHARACTER(len=30) :: filename

   REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:)
   REAL, ALLOCATABLE :: e(:), dos(:)

   COMPLEX, ALLOCATABLE :: mmpMat(:,:,:,:)
   ALLOCATE(mmpMat(atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins))

   ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jspins) )
   ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jspins) )

   CALL usdus%init(atoms,jspins)

   ALLOCATE(dos(gOnsite%ne))
   dos = 0.0

   IF(gOnsite%sigma.NE.0.0) THEN
      !construct an energy grid for smoothing (required by the function in m_smooth)
      ALLOCATE (e(gOnsite%ne))
      DO i = 1, gOnsite%ne
         e(i) = gOnsite%del * (i-1) + gOnsite%e_bot
      ENDDO
   END IF

   DO i_hia = 1, atoms%n_hia
      l = atoms%lda_hia(i_hia)%l
      n = atoms%lda_hia(i_hia)%atomType  


      IF(gOnsite%nr(i_hia).NE.1)CALL radfun(l,n,jspin,enpara%el0(l,n,jspin),vr(:,n,jspin),atoms,&
                                       f(:,:,l,jspin),g(:,:,l,jspin),usdus, nodeu,noded,wronk)

      DO m = -l, l
         DO mp = -l,l
            !calculate the imaginary part if we use the tetrahedron method
            IF(gOnsite%l_tetra) THEN
               IF(gOnsite%nr(i_hia).NE.1) THEN
                  DO jr = 1, gOnsite%nr(i_hia)
                     gOnsite%im_gmmpmat(jr,:,i_hia,m,mp,jspin) = gOnsite%im_gmmpmat(jr,:,i_hia,m,mp,jspin) + &
                                       gOnsite%uu(:,i_hia,m,mp,jspin) * (f(jr,1,l,jspin)*f(jr,1,l,jspin)+f(jr,2,l,jspin)*f(jr,2,l,jspin)) +&
                                       gOnsite%dd(:,i_hia,m,mp,jspin) * (g(jr,1,l,jspin)*g(jr,1,l,jspin)+g(jr,2,l,jspin)*g(jr,2,l,jspin)) +&
                                       gOnsite%du(:,i_hia,m,mp,jspin) * (f(jr,1,l,jspin)*g(jr,1,l,jspin)+f(jr,2,l,jspin)*g(jr,2,l,jspin))
                  ENDDO
               END IF
            ENDIF
            !
            !taking care of spin degeneracy
            !
            IF(jspins.EQ.1) gOnsite%im_gmmpMat(:,:,i_hia,m,mp,1) = 2.0 * gOnsite%im_gmmpMat(:,:,i_hia,m,mp,1)
            !
            !smooth the imaginary part using gaussian broadening 
            !
            IF(gOnsite%sigma.NE.0.0) THEN
               DO jr = 1, gOnsite%nr(i_hia)
                  CALL smooth(e(:),gOnsite%im_gmmpMat(jr,:,i_hia,m,mp,jspin),gOnsite%sigma,gOnsite%ne)
               ENDDO
            ENDIF

         ENDDO
      ENDDO
      !
      !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration 
      !
      CALL greensf_cutoff(gOnsite,atoms,jspins)

      CALL gOnsite%init_e_contour(ef)

      CALL timestart("On-Site: Kramer-Kronigs-Integration")
      DO m= -l,l
         DO mp= -l,l
            DO jr = 1, gOnsite%nr(i_hia)
               IF(gOnsite%mode.EQ.1) THEN
                  CALL kkintgr_real(gOnsite%nz,gOnsite%e(:),gOnsite%ne,gOnsite%sigma,gOnsite%del,gOnsite%e_bot,&
                                    gOnsite%im_gmmpMat(jr,:,i_hia,m,mp,jspin),gOnsite%gmmpMat(jr,:,i_hia,m,mp,jspin))
               ELSE IF(gOnsite%mode.EQ.2) THEN
                  CALL kkintgr_complex(gOnsite%nz,gOnsite%e(:),gOnsite%ne,gOnsite%sigma,gOnsite%del,gOnsite%e_bot,&
                                       gOnsite%im_gmmpMat(jr,:,i_hia,m,mp,jspin),gOnsite%gmmpMat(jr,:,i_hia,m,mp,jspin))
               END IF
            ENDDO
         ENDDO
      ENDDO
      CALL timestop("On-Site: Kramer-Kronigs-Integration")

      CALL gOnsite%calc_mmpmat(atoms,sym,jspins,mmpMat)

      !write density matrix to file 

      filename = "n_mmp_mat_g"

      OPEN (69,file=TRIM(ADJUSTL(filename)),status='replace',form='formatted')
      WRITE (69,'(14f14.8)') mmpMat(:,:,:,:)
      CLOSE (69)

   ENDDO


END SUBROUTINE calc_onsite


SUBROUTINE greensf_cutoff(gOnsite,atoms,jspins)
   !This Subroutine determines the cutoff energy for the kramers-kronig-integration
   !This cutoff energy is defined so that the integral over the fDOS up to this cutoff 
   !is equal to 14 (the number of states in the 4f shell) or not to small


   USE m_types
   USE m_intgr
   USE m_juDFT
   USE m_constants
   USE m_kkintgr
   
   IMPLICIT NONE

   TYPE(t_greensf),     INTENT(INOUT)  :: gOnsite
   TYPE(t_atoms),       INTENT(IN)     :: atoms
   INTEGER,             INTENT(IN)     :: jspins


   REAL, ALLOCATABLE :: fDOS(:)

   INTEGER i_hia, i,l, m,mp, ispin, j, n_c, kkintgr_cut, jr, n

   REAL integral

   REAL a,b, imag
   LOGICAL l_write

   ALLOCATE(fDOS(gOnsite%ne))

   l_write=.true.

   DO i_hia = 1, atoms%n_hia

      fDOS(:) = 0.0
      l = atoms%lda_hia(i_hia)%l
      n = atoms%lda_hia(i_hia)%atomType  

      !Calculate the trace over m,mp of the Greens-function matrix to obtain the fDOS 

      !n_f(e) = -1/pi * TR[Im(G_f(e))]
      DO ispin = 1, jspins
         DO m = -l , l
            DO j = 1, gOnsite%ne
               IF(gOnsite%nr(i_hia).NE.1) THEN
                  CALL intgr3(gOnsite%im_gmmpMat(:,j,i_hia,m,m,ispin),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imag)
               ELSE
                  imag = gOnsite%im_gmmpMat(1,j,i_hia,m,m,ispin)
               END IF
               fDOS(j) = fDOS(j) + imag
            ENDDO
         ENDDO
      ENDDO

      fDOS(:) = -1/pi_const * fDOS(:)


      IF(l_write) THEN

         open(unit=1337,file="fDOS.txt",status="replace", action="write")

         DO j = 1, gOnsite%ne

            WRITE(1337,*) (j-1) * gOnsite%del + gOnsite%e_bot, fDOS(j)

         ENDDO

         close(unit=1337)

      END IF       



      CALL trapz(fDOS(:), gOnsite%del, gOnsite%ne, integral)

      WRITE(*,*) "Integral over fDOS: ", integral

      kkintgr_cut = gOnsite%ne

      IF(integral.LT.13.5) THEN
         ! If the integral is to small we stop here to avoid problems
         CALL juDFT_error("fDOS-integral too small: make sure numbands is big enough", calledby="greensf_cutoff")
         
      ELSE IF((integral.GT.14).AND.((integral-14).GT.0.001)) THEN
         !IF the integral is bigger than 14, search for the cutoff using the bisection method   

         a = gOnsite%e_bot
         b = gOnsite%e_top


         DO

            n_c = INT(((a+b)/2.0-gOnsite%e_bot)/gOnsite%del)+1
            
            CALL trapz(fDOS(:),gOnsite%del,n_c,integral)

            IF((ABS(integral-14).LT.0.001).OR.(ABS(a-b)/2.0.LT.gonsite%del)) THEN

               kkintgr_cut = INT(((a+b)/2.0-gOnsite%e_bot)/gOnsite%del)+1
               EXIT

            ELSE IF((integral-14).LT.0) THEN
               a = (a+b)/2.0
            ELSE IF((integral-14).GT.0) THEN
               b = (a+b)/2.0
            END IF

         ENDDO

         CALL trapz(fDOS(:),gOnsite%del,kkintgr_cut,integral)

         WRITE(*,*) "CALCULATED CUTOFF: ", kkintgr_cut
         WRITE(*,*) "INTEGRAL OVER fDOS with cutoff: ", integral
      END IF


      !Now we set the imaginary part of the greens function to zero above this cutoff
      DO ispin = 1, jspins
         DO i = kkintgr_cut+1, gOnsite%ne
            DO m = -l, l
               DO mp = -l,l
                  DO jr = 1, gOnsite%nr(i_hia)
                     gOnsite%im_gmmpMat(jr,i,i_hia,m,mp,ispin) = 0.0e0
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   ENDDO



END SUBROUTINE greensf_cutoff


END MODULE m_gOnsite