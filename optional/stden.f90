      MODULE m_stden
      USE m_juDFT
!     ************************************************************
!     generate flapw starting density by superposition of
!     atomic densities. the non-spherical terms inside
!     the spheres are obtained by a least squares fit
!     and the interstitial and vacuum warping terms are calculated
!     by a fast fourier transform.
!     e. wimmer   nov. 1984       c.l.fu diagonized 1987
!     *************************************************************
      CONTAINS
        SUBROUTINE stden(mpi,&
             &                 sphhar,stars,atoms,sym,&
             &                 DIMENSION,vacuum,&
             &                 input,&
             &                 cell,&
             &                 xcpot,&
             &                 obsolete,&
             &                 oneD)
          USE m_sphpts
          USE m_constants
          USE m_enpara,    ONLY : w_enpara
          USE m_xcall,     ONLY : vxcall
          USE m_qsf
          USE m_checkdop
          USE m_cdnovlp
          USE m_cdn_io
          USE m_qfix
          USE m_atom2
          USE m_types
          USE m_cylpts
          USE m_points
          USE m_juDFT_init
          IMPLICIT NONE
          !     ..
          TYPE(t_mpi),INTENT(IN)      :: mpi
          TYPE(t_atoms),INTENT(IN)    :: atoms
          TYPE(t_dimension),INTENT(IN):: DIMENSION
          TYPE(t_sphhar),INTENT(IN)   :: sphhar
          TYPE(t_obsolete),INTENT(IN) :: obsolete
          TYPE(t_sym),INTENT(IN)      :: sym
          TYPE(t_stars),INTENT(IN)    :: stars
          TYPE(t_oneD),INTENT(IN)     :: oneD
          TYPE(t_input),INTENT(IN)    :: input
          TYPE(t_vacuum),INTENT(IN)   :: vacuum
          TYPE(t_cell),INTENT(IN)     :: cell
          TYPE(t_xcpot),INTENT(IN)    :: xcpot
          

          TYPE(t_enpara)   :: enpara
          CHARACTER*8 :: name(10)
          !-odim
          !+odim
          !     ..
          !     .. Local Scalars ..
          REAL d,del,fix,h,r,rnot,sign,z,bm,qdel
          REAL denz1(1),vacxpot(1),vacpot(1) 
          INTEGER i,iter,ivac,iza,j,jr,k,n,n1,npd,ispin 
          INTEGER nw,ilo,natot,icorr_dummy,nat 
          COMPLEX czero
          COMPLEX :: cdom(1),cdomvz(1,1),cdomvxy(1,1,1)
          !     ..
          !     .. Local Arrays ..
          COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
          REAL,    ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:),vbar(:,:)
          REAL,    ALLOCATABLE :: xp(:,:),rat(:,:),eig(:,:,:),sigm(:)
          REAL,    ALLOCATABLE :: rh(:,:,:),rh1(:,:,:),rhoss(:,:)
          REAL,    ALLOCATABLE :: vacpar(:)
          INTEGER lnum(DIMENSION%nstd,atoms%ntype),nst(atoms%ntype) 
          INTEGER jrc(atoms%ntype)
          LOGICAL l_found(0:3),llo_found(atoms%nlod),l_enpara,l_st
          CHARACTER*8 name_l(10)
          !     ..
          !     .. Intrinsic Functions ..
          INTRINSIC abs,exp
          !     ..
          !     .. Data statements ..
          DATA del/1.e-6/
          PARAMETER (czero=(0.0,0.0),l_st=.true.)
          !     ..
          !
          IF (input%jspins > DIMENSION%jspd)  CALL juDFT_error("input%jspins > dimension%jspd",calledby&
               &     ="stden")
          ALLOCATE ( rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,DIMENSION%jspd) )

          ALLOCATE ( qpw(stars%ng3,DIMENSION%jspd),rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,DIMENSION%jspd) )
          ALLOCATE ( xp(3,DIMENSION%nspd),rat(DIMENSION%msh,atoms%ntype),eig(DIMENSION%nstd,DIMENSION%jspd,atoms%ntype) )
          ALLOCATE ( rh(DIMENSION%msh,atoms%ntype,DIMENSION%jspd),rh1(DIMENSION%msh,atoms%ntype,DIMENSION%jspd) )
          ALLOCATE ( enpara%ello0(atoms%nlod,atoms%ntype,input%jspins),vacpar(2) )
          ALLOCATE ( enpara%el0(0:3,atoms%ntype,input%jspins))   
          ALLOCATE ( enpara%lchange(0:3,atoms%ntype,input%jspins))
          ALLOCATE ( enpara%skiplo(atoms%ntype,input%jspins))
          ALLOCATE ( enpara%llochg(atoms%nlod,atoms%ntype,input%jspins))
          ALLOCATE ( enpara%enmix(input%jspins))
          ALLOCATE ( enpara%evac0(2,dimension%jspd))
          ALLOCATE ( enpara%lchg_v(2,dimension%jspd))
          ALLOCATE ( rht(vacuum%nmzd,2,DIMENSION%jspd),vbar(2,atoms%ntype),sigm(vacuum%nmzd) )
          ALLOCATE ( rhoss(DIMENSION%msh,DIMENSION%jspd) )
          enpara%enmix=1.0
          rho = 0.0
          IF (mpi%irank == 0) THEN

             !--->    if sigma is not 0.0, then divide this charge among all atoms
             IF ( ABS(input%efield%sigma).LT. 1.e-6) THEN
                qdel = 0.0
             ELSE
                natot = 0
                DO n=1,atoms%ntype
                   IF (atoms%zatom(n).GE.1.0) natot = natot + atoms%neq(n)
                ENDDO
                qdel = 2.*input%efield%sigma/natot
             ENDIF
             !
             WRITE (6,FMT=8000)
8000         FORMAT (/,/,/,' superposition of atomic densities',/,/,&
                  &       ' original atomic densities:',/)
             DO n = 1,atoms%ntype
                r = atoms%rmsh(1,n)
                d = EXP(atoms%dx(n))
                !         DO i = 1, msh
                !            rat(i,n) = r
                !            r = r*d
                !         ENDDO
                jrc(n) = 0
                DO WHILE (r < atoms%rmt(n) + 20.0) 
                   IF ( jrc(n) > DIMENSION%msh )   CALL juDFT_error&
                        &           ("increase msh in fl7para!",calledby ="stden")
                   jrc(n) = jrc(n) + 1
                   rat(jrc(n),n) = r
                   r = r*d
                ENDDO
             ENDDO
             !
             ! Generate the atomic charge densities
             !
             DO n = 1,atoms%ntype
                z = atoms%zatom(n)
                r = atoms%rmt(n)
                h = atoms%dx(n)
                jr = atoms%jri(n)
                IF (input%jspins.EQ.2) THEN
                   bm = atoms%bmu(n)
                ELSE
                   bm = 0.
                ENDIF
                !--->    check whether this atom has been done already
                DO n1 = 1,n - 1
                   IF (ABS(z-atoms%zatom(n1)).GT.del) GO TO 40
                   IF (ABS(r-atoms%rmt(n1)).GT.del) GO TO 40
                   IF (ABS(h-atoms%dx(n1)).GT.del) GO TO 40
                   IF (ABS(bm-atoms%bmu(n1)).GT.del) GO TO 40
                   IF (jr.NE.atoms%jri(n1)) GO TO 40
                   DO ispin = 1, input%jspins
                      DO i = 1,jrc(n) ! dimension%msh
                         rh(i,n,ispin) = rh(i,n1,ispin)
                      ENDDO
                   ENDDO
                   nst(n) = nst(n1)
                   vbar(1,n) = vbar(1,n1)
                   vbar(input%jspins,n) = vbar(input%jspins,n1)
                   DO i = 1, nst(n1)
                      lnum(i,n)  = lnum(i,n1)
                      eig(i,1,n) = eig(i,1,n1)
                      eig(i,input%jspins,n) = eig(i,input%jspins,n1)
                   ENDDO
                   GO TO 70
40                 CONTINUE
                ENDDO
                !--->    new atom
                rnot = atoms%rmsh(1,n)
                IF (z.LT.1.0) THEN
                   DO ispin = 1, input%jspins
                      DO i = 1,jrc(n) ! dimension%msh
                         rh(i,n,ispin) = 1.e-10
                      ENDDO
                   ENDDO
                ELSE
                   CALL atom2(&
                        &                 DIMENSION,atoms,xcpot,input,n,jrc(n),rnot,qdel,&
                        &                 rhoss,nst(n),lnum(1,n),eig(1,1,n),vbar(1,n))
                   DO ispin = 1, input%jspins
                      DO i = 1, jrc(n) ! dimension%msh
                         rh(i,n,ispin) = rhoss(i,ispin)
                      ENDDO
                   ENDDO
                END IF
                !--->    list atomic density
                iza = atoms%zatom(n) + 0.0001
                WRITE (6,FMT=8030) namat_const(iza)
8030            FORMAT (/,/,' atom: ',a2,/)
8040            FORMAT (4 (3x,i5,f8.5,f12.6))
70              CONTINUE
             ENDDO

             !roa+
             !..use cdnovlp to generate total density out of atom densities...
             !
             DO ispin = 1, input%jspins
                nat = 1
                DO  n = 1,atoms%ntype
                   DO  i = 1, jrc(n)
                      rh1(i,n,ispin) = rh(i,n,ispin)*fpi_const*rat(i,n)**2
                   ENDDO
                   rh1(jrc(n):DIMENSION%msh,n,ispin) = 0.0
                   !..prepare spherical mt charge
                   DO i = 1,atoms%jri(n)
                      rho(i,0,n,ispin) = rh(i,n,ispin)*sfp_const*atoms%rmsh(i,n)**2
                   ENDDO
                   !..reset nonspherical mt charge
                   DO k = 1,sphhar%nlh(atoms%ntypsy(nat))
                      DO j = 1,atoms%jri(n)
                         rho(j,k,n,ispin) = 0.e0
                      ENDDO
                   ENDDO
                   nat = nat + atoms%neq(n)
                ENDDO
             ENDDO ! ispin

             qpw(:,:) = czero ; rht(:,:,:) = 0.e0 ; rhtxy(:,:,:,:) = czero

          ENDIF ! mpi%irank == 0
          DO ispin = 1, input%jspins
             CALL cdnovlp(mpi,&
                  &               sphhar,stars,atoms,sym,&
                  &               DIMENSION,vacuum,&
                  &               cell,&
                  &               input,oneD,l_st,&
                  &               ispin,rh1(:,:,ispin),&
                  &               qpw,rhtxy,rho,rht)
             !roa-
             !-spinloop
          ENDDO
          IF ( mpi%irank == 0 ) THEN
             !
             ! Check the normalization of total density
             !
             CALL qfix(&
                  &          stars,atoms,sym,vacuum,&
                  &          sphhar,input,cell,oneD,&
                  &          qpw,rhtxy,rho,rht,.FALSE.,.true.,&
                  &          fix)
             !
             ! Write superposed density onto density file
             !
             iter = 0
             CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,CDN_INPUT_DEN_const,&
                               1,-1.0,0.0,.TRUE.,iter,rho,qpw,rht,rhtxy,cdom,cdomvz,cdomvxy)
             !
             ! Check continuity
             !
             IF (input%vchk)THEN
                DO ispin = 1,input%jspins
                   WRITE (6,'(a8,i2)') 'spin No.',ispin
                   IF (input%film .AND. .NOT.oneD%odi%d1) THEN
                      !         ---> vacuum boundaries
                      npd = MIN(DIMENSION%nspd,25)
                      CALL points(xp,npd)
                      DO ivac = 1,vacuum%nvac
                         sign = 3. - 2.*ivac
                         DO j = 1,npd
                            xp(3,j) = sign*cell%z1/cell%amat(3,3)
                         ENDDO
                         CALL checkdop(&
                              &                    xp,npd,0,0,ivac,1,ispin,.TRUE.,DIMENSION,atoms,&
                              &                    sphhar,stars,sym,&
                              &                    vacuum,cell,oneD,&
                              &                    qpw,rho,rhtxy,rht)
                      ENDDO
                   ELSEIF (oneD%odi%d1) THEN
                      !-odim
                      npd = MIN(DIMENSION%nspd,25)
                      CALL cylpts(xp,npd,cell%z1)
                      CALL checkdop(&
                           &          xp,npd,0,0,vacuum%nvac,1,ispin,.TRUE.,DIMENSION,atoms,&
                           &          sphhar,stars,sym,&
                           &          vacuum,cell,oneD,&
                           &          qpw,rho,rhtxy,rht)
                      !+odim
                   END IF
                   !         ---> m.t. boundaries
                   nat = 1
                   DO n = 1,atoms%ntype
                      CALL sphpts(xp,DIMENSION%nspd,atoms%rmt(n),atoms%pos(1,nat))
                      CALL checkdop(&
                           &                   xp,DIMENSION%nspd,n,nat,0,-1,ispin,.TRUE.,&
                           &                   dimension,atoms,sphhar,stars,sym,&
                           &                   vacuum,cell,oneD,&
                           &                   qpw,rho,rhtxy,rht)
                      nat = nat + atoms%neq(n)
                   ENDDO
                ENDDO
             ENDIF
             l_enpara = .FALSE.
             INQUIRE (file='enpara',exist=l_enpara)
             l_enpara = l_enpara.OR.input%l_inpXML
             !
             ! set up parameters for enpara-file
             !
             IF ((juDFT_was_argument("-genEnpara")).AND..NOT.l_enpara) THEN
                OPEN (40,file='enpara',form='formatted',status='unknown')
                enpara%lchange = .TRUE.
                enpara%llochg = .TRUE.
                
                DO ispin = 1, input%jspins
                   !
                   ! vacpar is taken as highest occupied level from atomic eigenvalues
                   ! vacpar (+0.3)  serves as fermi level for film or bulk
                   !
                   vacpar(1) = -999.9
                   DO n = 1,atoms%ntype
                      vacpar(1) = MAX( vacpar(1),eig(nst(n),ispin,n) )
                   ENDDO
                   IF (.NOT.input%film) vacpar(1) = vacpar(1) + 0.4
                   vacpar(2) = vacpar(1)


                   DO n = 1,atoms%ntype
                      enpara%skiplo(n,ispin) = 0
                      DO i = 0, 3
                         l_found(i) = .FALSE.
                         enpara%el0(i,n,ispin) = vacpar(1)
                      ENDDO
                      DO i = 1, atoms%nlo(n)
                         llo_found(i) = .FALSE.
                         enpara%ello0(i,n,ispin) = vacpar(1) - 1.5
                      ENDDO
                      !
                      ! take energy-parameters from atomic calculation
                      !
                      DO i = nst(n), 1,-1 
                         IF (.NOT.input%film) eig(i,ispin,n) = eig(i,ispin,n) + 0.4

                         IF (.NOT.l_found(lnum(i,n)).AND.(lnum(i,n).LE.3)) THEN
                            enpara%el0(lnum(i,n),n,ispin) = eig(i,ispin,n)
                            IF (enpara%el0(lnum(i,n),n,ispin).LT.input%ellow) THEN
                               enpara%el0(lnum(i,n),n,ispin) = vacpar(1)
                               l_found(lnum(i,n))  = .TRUE.
                            ENDIF
                            IF (enpara%el0(lnum(i,n),n,ispin).LT.input%elup) THEN
                               l_found(lnum(i,n))  = .TRUE.
                            ENDIF
                         ELSE
                            IF (l_found(lnum(i,n)).AND.(atoms%nlo(n).GT.0)) THEN
                               DO ilo = 1,atoms%nlo(n) 
                                  IF (atoms%llo(ilo,n).EQ.lnum(i,n)) THEN
                                     IF ( .NOT.llo_found(ilo) ) THEN
                                        enpara%ello0(ilo,n,ispin) = eig(i,ispin,n)
                                        IF (( enpara%ello0(ilo,n,ispin).GT.input%elup).OR.&
                                             &                      ( enpara%ello0(ilo,n,ispin).LT.input%ellow)) THEN
                                           enpara%ello0(ilo,n,ispin)= vacpar(1)
                                        ELSEIF (atoms%l_dulo(ilo,n)) THEN
                                           enpara%ello0(ilo,n,ispin)= enpara%el0(atoms%llo(ilo,n),n,ispin)
                                        ELSE
                                           enpara%skiplo(n,ispin) = enpara%skiplo(n,ispin) + 2*atoms%llo(ilo,n)+1
                                        ENDIF
                                        llo_found(ilo) = .TRUE.
                                        IF ((enpara%el0(atoms%llo(ilo,n),n,ispin)-enpara%ello0(ilo,n,ispin).LT.0.5)&
                                                                 .AND.(atoms%llo(ilo,n).GE.0)) THEN
                                           enpara%el0(atoms%llo(ilo,n),n,ispin) = vacpar(1)
                                           IF (atoms%l_dulo(ilo,n)) &
                                                &                              enpara%ello0(ilo,n,ispin)= enpara%el0(atoms%llo(ilo,n),n,ispin)
                                        ENDIF
                                     ENDIF
                                  ENDIF
                               ENDDO
                            ENDIF
                         ENDIF

                      ENDDO
                      IF (obsolete%lepr.EQ.1) THEN
                         DO i = 0, 3
                            enpara%el0(i,n,ispin) = enpara%el0(i,n,ispin) - vbar(ispin,n)
                         ENDDO
                         DO ilo = 1,atoms%nlo(n)
                            enpara%ello0(ilo,n,ispin) = enpara%ello0(ilo,n,ispin) - vbar(ispin,n)
                         ENDDO
                      ENDIF
                   ENDDO  ! atom types

                   IF (input%film) THEN
                      ! get guess for vacuum parameters
                      ! YM : in 1D case should be modified, but also not that bad like this
                      !
                      !           generate coulomb potential by integrating inward to z1
                      !
                      DO ivac = 1, vacuum%nvac
                         icorr_dummy = MIN(MAX(xcpot%icorr,0),5)  ! use LDA for this purpose
                         DO i=1,vacuum%nmz
                            sigm(i) = (i-1)*vacuum%delz*rht(i,ivac,ispin)
                         ENDDO
                         CALL qsf(vacuum%delz,sigm,vacpar(ivac),vacuum%nmz,0)
                         denz1(1) = rht(1,ivac,ispin)          ! get estimate for potential at
                         CALL  vxcall(6,icorr_dummy,input%krla,1,&    !               vacuum boundary&
                         &                     1,1,denz1,&
                              &                     vacxpot,vacpot)
                         ! seems to be the best choice for 1D not to substract vacpar
                         IF (.NOT.oneD%odi%d1) THEN
                            vacpot(1) = vacpot(1) - fpi_const*vacpar(ivac)
                         END IF
                         IF (obsolete%lepr.EQ.1) THEN
                            vacpar(ivac) = -0.2 - vacpot(1)
                            WRITE (6,'(" vacuum",i2," reference energy =",f12.6)')&
                                 &                 ivac,vacpot
                         ELSE
                            vacpar(ivac) = vacpot(1)
                         ENDIF
                      ENDDO
                      IF (vacuum%nvac.EQ.1) vacpar(2) = vacpar(1)

                   ENDIF
                   IF (obsolete%lepr.EQ.1) THEN
                      enpara%enmix = 0.3
                   ELSE
                      enpara%enmix = 1.0
                   ENDIF
                   !
                   ! write enpara-file
                   !
                   enpara%evac0(:,ispin)=vacpar(:SIZE(enpara%evac0,1))
                   CALL w_enpara(&
                        &                  atoms,ispin,input%film,enpara,&
                        &                  16)

                ENDDO  ! ispin


                CLOSE (40)
             ENDIF
             DEALLOCATE ( qpw,rhtxy,rht,xp,rat,eig,rh,rh1 )
             DEALLOCATE ( rhoss,vacpar,vbar,sigm )
             DEALLOCATE ( enpara%ello0,enpara%el0,enpara%lchange)
             DEALLOCATE ( enpara%skiplo,enpara%llochg,enpara%enmix,enpara%evac0)
          ENDIF ! mpi%irank == 0
          DEALLOCATE ( rho )
          !
        END SUBROUTINE stden
      END MODULE m_stden
