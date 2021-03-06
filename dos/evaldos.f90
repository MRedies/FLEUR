      MODULE m_evaldos
      CONTAINS
      SUBROUTINE evaldos(eig_id,input,banddos,vacuum,kpts,atoms,sym,noco,oneD,cell,results,dos,&
                         efermiarg,bandgap,l_mcd,mcd,slab,orbcomp,jDOS)
!----------------------------------------------------------------------
!
!     vk: k-vectors
!     wk weight of k-point (not used)
!     nevk no of eigenvalues
!     ev eigenvalue
!     qal partial charges
!             partial charge interstitial: qal(lmax*ntype+1...
!             partial charge vacuum      : qal(lmax*ntype+2...
!     qlay,qstar read in vacuum charges
!     qval partial charges in the vacuum
!             qval(m*nstars,neigd,nkptd):charge in m'th layer
!             qval(m*nstars+nstars,... ):star-resolved charge of that layer
!             qval(layers*nstars+....  ):same for second vacuum
!     ntb=max(nevk)
!
!----------------------------------------------------------------------
      USE m_types
      USE m_juDFT
      USE m_constants
      USE m_tetrados
      USE m_dostetra !confusing names (TODO:change/cleanup)
      USE m_dosbin
      USE m_ptdos
      USE m_smooth
      USE m_cdn_io
      IMPLICIT NONE
      INTEGER,INTENT(IN)             :: eig_id

      TYPE(t_oneD),INTENT(IN)        :: oneD
      TYPE(t_banddos),INTENT(IN)     :: banddos
      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_vacuum),INTENT(IN)      :: vacuum
      TYPE(t_noco),INTENT(IN)        :: noco
      TYPE(t_sym),INTENT(IN)         :: sym
      TYPE(t_cell),INTENT(IN)        :: cell
      TYPE(t_results),INTENT(IN)     :: results
      TYPE(t_dos),INTENT(IN)         :: dos
      TYPE(t_mcd),INTENT(IN)         :: mcd
      TYPE(t_slab),INTENT(IN)        :: slab
      TYPE(t_orbcomp),INTENT(IN)     :: orbcomp
      TYPE(t_jDOS),INTENT(IN)        :: jDOS
      TYPE(t_kpts),INTENT(IN)        :: kpts
      TYPE(t_atoms),INTENT(IN)       :: atoms

      REAL,    INTENT(IN) :: efermiarg, bandgap
      LOGICAL, INTENT(IN) :: l_mcd

!    locals
      INTEGER, PARAMETER ::  lmax= 4, ned = 1301
      INTEGER  i,s,v,index,jspin,k,l,l1,l2,ln,n,nl,ntb,col,jj,iBand
      INTEGER  icore,qdim,n_orb,ncored,jsp,n_jDOS
      REAL     de,efermi,emax,emin,qmt,sigma,totdos,efermiPrev
      REAL     e_up,e_lo,e_test1,e_test2,fac,sumwei,dk,eFermiCorrection
      LOGICAL  l_orbcomp,l_error,l_jDOS

      REAL     kx(kpts%nkpt),vkr(3,kpts%nkpt),ldos(0:3)
      REAL     ev(input%neig,kpts%nkpt),e(ned),gpart(ned,atoms%ntype)
      REAL     e_grid(ned+1),spect(ned,3*atoms%ntype),ferwe(input%neig,kpts%nkpt)
      REAL,    ALLOCATABLE :: qal(:,:,:),qval(:,:,:),qlay(:,:,:),g(:,:)
      REAL,    ALLOCATABLE :: mcd_local(:,:,:)
      REAL,    ALLOCATABLE :: qvac(:,:)
      CHARACTER(len=2) :: spin12(2),ch_mcd(3)
      CHARACTER(len=8) :: chntype*2,chform*19
      DATA spin12/'.1' , '.2'/
      DATA ch_mcd/'.+' , '.-' , '.0'/

      ncored =  MAX(0,MAXVAL(mcd%ncore))
      qdim = lmax*atoms%ntype+3
      l_orbcomp = banddos%l_orb
      l_jDOS = banddos%l_jDOS
      n_orb = 0
      n_jDOS = 0
      IF (banddos%ndir.EQ.-3) THEN
        qdim = 2*slab%nsld
        n_orb = 0
        IF (banddos%l_orb) THEN
           n_orb = banddos%orbCompAtom
           WRITE (*,*) 'DOS: orbcomp',n_orb
           qdim = 23
        END IF
      ENDIF

      IF(l_jDOS) THEN
         qdim = 7
         n_jDOS = banddos%jDOSAtom
         WRITE(*,*) 'DOS: jDOS', n_jDOS
      ENDIF

      ALLOCATE( qal(qdim,input%neig,kpts%nkpt),&
     &          qval(vacuum%nstars*vacuum%layers*vacuum%nvac,input%neig,kpts%nkpt),&
     &          qlay(input%neig,vacuum%layerd,2))
      IF (l_mcd) THEN
         ALLOCATE(mcd_local(3*atoms%ntype*ncored,input%neig,kpts%nkpt) )
      ELSE
         ALLOCATE(mcd_local(0,0,0))
      ENDIF
!
! scale energies
      sigma = banddos%sig_dos*hartree_to_ev_const
      emin =min(banddos%e1_dos*hartree_to_ev_const,banddos%e2_dos*hartree_to_ev_const)
      emax =max(banddos%e1_dos*hartree_to_ev_const,banddos%e2_dos*hartree_to_ev_const)
      efermi = efermiarg*hartree_to_ev_const

      WRITE (oUnit,'(a)') 'DOS-Output is generated!'

      IF ( NINT((emax - emin)/sigma) > ned ) THEN
        WRITE(oUnit,*) 'sig_dos too small for DOS smoothing:'
        WRITE(oUnit,*) 'Reduce energy window or enlarge banddos%sig_dos!'
        WRITE(oUnit,*) 'For now: setting sigma to zero !'
        sigma = 0.0
      ENDIF

      WRITE (oUnit,*) 'sigma=   ' , sigma
      WRITE (oUnit,*) 'emax=   ' , emax
      WRITE (oUnit,*) 'emin=   ' , emin
      WRITE (oUnit,*) 'ef_inp=   ' , efermi
!
!     create energy grid
      emax = emax - efermi
      emin = emin - efermi
      de = (emax-emin)/(ned-1)
      DO i=1,ned
         e(i) = emin + (i-1)*de
      ENDDO

      IF ( l_mcd ) THEN ! create an energy grid for mcd-spectra
        e_lo =  9.9*10.0**9
        e_up = -9.9*10.0**9
        DO jspin = 1,input%jspins
          DO n = 1,atoms%ntype
            DO icore = 1 , mcd%ncore(n)
              e_lo = min(mcd%e_mcd(n,jspin,icore),e_lo)
              e_up = max(mcd%e_mcd(n,jspin,icore),e_up)
            ENDDO
          ENDDO
        ENDDO
        e_lo = e_lo*hartree_to_ev_const - efermi - emax
        e_up = e_up*hartree_to_ev_const - efermi
        de = (e_up-e_lo)/(ned-1)
        DO i=1,ned
          e_grid(i) = e_lo + (i-1)*de
          spect(i,:) = 0.0
        ENDDO
        e_grid(ned+1) = e_lo + ned*de
      ENDIF

      DO jspin = 1,MERGE(1,input%jspins,l_jDOS)
         ntb = 0
         DO k = 1,kpts%nkpt

            qal(:,:,k) = 0.0
            qval(:,:,k) = 0.0

            jsp = MERGE(1,jspin,noco%l_noco)
            ntb = max(ntb,results%neig(k,jsp))
            IF (l_mcd) mcd_local(:,:,k) = RESHAPE(mcd%mcd(:,1:ncored,:,k,jspin),(/3*atoms%ntype*ncored,input%neig/))
            IF (.NOT.l_orbcomp.AND..NOT.l_jDOS) THEN
               qal(1:lmax*atoms%ntype,:,k)=reshape(dos%qal(0:,:,:,k,jspin),(/lmax*atoms%ntype,size(dos%qal,3)/))
               qal(lmax*atoms%ntype+2,:,k)=dos%qvac(:,1,k,jspin) ! vacuum 1
               qal(lmax*atoms%ntype+3,:,k)=dos%qvac(:,2,k,jspin) ! vacuum 2
               qal(lmax*atoms%ntype+1,:,k)=dos%qis(:,k,jspin)    ! interstitial
            ELSE IF(l_orbcomp) THEN
               IF (n_orb == 0) THEN
                  qal(1:slab%nsld,:,k)             = slab%qintsl(:,:,k,jspin)
                  qal(slab%nsld+1:2*slab%nsld,:,k) = slab%qmtsl(:,:,k,jspin)
               ELSE
                  DO i = 1, 23
                     DO l = 1, results%neig(k,jsp)
                        qal(i,l,k) = orbcomp%comp(l,i,n_orb,k,jspin)*orbcomp%qmtp(l,n_orb,k,jspin)/10000.
                     END DO
                     DO l = results%neig(k,jsp)+1, input%neig
                        qal(i,l,k) = 0.0
                     END DO
                  END DO
               END IF
            ELSE IF(l_jDOS) THEN
               i = 0
               DO l= 0, 3
                  DO jj = 1, MERGE(1,2,l==0)
                     i = i+1
                     DO iBand = 1, results%neig(k,jsp)
                        qal(i,iBand,k) = jDOS%comp(iBand,l,jj,n_jDOS,k)*jDOS%qmtp(iBand,n_jDOS,k)/10000.
                     END DO
                     DO iBand = results%neig(k,jsp)+1, input%neig
                        qal(i,iBand,k) = 0.0
                     END DO
                  ENDDO
               ENDDO
            END IF
!
!     set vacuum partial charge zero, if bulk calculation
!     otherwise, write vacuum charge in correct arrays
!
            IF ((.NOT.input%film).AND.(banddos%ndir.NE.-3).AND..NOT.l_jDOS) THEN
               DO n = 1,input%neig
                  qal(lmax*atoms%ntype+2,n,k) = 0.0
                  qal(lmax*atoms%ntype+3,n,k) = 0.0
               ENDDO
            ELSEIF ( banddos%vacdos .and. input%film ) THEN
               DO i = 1,results%neig(k,jsp)
                  DO v = 1,vacuum%nvac
                     DO l = 1,vacuum%layers
                        index = (l-1)*vacuum%nstars + (v-1)*(vacuum%nstars*vacuum%layers) + 1
                        qval(index,i,k) = qlay(i,l,v)
                        DO s = 1,vacuum%nstars - 1
                           qval(index+s,i,k) = real(dos%qstars(s,i,l,v,k,jspin))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
!
!     calculate interstitial dos if not noco
!     in the noco case, qis has been calculated in pwden and is read in from tmp_dos
!
            IF ((.NOT.noco%l_noco).AND.(banddos%ndir.NE.-3).AND..NOT.l_jDOS) THEN
               DO i = 1 , input%neig
                  qal(lmax*atoms%ntype+1,i,k) = 1.
                  DO nl = 1 , atoms%ntype
                     l1 = lmax*(nl-1) + 1
                     l2 = lmax*nl
                     qmt=0.0
                     DO l = l1 , l2
                        qmt = qmt + qal(l,i,k)*atoms%neq(nl)
                     ENDDO
                     qal(lmax*atoms%ntype+1,i,k) = qal(lmax*atoms%ntype+1,i,k) - qmt
                  ENDDO
                  qal(lmax*atoms%ntype+1,i,k) = qal(lmax*atoms%ntype+1,i,k)&
                       -qal(lmax*atoms%ntype+2,i,k)*(3-vacuum%nvac) -qal(lmax*atoms%ntype+3,i,k)*(vacuum%nvac-1)
               ENDDO
            ENDIF
!
!---- >     convert eigenvalues to ev and shift them by efermi
!
            DO i = 1 , results%neig(k,jsp)
               ev(i,k) = results%eig(i,k,jsp)*hartree_to_ev_const - efermi
            ENDDO
            DO i = results%neig(k,jsp) + 1, input%neig
               ev(i,k) = 9.9e+99
            ENDDO
!
!
         ENDDO                                                 ! end of k-point loop

         IF ( .not.l_mcd ) THEN
            ALLOCATE (g(ned,qdim))
         ELSE
            ALLOCATE (g(ned,3*atoms%ntype*ncored))
         ENDIF

         IF ( (input%bz_integration.EQ.2 .OR. input%bz_integration.EQ.3) .and. &
             .not.l_mcd .and. banddos%ndir.NE.-3 .and. .not.l_jDOS) THEN
            !
            ! DOS calculation: use triangular method!!
            !
            IF(input%bz_integration.EQ.2) THEN
               IF ( input%film ) THEN
                  CALL ptdos(input%jspins,ned,qdim,ntb,kpts,ev,qal,e, g)
               ELSE
                  write(*,*) efermi
                  CALL tetra_dos(qdim,input%neig,ned,kpts,efermi,e,results%neig(:,jsp),ev,qal,g)
                  IF (input%jspins.EQ.1) g = 2.0 * g
               ENDIF
            ELSEIF(input%bz_integration.EQ.3) THEN
               !Alternative tetrahedron method
               CALL dostetra(kpts,input,qdim,ned,e,results%neig(:,jsp),ev,qal,g)
            ENDIF
         ELSE
            !
            ! DOS calculation: use histogram method
            !
            IF ( .not.l_mcd ) THEN
               CALL dos_bin(input%jspins,qdim,ned,emin,emax,input%neig,kpts%nkpt,&
                            results%neig(:,jsp),kpts%wtkpt(1:kpts%nkpt),ev,qal, g)
            ELSE
               CALL dos_bin(input%jspins,3*atoms%ntype*ncored,ned,emin,emax,ntb,kpts%nkpt,&
                            results%neig(:,jsp),kpts%wtkpt(1:kpts%nkpt),ev(1:ntb,1:kpts%nkpt),&
                            mcd_local(1:3*atoms%ntype*ncored,1:ntb,1:kpts%nkpt), g)
            ENDIF
         ENDIF
!
!---- >     smoothening
!
         IF ( .not.l_mcd ) THEN
            IF ( sigma.GT.0.0 ) THEN
              DO ln = 1 , qdim
                CALL smooth(e,g(1,ln),sigma,ned)
              ENDDO
            ENDIF

!*** sum up for all atoms

         IF (banddos%ndir.NE.-3.AND..NOT.l_jDOS) THEN
            DO l = 1 , atoms%ntype
               l1 = lmax*(l-1) + 1
               l2 = lmax*l
               DO i = 1 , ned
                  gpart(i,l) = 0.0
                  DO nl = l1 , l2
                     gpart(i,l) = gpart(i,l) + g(i,nl)
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF (n_orb == 0.AND..not.l_jDOS) THEN
            DO l = 1, slab%nsld
               nl = slab%nsld+l
               DO i = 1 , ned
                  gpart(i,l) = g(i,l) + g(i,nl)
               ENDDO
            ENDDO
         ENDIF

!**** write out DOS
         OPEN (18,FILE='DOS'//spin12(jspin))
         IF (atoms%ntype >= 20.AND.banddos%projdos.NE.0) OPEN (1337,FILE="PROJDOS"//spin12(jspin))

         DO i = 1 , ned
            totdos = 0.0
            IF (banddos%ndir.NE.-3.AND..NOT.l_jDOS) THEN
               DO nl = 1 , atoms%ntype
                  totdos = totdos + gpart(i,nl)*atoms%neq(nl)
               ENDDO
               totdos = totdos + g(i,lmax*atoms%ntype+1) + g(i,lmax*atoms%ntype+2) *&
                    (3 - vacuum%nvac) + g(i,lmax*atoms%ntype+3)*(vacuum%nvac - 1)
               IF (atoms%ntype < 20) THEN
                  WRITE (18,99001)  e(i),totdos,g(i,lmax*atoms%ntype+1), &
                       g(i,lmax*atoms%ntype+2),g(i,lmax*atoms%ntype+3),&
                       (gpart(i,l),l=1,atoms%ntype), (g(i,l),l=1,atoms%ntype*lmax)
               ELSE
                  WRITE (18,99001)  e(i),totdos,g(i,lmax*atoms%ntype+1), &
                       g(i,lmax*atoms%ntype+2),g(i,lmax*atoms%ntype+3), (gpart(i,l),l=1,atoms%ntype)
                  IF(banddos%projdos.NE.0) WRITE (1337,99001)  e(i),(g(i,l),l=lmax*(banddos%projdos-1)+1,lmax*banddos%projdos)
               ENDIF
             ELSEIF (l_orbcomp) THEN
                IF(n_orb == 0) THEN
                   DO nl = 1, slab%nsld
                      totdos = totdos + gpart(i,nl)
                   ENDDO
                   WRITE (18,99001)  e(i),totdos,(gpart(i,nl),nl=1,slab%nsld), (g(i,l),l=1,2*slab%nsld)
                ELSE
                   DO nl = 1 , 23
                      totdos = totdos + g(i,nl)
                   ENDDO
                   WRITE (18,99001)  e(i),totdos,(g(i,l),l=1,23)
                ENDIF
             ELSEIF(l_jDOS) THEN
                DO nl = 1 , 7
                   totdos = totdos + g(i,nl)
                ENDDO
                ldos= 0.0
                col = 0
                DO l = 0, 3
                  DO jj = 1, MERGE(1,2,l==0)
                    col = col +1
                    ldos(l) = ldos(l) + g(i,col)
                  ENDDO
                ENDDO
                WRITE(18,99001) e(i),totdos,(ldos(l),l=0,3),(g(i,l),l=2,7)
             ENDIF
         ENDDO
         CLOSE (18)
         IF (atoms%ntype >= 20.AND.banddos%projdos.NE.0) CLOSE (1337)

         ELSE
           write(*,'(4f15.8)') ((mcd%e_mcd(n,jspin,i),n=1,atoms%ntype),i=1,ncored)
           write(*,*)
           write(*,'(4f15.8)') (g(800,n),n=1,3*atoms%ntype*ncored)
           write(*,*)
           write(*,'(4f15.8)') (mcd_local(n,10,8),n=1,3*atoms%ntype*ncored)
           DO n = 1,atoms%ntype
             DO l = 1 , ned
               DO icore = 1 , mcd%ncore(n)
                 DO i = 1 , ned-1
                   IF (e(i).GT.0) THEN     ! take unoccupied part only
                   e_test1 = -e(i) - efermi +mcd%e_mcd(n,jspin,icore)*hartree_to_ev_const
                   e_test2 = -e(i+1)-efermi +mcd%e_mcd(n,jspin,icore)*hartree_to_ev_const
                   IF ((e_test2.LE.e_grid(l)).AND. (e_test1.GT.e_grid(l))) THEN
                     fac = (e_grid(l)-e_test1)/(e_test2-e_test1)
                     DO k = 3*(n-1)+1,3*(n-1)+3
                       spect(l,k) = spect(l,k)+ g(i,3*atoms%ntype*(icore-1)+k)&
                           *(1.-fac) + fac * g(i+1,3*atoms%ntype*(icore-1)+k)
                     ENDDO
                   ENDIF
                   ENDIF
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
           CLOSE (18)
         ENDIF
         DEALLOCATE (g)
!
!------------------------------------------------------------------------------
!     now calculate the VACOS
!------------------------------------------------------------------------------

         IF ( banddos%vacdos .and. input%film ) THEN
            ALLOCATE(g(ned,vacuum%nstars*vacuum%layers*vacuum%nvac))
            IF(kpts%ntet.EQ.0) THEN
               CALL juDFT_error("VACDOS requires a kpoint set with generated triangles",calledby="evaldos")
            ENDIF
!            CALL ptdos(
!     >                 emin,emax,jspins,ned,nstars*nvac*layers,neigd,
!     >                 ntria,as,atr,2*nkpt,itria,nkptd,ev,qval,e,
!     <                 g)
            CALL ptdos(input%jspins,ned,vacuum%nstars*vacuum%nvac*vacuum%layers,&
                       ntb,kpts,ev,qval,e,g)

!---- >     smoothening
            IF ( sigma.GT.0.0 ) THEN
               DO ln = 1 , vacuum%nstars*vacuum%nvac*vacuum%layers
                  CALL smooth(e,g(1,ln),sigma,ned)
               ENDDO
            ENDIF

!     write VACDOS

            OPEN (18,FILE='VACDOS'//spin12(jspin))
!            WRITE (18,'(i2,25(2x,i3))') Layers , (Zlay(l),l=1,Layers)
            DO i = 1 , ned
             WRITE (18,99001) e(i) , (g(i,l),l=1,vacuum%Layers*vacuum%Nstars*vacuum%Nvac)
            ENDDO
            CLOSE (18)
            DEALLOCATE(g)
         ENDIF
!
!------------------------------------------------------------------------------
!     for bandstructures
!------------------------------------------------------------------------------

         IF (banddos%ndir == -4) THEN
            eFermiCorrection = 0.0
            IF(bandgap.LT.(8.0*input%tkb*hartree_to_ev_const)) THEN
               CALL readPrevEFermi(eFermiPrev,l_error)
               IF(.NOT.l_error) THEN
                  WRITE(*,*) 'Fermi energy is automatically corrected in bands.* files.'
                  WRITE(*,*) 'It is consistent with last calculated density!'
                  WRITE(*,*) 'No manual correction (e.g. in band.gnu file) required.'
                  eFermiCorrection = (eFermiPrev-efermiarg)*hartree_to_ev_const
               ELSE
                  WRITE(*,*) 'Fermi energy in bands.* files may not be consistent with last density.'
                  WRITE(*,*) 'Please correct it manually (e.g. in band.gnu file).'
               END IF
            END IF

            OPEN (18,FILE='bands'//spin12(jspin))
            ntb = minval(results%neig(:,jsp))
            kx(1) = 0.0
            vkr(:,1)=matmul(kpts%bk(:,1),cell%bmat)
            DO k = 2, kpts%nkpt

               vkr(:,k)=matmul(kpts%bk(:,k),cell%bmat)
               dk = (vkr(1,k)-vkr(1,k-1))**2 + (vkr(2,k)-vkr(2,k-1) )**2 + &
                    (vkr(3,k)-vkr(3,k-1))**2
               kx(k) = kx(k-1) + sqrt(dk)
            ENDDO
            DO i = 1, ntb
               DO k = 1, kpts%nkpt
                  write(18,'(2f15.9)') kx(k),ev(i,k)-eFermiCorrection
               ENDDO
            ENDDO
            CLOSE (18)
         ENDIF

      ENDDO
!
!------------------------------------------------------------------------------
!     for MCD calculations ...
!------------------------------------------------------------------------------

      IF (l_mcd) THEN
        WRITE (chntype,'(i2)') atoms%ntype+1
        chform = '('//chntype//'f15.8)'
        IF ( sigma.GT.0.0 ) THEN
           IF ( l_mcd ) THEN
             DO ln = 1 , 3*atoms%ntype
               CALL smooth(e_grid,spect(1,ln),sigma,ned)
             ENDDO
           ENDIF
        ENDIF
        DO l = 1,3
          OPEN (18,FILE='MCD_SPEC'//ch_mcd(l))
          DO i = 1 , ned
          WRITE (18,FMT=chform) e_grid(i),(spect(i,3*(n-1)+l),n=1,atoms%ntype)
          ENDDO
          CLOSE (18)
        ENDDO
      ENDIF

      DEALLOCATE(qal,qval,qlay)
      IF (l_mcd) DEALLOCATE( mcd_local )
99001 FORMAT (f10.5,110(1x,e10.3))

      END SUBROUTINE evaldos
      END MODULE m_evaldos
