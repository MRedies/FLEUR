      MODULE m_inped
      USE m_juDFT
!     *******************************************************
!     read in input parameters
!     modified to include also empty spheres (z=1.e-10)
!     r.p. aug. 1990
!     now includes readin of k-points          * shz Apr.96 *
!     modified to include all exchange correlation potentials
!     and relativistic correction to vxc
!     r.pentcheva dec. 1995
!
!ta+
!     igrd=0: no gradient correction.
!     igrd=1: pw-91. icorr=6.

!     ndvgrd: number of points used to calculate the numerical
!c           derivatives. 6 is biggest and presumably the best.
!     ntimec: if ntime ge ntimec, iteration stops with code=2.
!     distc : distance of convergence in charge criterion.
!     tendfc: read-in in mhtree.
!c            if tendf (total energy diff. in mili-hartree from former
!c            tenrg) becomes less than tendfc, ntime=ntime+1.
!     chng  : charge-negative.
!c             if(ro.lt.chng) ineg=1 and stop.
!ta-
!     *******************************************************
!
      CONTAINS
        SUBROUTINE inped(atoms,vacuum,input,banddos,xcpot,sym,&
                         cell,sliceplot,noco,&
                         stars,oneD,hybinp,kpts,a1,a2,a3,namex,relcor,latnam,namgrp,grid)
          USE m_rwinp
          !USE m_chkmt
          USE m_inpnoco
          USE m_constants
          USE m_types_atoms
          USE m_types_vacuum
          USE m_types_input
          USE m_types_banddos
          USE m_types_xcpot_inbuild_nofunction
          USE m_types_sym
          USE m_types_cell
          USE m_types_sliceplot
          USE m_types_noco
          USE m_types_stars
          USE m_types_oneD
          USE m_types_hybinp
          USE m_types_kpts
          USE m_constants
          USE m_setlomap
          IMPLICIT NONE
          !     ..
          !     .. Scalar Arguments ..
          TYPE(t_atoms),     INTENT(INOUT) :: atoms
          TYPE(t_vacuum),    INTENT(INOUT) :: vacuum
          TYPE(t_input),     INTENT(INOUT) :: input
          TYPE(t_banddos),   INTENT(INOUT) :: banddos
          TYPE(t_xcpot_inbuild_nf),     INTENT(INOUT) :: xcpot
          TYPE(t_sym),       INTENT(INOUT) :: sym
          TYPE(t_cell),      INTENT(INOUT) :: cell
          TYPE(t_sliceplot), INTENT(INOUT) :: sliceplot
          TYPE(t_noco),      INTENT(INOUT) :: noco
          TYPE(t_stars),     INTENT(INOUT) :: stars
          TYPE(t_oneD),      INTENT(INOUT) :: oneD
          TYPE(t_hybinp),    INTENT(INOUT) :: hybinp
          TYPE(t_kpts),      INTENT(INOUT) :: kpts
          REAL,              INTENT(OUT)   :: a1(3)
          REAL,              INTENT(OUT)   :: a2(3)
          REAL,              INTENT(OUT)   :: a3(3)
          CHARACTER(len=4),  INTENT(OUT)   :: namex
          CHARACTER(len=12), INTENT(OUT)   :: relcor
          CHARACTER(len=4),INTENT(OUT)     ::latnam,namgrp
          INTEGER,INTENT(OUT):: grid(3)
          !     .. Local Scalars ..
          REAL dr,dtild,r,kmax1,dvac1,zp
          INTEGER i,iz,j,n,n1,na,ntst,nn,ios
          LOGICAL l_gga,l_test,l_vca
          CHARACTER(len=2)  :: str_up,str_do

          !     ..
          !     .. Local Arrays ..
          CHARACTER(3) noel(atoms%ntype)
          CHARACTER(8) llr(0:1)
          INTEGER  jri1(atoms%ntype),lmax1(atoms%ntype)
          REAL    rmt1(atoms%ntype),dx1(atoms%ntype)
          real     ::scalecell
          !     ..
          !     .. Data statements ..
          DATA llr(0)/'absolute'/,llr(1)/'floating'/

          !

          a1(:) = 0
          a2(:) = 0
          a3(:) = 0

          na = 0

          CALL rw_inp('r',atoms,vacuum,input,stars,sliceplot,banddos,&
               cell,sym,xcpot,noco,oneD,hybinp,kpts, noel,namex,relcor,a1,a2,a3,latnam,grid,namgrp,scalecell)

          input%l_core_confpot=.TRUE. !this is the former CPP_CORE switch!
          input%l_useapw=.FALSE.      !this is the former CPP_APW switch!
          atoms%lapw_l(:) = -1
          IF(ANY(atoms%lnonsph(:).GT.39)) THEN
             WRITE(*,*) 'lnonsph > 39 detected. This is interpreted as an intention to use APW+lo.'
             WRITE(*,*) 'Switching to APW+lo!'
             input%l_useapw=.TRUE.
          END IF
          !---> pk non-collinear
          !---> read the angle information from nocoinf
          noco%qss_inp(:) = 0.0
          IF (noco%l_noco) THEN
             CALL inpnoco(atoms,input,sym,vacuum,noco)
          ELSE
             noco%l_ss = .FALSE.
             noco%l_mperp = .FALSE.
             noco%l_constr = .FALSE.
             noco%mix_b = 0.0
             noco%l_relax(:) = .FALSE.
             noco%alph_inp(:) = 0.0
             noco%beta_inp(:) = 0.0
          ENDIF
          !---> pk non-collinear

8010      FORMAT (/,/,4x,10a8,/,/)
          !--->    the menu for namgrp can be found in subroutine spgset
          WRITE (oUnit,FMT=8030) latnam,namgrp,sym%invs,sym%zrfs,sym%invs2,input%jspins
8030      FORMAT (' lattice=',a3,/,' name of space group=',a4,/,' inversion symmetry=   ',l1&
               ,/,' z-reflection symmetry=',l1,/,' vacuum-inversion symm=',l1,/,' jspins=',i1)

          IF (input%film.AND.(sym%invs.OR.sym%zrfs)) THEN
             IF ( (sym%invs.AND.sym%zrfs).NEQV.sym%invs2 ) THEN
                WRITE (oUnit,*) 'Settings of inversion and z-reflection symmetry='
                WRITE (oUnit,*) 'are inconsistent with vacuum-inversion symmetry!'
                CALL juDFT_error("invs, zrfs and invs2 do not match!",calledby ="inped")
             ENDIF
          ENDIF


          IF (ALL(a1.EQ.0.)) THEN
             WRITE (oUnit,'(a4,3f10.5,a8,a4)') 'a1 =',a1(:),' latnam=',latnam
             CALL juDFT_error("latnam",calledby ="inped")
          ENDIF
          dtild=a3(3)
          IF (scaleCell.EQ.0.0) scaleCell = 1.0
          vacuum%dvac = scaleCell*vacuum%dvac
          dtild = scaleCell*dtild
          !+odim
          IF (.NOT.oneD%odd%d1) THEN
             IF ((dtild-vacuum%dvac.LT.0.0).AND.input%film) THEN
                WRITE(oUnit,'(2(a7,f10.5))') 'dtild:',dtild,' dvac:',vacuum%dvac
                CALL juDFT_error("dtild < dvac",calledby="inped")
             ENDIF
          ELSE
             IF (vacuum%dvac.GE.SQRT(a1(1)**2 + a1(2)**2).OR. vacuum%dvac.GE.SQRT(a2(1)**2 + a2(2)**2)) THEN
                CALL juDFT_error("one-dim: dvac >= amat(1,1) or amat(2,2)" ,calledby ="inped")
             END IF
          ENDIF
          !-odim
          vacuum%nvac = 2
          IF (sym%zrfs .OR. sym%invs) vacuum%nvac = 1
          IF (oneD%odd%d1) vacuum%nvac = 1
          cell%z1 = vacuum%dvac/2
          vacuum%nmz = vacuum%nmzd
          vacuum%delz = 25.0/vacuum%nmz
          IF (oneD%odd%d1) vacuum%delz = 20.0/vacuum%nmz
          IF (vacuum%nmz>vacuum%nmzd)  CALL juDFT_error("nmzd",calledby ="inped")
          vacuum%nmzxy = vacuum%nmzxyd
          IF (vacuum%nmzxy>vacuum%nmzxyd)  CALL juDFT_error("nmzxyd",calledby ="inped")
          a1(:) = scaleCell*a1(:)
          a2(:) = scaleCell*a2(:)
          a3(:) = scaleCell*a3(:)
          WRITE (oUnit,FMT=8050) scaleCell
8050      FORMAT (' unit cell scaled by    ',f10.6)
          WRITE (oUnit,FMT=8060) cell%z1
8060      FORMAT (' the vacuum begins at z=',f10.6)
          WRITE (oUnit,FMT=8070) dtild/2.
8070      FORMAT (' dtilda/2=              ',f10.6)
          !     set up bravais matrices of real and reciprocal lattices
          cell%amat(:,1) = a1(:)
          cell%amat(:,2) = a2(:)
          cell%amat(:,3) = a3(:)
          call cell%init(dot_product(atoms%volmts(:),atoms%neq(:)))
          !CALL inv3(cell%amat,cell%bmat,cell%omtil)
          !cell%bmat(:,:) = tpi_const*cell%bmat(:,:)
          !cell%bbmat=MATMUL(cell%bmat,TRANSPOSE(cell%bmat))
          !cell%omtil = ABS(cell%omtil)

          IF (input%film .AND. .NOT.oneD%odd%d1) THEN
             cell%vol = cell%omtil/cell%amat(3,3)*vacuum%dvac
             cell%area = cell%omtil/cell%amat(3,3)
             !-odim
          ELSEIF (oneD%odd%d1) THEN
             cell%area = tpi_const*cell%amat(3,3)
             cell%vol = pi_const*(vacuum%dvac**2)*cell%amat(3,3)/4.
             !+odim
          ELSE
             cell%vol = cell%omtil
             cell%area = cell%amat(1,1)*cell%amat(2,2)-cell%amat(1,2)*cell%amat(2,1)
             IF (cell%area.LT.1.0e-7) THEN
                IF (latnam.EQ.'any') THEN
                   cell%area = 1.
                ELSE
                   CALL juDFT_error("area = 0",calledby ="inped")
                ENDIF
             ENDIF
          ENDIF


          WRITE (oUnit,FMT=8080)
8080      FORMAT (/,/,1x,'bravais matrices of real and reciprocal lattices', /)
          DO i = 1,3
             WRITE (oUnit,FMT=8090) (cell%amat(i,j),j=1,3), (cell%bmat(i,j),j=1,3)
          ENDDO
8090      FORMAT (3x,3f10.6,3x,3f10.6)
          WRITE (oUnit,FMT=8100) cell%omtil,cell%vol,cell%area
8100      FORMAT (/,4x,'the volume of the unit cell omega-tilda=',f12.6,/, 10x,'the volume of the unit cell omega=',&
               f12.6,/,2x, 'the area of the two-dimensional unit cell=',f12.6)
          WRITE (oUnit,FMT=8120) namex,relcor
8120      FORMAT (1x,'exchange-correlation: ',a4,2x,a12,1x,'correction')
          xcpot%l_inbuild=.true.
          xcpot%inbuild_name=namex
          xcpot%l_relativistic=relcor.EQ.'relativistic'
          CALL xcpot%init(atoms%ntype)
!!$          xcpot%icorr = -99
!!$
!!$          !     l91: lsd(igrd=0) with dsprs=1.d-19 in pw91.
!!$          IF (namex.EQ.'exx ') xcpot%icorr = icorr_exx
!!$          IF (namex.EQ.'hf  ') xcpot%icorr = icorr_hf
!!$          IF (namex.EQ.'l91 ') xcpot%icorr = -1
!!$          IF (namex.EQ.'x-a ') xcpot%icorr =  0
!!$          IF (namex.EQ.'wign') xcpot%icorr =  1
!!$          IF (namex.EQ.'mjw')  xcpot%icorr =  2
!!$          IF (namex.EQ.'hl')   xcpot%icorr =  3
!!$          IF (namex.EQ.'bh')   xcpot%icorr =  3
!!$          IF (namex.EQ.'vwn')  xcpot%icorr =  4
!!$          IF (namex.EQ.'pz')   xcpot%icorr =  5
!!$          IF (namex.EQ.'pw91') xcpot%icorr =  6
!!$          !     pbe: easy_pbe [Phys.Rev.Lett. 77, 3865 (1996)]
!!$          !     rpbe: rev_pbe [Phys.Rev.Lett. 80, 890 (1998)]
!!$          !     Rpbe: Rev_pbe [Phys.Rev.B 59, 7413 (1999)]
!!$          IF (namex.EQ.'pbe')  xcpot%icorr =  7
!!$          IF (namex.EQ.'rpbe') xcpot%icorr =  8
!!$          IF (namex.EQ.'Rpbe') xcpot%icorr =  9
!!$          IF (namex.EQ.'wc')   xcpot%icorr = 10
!!$          !     wc: Wu & Cohen, [Phys.Rev.B 73, 235116 (2006)]
!!$          IF (namex.EQ.'PBEs') xcpot%icorr = 11
!!$          !     PBEs: PBE for solids ( arXiv:0711.0156v2 )
!!$          IF (namex.EQ.'pbe0') xcpot%icorr = icorr_pbe0
!!$          !     hse: Heyd, Scuseria, Ernzerhof, JChemPhys 118, 8207 (2003)
!!$          IF (namex.EQ.'hse ') xcpot%icorr = icorr_hse
!!$          IF (namex.EQ.'vhse') xcpot%icorr = icorr_vhse
!!$          ! local part of HSE
!!$          IF (namex.EQ.'lhse') xcpot%icorr = icorr_hseloc
!!$
!!$          IF (xcpot%icorr == -99) THEN
!!$             WRITE(oUnit,*) 'Name of XC-potential not recognized. Use one of:'
!!$             WRITE(oUnit,*) 'x-a,wign,mjw,hl,bh,vwn,pz,l91,pw91,pbe,rpbe,Rpbe,wc,PBEs,pbe0,hf,hse,lhse'
!!$             CALL juDFT_error("Wrong name of XC-potential!",calledby="inped")
!!$          ENDIF
!!$          xcpot%krla = 0
!!$          IF (relcor.EQ.'relativistic') THEN
!!$             xcpot%krla = 1
!!$
!!$          ENDIF

!!$          IF (xcpot%icorr.EQ.0) WRITE(oUnit,*) 'WARNING: using X-alpha for XC!'
!!$          IF (xcpot%icorr.EQ.1) WRITE(oUnit,*) 'INFO   : using Wigner  for XC!'
!!$          IF ((xcpot%icorr.EQ.2).AND.(namex.NE.'mjw')) WRITE(oUnit,*) 'CAUTION: using MJW(BH) for XC!'
!!$
!!$          !+guta
!!$          IF ((xcpot%icorr.EQ.-1).OR.(xcpot%icorr.GE.6)) THEN

          !-guta
          !     specification of atoms

          cell%volint = cell%vol

          DO  n = 1,atoms%ntype
             IF (TRIM(ADJUSTL(noel(n))).NE.TRIM(ADJUSTL(namat_const(atoms%nz(n))))) THEN
                CALL trans(namat_const(n),str_up,str_do)
                IF ( (TRIM(ADJUSTL(noel(n))).NE.TRIM(ADJUSTL(str_up))) .OR.&
                     &        (TRIM(ADJUSTL(noel(n))).NE.TRIM(ADJUSTL(str_do))) ) THEN
                   WRITE(oUnit,*) 'Element ',noel(n),' does not match Z = ',atoms%nz(n)
                   CALL juDFT_warn ("Element name and nuclear number do not match!" ,calledby ="inped")
                ENDIF
             ENDIF
             WRITE (oUnit,8140) noel(n),atoms%nz(n),atoms%econf(n)%num_core_states,atoms%lmax(n),atoms%jri(n),atoms%rmt(n),atoms%dx(n)
8140         FORMAT (a3,i3,3i5,2f10.6)
             IF (atoms%jri(n)>atoms%jmtd)  CALL juDFT_error("jmtd",calledby ="inped")
             atoms%zatom(n) = atoms%nz(n)
             IF (atoms%nz(n).EQ.0) atoms%zatom(n) = 1.e-10
             !
             ! check for virtual crystal approximation
             !
             l_vca = .FALSE.
             INQUIRE (file="vca.in", exist=l_vca)
             IF (l_vca) THEN
                OPEN (17,file='vca.in',form='formatted')
                DO nn = 1, n
                   READ (17,*,IOSTAT=ios) ntst,zp
                   IF (ios /= 0) EXIT
                   IF (ntst == n) THEN
                      atoms%zatom(n) = atoms%zatom(n) + zp
                   ENDIF
                ENDDO
                CLOSE (17)
             ENDIF
             !
             r = atoms%rmt(n)*EXP(atoms%dx(n)* (1-atoms%jri(n)))
             dr = EXP(atoms%dx(n))
             DO i = 1,atoms%jri(n)
                atoms%rmsh(i,n) = r
                r = r*dr
             ENDDO
             atoms%volmts(n) = fpi_const/3.*atoms%rmt(n)**3
             cell%volint = cell%volint - atoms%volmts(n)*atoms%neq(n)

             DO n1 = 1,atoms%neq(n)
                na = na + 1
                IF (na>atoms%nat)  CALL juDFT_error("natd too small",calledby ="inped")
                !
                !--->    the in-plane coordinates refer to the lattice vectors a1 and a2,
                !--->    i.e. they are given in internal units scaled by 'scpos'
                !
                WRITE (oUnit,FMT=8170) (atoms%taual(i,na),i=1,3),1.0
8170            FORMAT (4f10.6)
                !
                !--->   for films, the z-coordinates are given in absolute values:
                !
                IF (input%film) atoms%taual(3,na) = scaleCell*atoms%taual(3,na)/a3(3)
                !
                ! Transform intern coordinates to cartesian:
                !
                atoms%pos(:,na)=MATMUL(cell%amat,atoms%taual(:,na))
             ENDDO  ! l.o. equivalent atoms (n1)
          ENDDO     ! loop over atom-types (n)

          IF (input%film .AND. .NOT.oneD%odd%d1) THEN
             !Check if setup is roughly centered
             IF (ABS(MAXVAL(atoms%pos(3,:))+MINVAL(atoms%pos(3,:)))>2.0) &
                  CALL juDFT_warn("Film setup not centered", hint= "The z = 0 plane is the center of the film",calledby="inped")
          ENDIF

          !
          !  check muffin tin radii
          !
          l_gga= xcpot%needs_grad()
          l_test = .TRUE.                  ! only checking, dont use new parameters
          !CALL chkmt(atoms,input,vacuum,cell,oneD,l_test,l_gga,noel, kmax1,dtild,dvac1,lmax1,jri1,rmt1,dx1)

          WRITE (oUnit,FMT=8180) cell%volint
8180      FORMAT (13x,' volume of interstitial region=',f12.6)
          atoms%nat = na
          !--->    evaluate cartesian coordinates of positions
          WRITE (oUnit,FMT=8190) atoms%ntype,atoms%nat
8190      FORMAT (/,/,' number of atom types=',i3,/, ' total number of atoms=',i4,/,/,t3,'no.',t10,'type',&
               &       t21,'int.-coord.',t49,'cart.coord.',t76,'rmt',t84, 'jri',t92,'dx',t98,'lmax',/)
          na = 0
          DO  n = 1,atoms%ntype
             DO n1 = 1,atoms%neq(n)
                na = na + 1
                iz = NINT(atoms%zatom(n))
                WRITE (oUnit,FMT=8200) na,namat_const(iz),n, (atoms%taual(i,na),i=1,3), (atoms%pos(i,na),i=1,3),&
                     atoms%rmt(n),atoms%jri(n), atoms%dx(n),atoms%lmax(n)
8200            FORMAT (1x,i3,4x,a2,t12,i3,2x,3f6.2,3x,3f10.6,3x, f10.6,i6,3x,f6.4,3x,i2)
             ENDDO
          ENDDO
          !
          !
          IF (input%l_useapw) THEN

             DO n = 1,atoms%ntype
                !+APW
                atoms%lapw_l(n) = (atoms%lnonsph(n) - MOD(atoms%lnonsph(n),10) )/10
                atoms%lnonsph(n) = MOD(atoms%lnonsph(n),10)
                !-APW
                IF (atoms%lnonsph(n).EQ.0) atoms%lnonsph(n) = atoms%lmax(n)
                atoms%lnonsph(n) = MIN(atoms%lnonsph(n),atoms%lmax(n))
             ENDDO
          ENDIF

          !--->    nwd = number of energy windows; lepr = 0 (1) for energy
          !--->    parameters given on absolute (floating) scale
          WRITE (oUnit,FMT=8320) input%l_f,input%eonly,1,llr(0)
          WRITE (oUnit,FMT=8330) atoms%ntype, (atoms%lnonsph(n),n=1,atoms%ntype)
8320      FORMAT (1x,/,/,/,' input of parameters for eigenvalues:',/,t5,&
               &       'calculate Pulay-forces = ',l1,/,t5,'eigenvalues ',&
               &       'only = ',l1,/,t5,'number of energy windows =',i2,/,t5,&
               &       'energy parameter mode: ',a8,/,/)
8330      FORMAT (t5,'max. l value in wavefunctions for atom type(s) 1 to',&
               &       i3,':',16i3,/, (t59,16i3,/))
          !
          !--->    input information  for each window
          !
          !--->    energy window

          !--->    for floating energy parameters, the window will be given relative
          !--->    to the highest/lowest energy parameters. a sanity check is made here
          !
          WRITE (oUnit,FMT=8350) input%ellow,input%elup,input%zelec
8350      FORMAT (1x,/,/,' energy window from',f8.3,' to', f8.3,' hartrees; nr. of electrons=',f6.1)
          !--->    input of wavefunction cutoffs: input is a scaled quantity
          !--->    related to the absolute value by rscale (e.g. a muffin-tin
          !--->    radius)
          WRITE (oUnit,FMT=8290) input%rkmax
8290      FORMAT (1x,/,' wavefunction cutoff =',f10.5)
          !
          WRITE (oUnit,FMT=8230) input%bz_integration==1,input%tkb
          WRITE (oUnit,FMT=8240) input%zelec,input%tkb
8230      FORMAT (/,10x,'gauss-integration is used  =',3x,l1,/,10x, 'gaussian half width        =',f10.5)
8240      FORMAT (/,10x,'number of valence electrons=',f10.5,/,10x, 'temperature broadening     =',f10.5)
          WRITE (oUnit,FMT=*) 'itmax=',input%itmax,' broy_sv=',input%maxiter,' imix=',input%imix
          WRITE (oUnit,FMT=*) 'alpha=',input%alpha,' spinf=',input%spinf

          IF ((.NOT.sym%invs).AND.input%secvar) THEN
             WRITE(oUnit,*)'The second variation is not implemented in the'
             WRITE(oUnit,*)'complex version of the program.'
             CALL juDFT_error ("second variation not implemented in complex version" ,calledby ="inped")
          ENDIF

          IF ( (input%jspins.EQ.1).AND.(input%kcrel.EQ.1) )  THEN
             WRITE (oUnit,*) 'WARNING : in a non-spinpolarized calculation the'
             WRITE (oUnit,*) 'coupled-channel relativistic coreprogram (kcrel=1)'
             WRITE (oUnit,*) 'makes no sense; **** setting kcrel = 0 ****'
             input%kcrel = 0
          ENDIF

          WRITE (oUnit,'(a7,l1)') 'swsp = ',input%swsp
          WRITE (oUnit,'(15f6.2)') (atoms%bmu(i),i=1,atoms%ntype)
          IF (vacuum%layers>vacuum%layerd)  CALL juDFT_error("too many layers",calledby ="inped")
          IF (sliceplot%slice) THEN
             input%cdinf = .FALSE.
             WRITE (oUnit,FMT=8390) sliceplot%kk,sliceplot%e1s,sliceplot%e2s
          END IF
8390      FORMAT (' slice: k=',i3,' e1s=',f10.6,' e2s=',f10.6)
          !
          ! Check the LO stuff:
          !
          DO n=1,atoms%ntype
             IF (atoms%nlo(n).GE.1) THEN
                IF (input%secvar)          CALL juDFT_error ("LO + sevcar not implemented",calledby ="inped")
                IF (atoms%nlo(n).GT.atoms%nlod) THEN
                   WRITE (oUnit,*) 'nlo(n) =',atoms%nlo(n),' > nlod =',atoms%nlod
                   CALL juDFT_error("nlo(n)>nlod",calledby ="inped")
                ENDIF
                DO j=1,atoms%nlo(n)
                   IF (.NOT.input%l_useapw) THEN
                      IF (atoms%llo(j,n).LT.0) THEN ! CALL juDFT_error("llo<0 ; compile with DCPP_APW!",calledby="inped")
                         WRITE(oUnit,'(A)') 'Info: l_useapw not set.'
                         WRITE(oUnit,'(A,I2,A,I2,A)') '      LO #',j,' at atom type',n, ' is an e-derivative.'
                      ENDIF
                   ENDIF
                   IF ( (atoms%llo(j,n).GT.atoms%llod).OR.(MOD(-atoms%llod,10)-1).GT.atoms%llod ) THEN
                      WRITE (oUnit,*) 'llo(j,n) =',atoms%llo(j,n),' > llod =',atoms%llod
                      CALL juDFT_error("llo(j,n)>llod",calledby ="inped")
                   ENDIF
                ENDDO
                CALL setlomap(n, input%l_useapw,atoms)
                WRITE (oUnit,*) 'atoms%lapw_l(n) = ',atoms%lapw_l(n)
             ENDIF
          ENDDO
          !
          ! Check for LDA+U:
          !
          DO i = 1, atoms%n_u
             n = atoms%lda_u(i)%atomType
             IF (atoms%nlo(n).GE.1) THEN
                DO j = 1, atoms%nlo(n)
                   IF ((ABS(atoms%llo(j,n)).EQ.atoms%lda_u(i)%l) .AND. (.NOT.atoms%l_dulo(j,n)) ) &
                        WRITE (*,*) 'LO and LDA+U for same l not implemented'
                END DO
             END IF
          END DO
          IF (atoms%n_u.GT.0) THEN
             IF (input%secvar)          CALL juDFT_error ("LDA+U and sevcar not implemented",calledby ="inped")
             IF (noco%l_mperp)         CALL juDFT_error ("LDA+U and l_mperp not implemented",calledby ="inped")
          ENDIF
          !
          !     check all the dos-related switches!
          !
          IF (banddos%ndir.LT.0.AND..NOT.banddos%dos) THEN
             CALL juDFT_error('STOP banddos: the inbuild dos-program  <0 can only be used if dos = true',calledby ="inped")
          ENDIF

          IF (banddos%ndir.LT.0.AND.banddos%dos) THEN
             IF (banddos%e1_dos-banddos%e2_dos.LT.1e-3) THEN
                CALL juDFT_error("STOP banddos: no valid energy window for internal dos-program",calledby ="inped")
             ENDIF
             IF (banddos%sig_dos.LT.0) THEN
                CALL juDFT_error ("STOP DOS: no valid broadening (sig_dos) for internal dos-PROGRAM",calledby ="inped")
             ENDIF
          ENDIF

          IF (banddos%vacdos) THEN
             IF (.NOT. banddos%dos) THEN
                CALL juDFT_error ("STOP DOS: only set vacdos = .true. if dos = .true." ,calledby ="inped")
             ENDIF
             IF (.NOT.vacuum%starcoeff.AND.(vacuum%nstars.NE.1))THEN
                CALL juDFT_error("STOP banddos: if stars = f set vacuum=1" ,calledby ="inped")
             ENDIF
             IF (vacuum%layers.LT.1) THEN
                CALL juDFT_error("STOP DOS: specify layers if vacdos = true" ,calledby ="inped")
             ENDIF
             DO i=1,vacuum%layers
                IF (vacuum%izlay(i,1).LT.1) THEN
                   CALL juDFT_error("STOP DOS: all layers must be at z>0" ,calledby ="inped")

                ENDIF
             ENDDO
          ENDIF

          RETURN
        END SUBROUTINE inped
!--------------------------------------------------------------
      SUBROUTINE trans(string, str_up,str_do)

      IMPLICIT NONE
      CHARACTER(len=2), INTENT(IN)  :: string
      CHARACTER(len=2), INTENT(OUT) :: str_up,str_do

      INTEGER offs,i,n
      CHARACTER(len=2) :: str_in
      CHARACTER(len=1) :: st(2)

      str_up='  ' ; str_do='  ' ; st(:)=' '
      offs = IACHAR('A') - IACHAR('a')
      str_in = TRIM(ADJUSTL(string))
      n = LEN_TRIM(str_in)
      st = (/(str_in(i:i),i=1,n)/)
      DO i=1,n
        IF (IACHAR(st(i)) > IACHAR('Z')) THEN ! lowercase
          str_up(i:i) = CHAR( IACHAR(st(i)) + offs)
        ELSE
          str_up(i:i) = st(i)
        ENDIF
      ENDDO
      DO i=1,n
        str_do(i:i) = CHAR( IACHAR(str_up(i:i)) - offs)
      ENDDO
      END SUBROUTINE trans

      END MODULE
