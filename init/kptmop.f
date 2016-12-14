      MODULE m_kptmop
      use m_juDFT
!-----------------------------------------------------------------------
! ---> This program generates k-points
!           in irreducible wedge of BZ  (for nreg=0)
!           in total BZ                 (for nreg=1)
!      (BZ = 1. Brillouin-zone) for all canonical Bravais lattices
!      in 2 and 3 dimensions,
!      using the basis vectors of the reciprocal lattice
!      and the bordering planes of the irreducible wedge.
!
!      The k-points are generated by the Monkhorst--Pack method.
!      The information on the bordering planes of the irr wedge
!      is taken from BRZONE.
!
!      The program checks the compatibility of the dimension and
!      symmetry and of the provided Monkhorst-Pack-parameters.
!-----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE kptmop(
     >                  iofile,iokpt,kpri,ktest,
     >                  idsyst,idtype,nmop,ikzero,kzero,
     >                  rltv,bltv,nreg,nfulst,nbound,idimens,
     >                  xvec,fnorm,fdist,ncorn,nface,nedge,cpoint,
     >                  nsym,ccr,rlsymr,talfa,mkpt,mface,mdir,
     <                  nkpt,divis,vkxyz,nkstar,wghtkp)
c
c    Meaning of variables:
c    INPUT:
c
c    Symmetry of lattice:
c    idsyst   : crystal system identification in MDDFT programs
c    idtype   : lattice type identification in MDDFT programs
c    bltv     : cartesian coordinates of basis vectors for
c               Bravais lattice: bltv(ix,jn), ix=1,3; jn=1,3
c    rltv     : cartesian coordinates of basis vectors for
c               reciprocal lattice: rltv(ix,jn), ix=1,3; jn=1,3
c    nsym     : number of symmetry elements of points group
c    ccr     : rotation matrix for symmetry element
c                   in cartesian representation
c    rlsymr   : rotation matrix for symmetry element
c                   in reciprocal lattice basis representation
c    talfa    : translation vector associated with (non-symmorphic)
c               symmetry elements in Bravais lattice representation
c
c    representation of the irreducible part of the BZ:
c    fnorm    : normal vector of the planes bordering the irrBZ
c    fdist    : distance vector of the planes bordering the irrBZ
c    ncorn    : number of corners of the irrBZ
c    nface    : number of faces of the irrBZ
c    nedge    : number of edges of the irrBZ
c    xvec     : arbitrary vector lying in the irrBZ (FOR SURE!!)
c               components are:
c
c    characterization of the Monkhorst-Pack k-point set:
c    idimens  : number of dimensions for k-point set (2 or 3)
c    nreg     : 1 kpoints in full BZ; 0 kpoints in irrBZ
c    nfulst   : 1 kpoints ordered in full stars 
c                  (meaningful only for nreg =1; full BZ)
c    nbound   : 0 no primary points on BZ boundary; 
c               1 with boundary points (not for BZ integration!!!)
c    ikzero   : 0 no shift of k-points; 
c               1 shift of k-points for better use of sym in irrBZ
c    kzero    : shifting vector to bring one k-point to (0,0,0)(for even nmop)
c                                             away from (0,0,0)(for odd nmop)  
c    nmop     : integer number triple: nmop(i), i=1,3; nmop(i)
c               determines number of k-points in direction of rltv(ix,i)
c
c    OUTPUT: k-point set
c    nkpt     : number of k-points generated in set
c    vkxyz    : vector of kpoint generated; in cartesian representation
c    wghtkp   : weight associated with k-points for BZ integration
c    divis    : integer triple divis(i); i=1,4.
c               Used to find more accurate representation of k-points
c               vklmn(i,kpt)/divis(i) and weights as wght(kpt)/divis(4)
c    nkstar   : number of stars for k-points generated in full stars
c
c-----------------------------------------------------------------------
      USE m_constants, ONLY : pimach
      USE m_ordstar
      USE m_fulstar
      IMPLICIT NONE
C
C-----> PARAMETER STATEMENTS
C
      INTEGER, INTENT (IN) :: mkpt,mface,mdir
c
c ---> file number for read and write
c
      INTEGER, INTENT (IN) :: iofile,iokpt
c
c ---> running mode parameter
c
      INTEGER, INTENT (IN) :: kpri,ktest
C
C----->  Symmetry information
C
      INTEGER, INTENT (IN) :: nsym,idsyst,idtype
      REAL,    INTENT (IN) :: ccr(3,3,48)
      REAL,    INTENT (IN) :: rlsymr(3,3,48),talfa(3,48)
C
C----->  BRAVAIS LATTICE INFORMATION
C
      REAL,    INTENT (IN) :: bltv(3,3),cpoint(3,mface)
C
C----->  RECIPROCAL LATTICE INFORMATION
C
      INTEGER, INTENT (IN) :: ncorn,nface,nedge
      REAL,    INTENT (IN) :: xvec(3),rltv(3,3)
      REAL,    INTENT (IN) :: fnorm(3,mface),fdist(mface)
C
C----->  BRILLOUINE ZONE INTEGRATION
C
      INTEGER, INTENT (IN) :: nreg,ikzero,nfulst,nbound,idimens
      INTEGER, INTENT (INOUT) :: nmop(3)
      REAL,    INTENT (INOUT) :: kzero(3)
      INTEGER, INTENT (OUT):: nkpt,nkstar
      REAL,    INTENT (OUT):: vkxyz(3,mkpt),wghtkp(mkpt),divis(4)
C
C --->  local variables
c
      CHARACTER*80 blank
      INTEGER  i,idim,i1,i2,i3,ii,ij,ik,is,isym,ifac, iik,iiik
      INTEGER  ikc, i1red,nred,isumkpt,nleft,nirrbz
      INTEGER  dirmin,dirmax,ndir1,ndir2,idir
      INTEGER  kpl,kpm,kpn,nstar(mdir),nstnew
      INTEGER  iplus,iminus,nc2d,n
      REAL     invtpi,zero,one,half,eps,eps1,orient,sum,denom,aivnkpt

      INTEGER  nfract(3),lim(3),isi(3)
      REAL     cp2d(3,mface)
      REAL     vktra(3),vkstar(3,48),ainvnmop(3),fsig(2),vktes(3)
      INTEGER, ALLOCATABLE :: ikpn(:,:),irrkpn(:),nkrep(:)
      INTEGER, ALLOCATABLE :: iside(:),iostar(:)
      REAL,    ALLOCATABLE :: fract(:,:), vkrep(:,:)
C
C --->  intrinsic functions
c
      INTRINSIC   real,abs
C
C --->  save and data statements
c
      SAVE     one,zero,half,eps,eps1,iplus,iminus
      DATA     zero/0.00/,one/1.00/,half/0.50/,
     +         eps/1.0e-8/,eps1/1.0e-6/,iplus/1/,iminus/-1/
c
c-----------------------------------------------------------------------
c
      ALLOCATE (fract(mkpt,3),vkrep(3,mkpt),ikpn(48,mkpt),irrkpn(mkpt))
      ALLOCATE (nkrep(mkpt),iostar(mkpt),iside(mface))
      if (kpri .ge. 1) then
c       write(iofile,'(/)')
        write(iofile,'(3x,'' *<* kptmop *>* '')')
        write(iofile,'(3x,'' generate k-vectors'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~'')')
        write(iofile,'(3x,'' by Monkhorst-Pack-method'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~'')')
        if (idimens .eq. 2) then
            write(iofile,'(3x,'' in 2 dimensions'')')
            write(iofile,'(3x,'' ~~~~~~~~~~~~~~~'')')
            write(iokpt,'(''#'',i4,20x,'' idimens: '',
     >      ''k-points for 2-dimensional lattice'')') idimens
        else if (idimens .lt. 2 .or. idimens.gt.3) then
            write(iofile,'(1x,i4,'' idimens wrong, choose 2 or 3'')')
     >                                                       idimens
             CALL juDFT_error("idimens",calledby="kptmop")
        end if
        if (nreg .eq. 1) then
           if (nfulst .eq. 1) then
            write(iofile,'(3x,'' in 1. Brillouin zone'')')
            write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~'')')
            write(iofile,'(3x,'' full stars generated'')')
            write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~'')')
            write(iokpt,'(''#'',2(i4,1x),14x,'' nreg,nfulst: '',
     >      ''k-points in totBZ, ordered in full stars'')') nreg,nfulst
           else if (nfulst .eq. 0) then
            write(iofile,'(3x,'' in parallelepiped'')')
            write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~'')')
            write(iofile,'(3x,'' full stars not generated'')')
            write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~'')')
            write(iokpt,'(''#'',2(i4,1x),14x,'' nreg,nfulst: '',
     >      ''k-points in par-epiped, not in full stars'')') nreg,nfulst
           else
            write(iofile,'(2(1x,i4),15x,'' nreg,nfulst: wrong choice '',
     >       ''of parameters; allowed combinations: (1,1); (1,0)'')' )
     >                                                      nreg,nfulst
              CALL juDFT_error("nfulst",calledby="kptmop")
           end if
        else if (nreg .eq. 0) then
        write(iofile,'(3x,'' in irred wedge of 1. Brillouin zone'')')
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')')
        write(iokpt,'(''#'',i4,21x,''nreg: k-points in irrBZ'')') nreg
 
        else
        write(iofile,'(3x,'' wrong choice of nreg: '', i4)') nreg
 
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')')
         CALL juDFT_error("nreg",calledby="kptmop")
        end if
        write(iokpt,'(''#'',3(i4,1x),10x,''nmop(1), nmop(2), nmop(3)'',
     +     '' as read in'')')                      (nmop(i), i=1,3)
        write(iokpt,'(''#'',2(i4,1x),14x,'' ikzero, nbound'')')
     +                                               ikzero, nbound
        write(iofile,'(/)')
        if (nbound .eq. 1) then
        write(iofile,'(3x,'' k-points on boundary included'')')
        write(iofile,'(3x,'' irregular Monkhorst-Pack-ratios'')')
        write(iofile,'(3x,'' cannot be used for BZ-integration'')')
        else if (nbound.eq. 0) then
        write(iofile,'(3x,'' no k-points on boundary of BZ'')')
        write(iofile,'(3x,'' regular Monkhorst-Pack-ratios'')')
        write(iofile,'(3x,'' can be used for BZ-integration'')')
        else
        write(iofile,'(3x,'' wrong choice of nbound: '', i4)') nbound
        write(iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')')
         CALL juDFT_error("nbound",calledby="kptmop")
        end if
c
        write(iofile,'(1x,i4,10x,''iofile'')') iofile
        write(iofile,'(1x,i4,10x,''iokpt'')')  iokpt
        write(iofile,'(1x,i4,10x,''kpri'')')  kpri
        write(iofile,'(1x,i4,10x,''ktest'')')  ktest
        write(iofile,'(1x,3(f10.7,1x),10x,''xvec'')') (xvec(ii),ii=1,3)
        write(iofile,'(1x,i4,10x,''ncorn'')')  ncorn
        write(iofile,'(1x,i4,10x,''nedge'')')  nedge
        write(iofile,'(1x,i4,10x,''nface'')')  nface
        do  5 ifac = 1,nface
        write(iofile,'(1x,i4,1x,3(f10.7,1x),10x,''fnorm'')')
     +                 ifac,(fnorm(ii,ifac),ii=1,3)
        write(iofile,'(1x,i4,1x,f10.7,1x,10x,''fdist'')')
     +                 ifac,fdist(ifac)
  5     continue
c
      end if
!
! --->   for 2 dimensions only the following Bravais lattices exist:
!          TYPE                    EQUIVALENT 3-DIM        idsyst/idtype
!         square               = p-tetragonal ( 1+2 axis )      2/1
!         rectangular          = p-orthorhomb ( 1+2 axis )      3/1
!         centered rectangular = c-face-orthorhomb( 1+2 axis)   3/6
!         hexagonal            = p-hexagonal  ( 1+2 axis )      4/1
!         oblique              = p-monoclinic ( 1+2 axis )      6/1
!
      IF (idimens .EQ. 2) THEN
!
! --->   identify the allowed symmetries
!        and check the consistency of the Monkhorst-Pack-parameters
!
        IF (idsyst.EQ.2 .OR. idsyst.EQ.4) THEN
          IF (idtype.EQ.1) THEN
             IF (nmop(1).NE.nmop(2) .OR. nmop(3).NE.0) THEN
                nmop(2) = nmop(1)
                nmop(3) = 0
                WRITE (iofile,'(1x,''WARNING!!!!!!!'',/,
     +            ''nmop-Parameters not in accordance with symmetry'',/,
     +            2(1x,i4),/,
     +            '' we have set nmop(2) = nmop(1)'',/,
     +            '' and/or nmop(3) = 0'')') idsyst, idtype
                WRITE (iofile,'(3(1x,i4),'' new val for nmop: '')')
     +                                              (nmop(i),i=1,3)
             ELSE
                WRITE (iofile,'('' values accepted unchanged'')')
                WRITE (iofile,'(3(1x,i4),14x,''nmop(i),i=1,3'')')
     +                                            (nmop(i),i=1,3)
             ENDIF
          ENDIF
        ELSEIF (idsyst.EQ.3) THEN
          IF (idtype.EQ.1 .OR. idtype.EQ.6) THEN
             IF (nmop(3).NE.0) THEN
                nmop(3) = 0
                WRITE (iofile,'(1x,''WARNING!!!!!!!'',/,
     +            ''nmop-Parameters not in accordance with symmetry'',/,
     +            2(1x,i4),/,
     +            '' we have set nmop(3) = 0'')') idsyst, idtype
                WRITE (iofile,'(3(1x,i4),'' new val for nmop: '')')
     +                                              (nmop(i),i=1,3)
             ELSE
                WRITE (iofile,'('' values accepted unchanged'')')
                WRITE (iofile,'(3(1x,i4),14x,''nmop(i),i=1,3'')')
     +                                            (nmop(i),i=1,3)
             ENDIF
          ENDIF
        ELSEIF (idsyst.EQ.6) THEN
          IF (idtype.EQ.1) THEN
             IF (nmop(3).NE.0) THEN
                nmop(3) = 0
                WRITE (iofile,'(1x,''WARNING!!!!!!!'',/,
     +            ''nmop-Parameters not in accordance with symmetry'',/,
     +            2(1x,i4),/,
     +            '' we have set nmop(3) = 0'')') idsyst, idtype
                WRITE (iofile,'(3(1x,i4),'' new val for nmop: '')')
     +                                              (nmop(i),i=1,3)
             ELSE
                WRITE (iofile,'('' values accepted unchanged'')')
                WRITE (iofile,'(3(1x,i4),14x,''nmop(i),i=1,3'')')
     +                                             (nmop(i),i=1,3)
             ENDIF
          ENDIF
        ELSE
!
! --->   in all other cases:
!
          WRITE (iofile,'(3(1x,i4),20x,'' idimens,idsyst,idtype: '',
     >      ''wrong choice for 2-dimensional crystal structure'')')
     >                                 idimens,idsyst,idtype
           CALL juDFT_error("2-dim crystal",calledby="kptmop")
        ENDIF 
!
! --->   check consistency of nmop-parameters with crystal symmetry
!
      ELSEIF (idimens.EQ.3) THEN
         IF (idsyst.EQ.1 .OR. idsyst.EQ.5) THEN
             IF (nmop(1).NE.nmop(2) .OR. nmop(1).NE.nmop(3)
     +                              .OR. nmop(2).NE.nmop(3)) THEN
                nmop(3) = nmop(1)
                nmop(2) = nmop(1)
                WRITE (iofile,'(1x,''WARNING!!!!!!!'',/,
     +          ''nmop-Parameters not in accordance with symmetry'',/,
     +          2(1x,i4),/,
     +          '' we have set all nmop(i) = nmop(1)'')') idsyst, idtype
                WRITE (iofile,'(3(1x,i4),'' new val for nmop(i): '')')
     +                                                 (nmop(i),i=1,3)
             ELSE
                WRITE (iofile,'('' values accepted unchanged'')')
                WRITE (iofile,'(3(1x,i4),14x,''nmop(i),i=1,3'')')
     +                                                 (nmop(i),i=1,3)
             ENDIF
         ELSEIF (idsyst.EQ.2 .OR. idsyst.eq.4) THEN
             if((nmop(3).eq.nmop(2)).and.idsyst.eq.2)then
                WRITE (iofile,'('' values accepted unchanged'')')
                WRITE (iofile,'(3(1x,i4),14x,''nmop(i),i=1,3'')')
     +                                            (nmop(i),i=1,3)
             elseif(nmop(1).NE.nmop(2)) THEN
                nmop(2) = nmop(1)
                WRITE (iofile,'(1x,''WARNING!!!!!!!'',/,
     +            ''nmop-Parameters not in accordance with symmetry'',/,
     +            2(1x,i4),/,
     +            '' we have set nmop(2) = nmop(1)'')') idsyst, idtype
                WRITE (iofile,'(3(1x,i4),'' new val for nmop: '')')
     +                                              (nmop(i),i=1,3)
             ELSE
                WRITE (iofile,'('' values accepted unchanged'')')
                WRITE (iofile,'(3(1x,i4),14x,''nmop(i),i=1,3'')')
     +                                            (nmop(i),i=1,3)
             ENDIF
         ELSEIF (idsyst.LT.1 .OR. idsyst.GT.7) THEN
             WRITE (iofile,'(1x,''wrong choice of symmetry'',/,
     +                               2(1x,i4))') idsyst, idtype
             WRITE (iofile,'(''only values 1.le.idsyst.le.7 allowed'')')
              CALL juDFT_error("wrong idsyst",calledby="kptmop")
         ELSE
             WRITE (iofile,'('' values accepted unchanged'')')
             WRITE (iofile,'(3(1x,i4),11x,''nmop(i),i=1,3'')')
     +                                          (nmop(i),i=1,3)
         ENDIF
      ELSE
          CALL juDFT_error("idimens =/= 2,3 ",calledby="kptmop")
      ENDIF
!
! --->   start calculation
! =====================================================================
!
! ---> set sign constants
       isi(1) = 0
       isi(2) = iminus
       isi(3) = iplus
!
! ---> calc orientation of boundary faces of irr wedge of BZ
!      characterized by
!          iside(i)= sign( (xvec,fnorm(i))-fdist(i) ) ;(i=1,nface )
!
      WRITE (iofile,'(1x,''orientation of boundary faces'')')
      DO ifac = 1, nface
         orient = zero
         iside(ifac) = iplus
         DO ii = 1, 3
            orient = orient + xvec(ii)*fnorm(ii,ifac)
         ENDDO
         orient = orient - fdist(ifac)
         IF (orient .LT. 0) iside(ifac) = iminus
         WRITE (iofile,'(1x,2(i4,2x),f10.7,10x,''ifac,iside,orient'',
     +                       '' for xvec'')') ifac,iside(ifac),orient
      ENDDO

      invtpi = one / ( 2.0 * pimach() )

      WRITE (iofile,'(''Bravais lattice vectors'')' )
      DO ii = 1, 3
         WRITE (iofile,'(43x,3(1x,f11.6))') (bltv(ii,ikc), ikc=1,3)
      ENDDO
      WRITE (iofile,'(''reciprocal lattice vectors'')' )
      DO ii = 1, 3
         WRITE (iofile,'(43x,3(1x,f11.6))' ) (rltv(ii,ikc), ikc=1,3)
      ENDDO
!
! ---> nmop(i) are Monkhorst-Pack parameters; they determine the
!                           number of k-points in i-direction
!      if basis vector lengths are not related by symmetry,
!      we can use independent fractions for each direction
!
      WRITE (iofile,'(3(1x,i4),10x,'' Monkhorst-Pack-parameters'')')
     +                                             (nmop(i1),i1=1,3)

      DO idim = 1, idimens
         IF (nmop(idim).GT.0) THEN
            ainvnmop(idim) = one/ real(nmop(idim))
         ELSE
            WRITE (iofile,'('' nmop('',i4,'') ='',i4,
     +                     '' not allowed'')') idim, nmop(idim)
             CALL juDFT_error("nmop wrong",calledby="kptmop")
         ENDIF
      ENDDO
!
! --->nreg determines region of BZ in which k-points are generated
!
      IF ( nreg .EQ. 1) THEN
        IF (nfulst .EQ. 1) THEN
          WRITE (iofile,'(2(1x,i4),15x,'' nreg,nfulst: '',
     >      ''k-points in totBZ, ordered in full stars'')') nreg,nfulst
        ELSEIF (nfulst .EQ. 0) then
          WRITE (iofile,'(2(1x,i4),15x,'' nreg,nfulst: '',
     >     ''k-points in par-epiped, not in full stars'')') nreg,nfulst
        ENDIF
      ELSEIF ( nreg .EQ. 0) THEN
        WRITE (iofile,'(1x,i4,10x,''nreg; k-points in'',
     +                '' irreducible wedge of BZ'')' ) nreg
      ENDIF

      WRITE (iofile,'(1x,''Monkhorst-Pack-fractions'')' )
!
! ---> nbound=1: k-points are generated on boundary of BZ
!        include  fract(1) =       -1/2
!            and  fract(2*nmop+1) = 1/2     for surface points of BZ
!
      IF ( nbound .EQ. 1) THEN
        WRITE (iofile,'(1x,i4,10x,''nbound; k-points on boundary'',
     +   '' of BZ included'')' ) nbound
!
! ---> irregular Monkhorst--Pack--fractions
!                                   fract(r) = r / (2*nmop)
!
        DO idim = 1,idimens
           denom = half*ainvnmop(idim)
           divis(idim) = one / denom

           DO kpn = -nmop(idim),nmop(idim)
             fract(kpn+nmop(idim)+1,idim) = denom * real (kpn)
             WRITE (iofile,'(10x,f10.7)' ) fract(kpn+nmop(idim)+1,idim)
           ENDDO
           nfract(idim) = 2*nmop(idim) + 1
        ENDDO
        IF (idimens .eq. 2) THEN
           nfract(3) = 1
           fract(1,3) = 0
           divis(3) = one
        END IF
!
! ---> nbound=0: k-points are NOT generated on boundary of BZ
!                This is the regular Monkhorst-Pack-method
!
      ELSEIF ( nbound .eq. 0) then
         WRITE (iofile,'(1x,i4,10x,''nbound; no k-points '',
     +                 '' on boundary of BZ'')' ) nbound
!
! --->   regular Monkhorst--Pack--fractions
!                                   fract(r) =(2*r-nmop-1) / (2*nmop)
!
         DO idim = 1,idimens
            denom = half*ainvnmop(idim)
            divis(idim) = one / denom
            WRITE(iofile,'(5x,i4,5x,''idim'')' ) idim
            DO  kpn = 1,nmop(idim)
              fract(kpn,idim) = denom * real (2*kpn -nmop(idim)-1)
              write(iofile,'(10x,f10.7)' ) fract(kpn,idim)
            ENDDO
            nfract(idim) = nmop(idim)
         ENDDO

         IF (idimens .EQ. 2) THEN
            nfract(3) = 1
            fract(1,3) = 0
            divis(3) = one
         ENDIF

      ELSE
         WRITE (iofile,'(3x,'' wrong choice of nbound:'', i4)') nbound
         WRITE (iofile,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')')
          CALL juDFT_error("nbound",calledby="kptmop")
      ENDIF
!
! ---> set kzero
!
      kzero(1) = zero
      kzero(2) = zero
      kzero(3) = zero

      IF (ikzero.NE.0) THEN
!
! --->   set kzero .ne. zero
!        for nmop even and non-orthogonal  rec lattice vectors
!          shift all k-points to bring one to (0,0,0) 
!        for nmop odd and non-orthogonal  rec lattice vectors
!          shift k-points to avoid (0,0,0) 
!            implemented for hexagon   (idsyst =4, idtype arbitrary)
!                        and ort(2-dim)(idsyst =3, idtype =1)
!                        and fcc       (idsyst =1, idtype =3)
!                        and p-monocli (idsyst =6, idtype =1)
!
        denom = half*ainvnmop(1)

        IF( idsyst .EQ. 4) THEN    !  hexagonal Bravais lattice

          DO ii = 1, 3
            kzero(ii) = denom * (rltv(ii,1)+rltv(ii,2))
          ENDDO

        ELSEIF ( idsyst .EQ. 3 .and. idtype .EQ. 1) THEN ! orthorhombic (2D)

          IF (nmop(1).EQ.1) THEN
              denom = half*ainvnmop(3)
            DO ii = 1, 3
              kzero(ii) = denom * (rltv(ii,2)+rltv(ii,3))
            ENDDO
          ELSE
            WRITE (iofile,'(''ikzero.ne.0 is NOT applied '',
     >                     ''for present choice of parameters'')')
            WRITE (iofile,'(5(1x,i4),
     >               ''; idsyst,idtype,nmop(1),nmop(2),nmop(3)'')')
     >                       idsyst,idtype,nmop(1),nmop(2),nmop(3)
          ENDIF

        ELSEIF ( idsyst .eq. 1 .and. idtype .eq. 3) then !  face centered cubic

          DO ii = 1, 3
            kzero(ii) = denom * (rltv(ii,1)+rltv(ii,2)+rltv(ii,3))
          ENDDO

        ELSEIF ( idsyst .eq. 6 .and. idtype .eq. 1) then ! p-monoclinic 

          DO ii = 1, 3
            kzero(ii) = half*ainvnmop(1) * rltv(ii,1)
     +                + half*ainvnmop(2) * rltv(ii,2)
          ENDDO

        ELSE
          WRITE (iofile,'(''ikzero.ne.0 is NOT applied '',
     >                   ''for this system'')')
          WRITE (iofile,'(5(1x,i4),
     >             ''; idsyst,idtype,nmop(1),nmop(2),nmop(3)'')')
     >                     idsyst,idtype,nmop(1),nmop(2),nmop(3)
        ENDIF
      ENDIF
!
!
! --->   initialize k-points = zero and weights = 1.0
!
      DO  kpn = 1,mkpt
         vkxyz(1,kpn) = zero
         vkxyz(2,kpn) = zero
         vkxyz(3,kpn) = zero
         wghtkp(kpn) = one
      ENDDO
!
! ---> generate equidistant k-vectors in cartesian coordinates
!                                   with off-set kzero
!
      nkpt = 0
      DO i3 = 1,nfract(3)
         DO i2 = 1,nfract(2)
            DO i1 = 1,nfract(1)
              nkpt = nkpt + 1
              IF (nkpt>mkpt)  CALL juDFT_error("nkpt > mkpt",calledby
     +             ="kptmop")
              vkxyz(1,nkpt) = kzero(1) + rltv(1,1)*fract(i1,1)
     +                                 + rltv(1,2)*fract(i2,2)
     +                                 + rltv(1,3)*fract(i3,3)
              vkxyz(2,nkpt) = kzero(2) + rltv(2,1)*fract(i1,1)
     +                                 + rltv(2,2)*fract(i2,2)
     +                                 + rltv(2,3)*fract(i3,3)
              vkxyz(3,nkpt) = kzero(3) + rltv(3,1)*fract(i1,1)
     +                                 + rltv(3,2)*fract(i2,2)
     +                                 + rltv(3,3)*fract(i3,3)
            ENDDO
         ENDDO
      ENDDO
!
! --->   calculate weights of k-points and print out k-points
!         wghtkp = 1/nkpt
!              ( = 1/(nmop(1)*nmop(2)*nmop(3)) for reg Monk-Pack-method)
!
      divis(4) = real(nkpt)
      aivnkpt  = one/real(nkpt)

      DO  kpn= 1,nkpt
        wghtkp(kpn) = wghtkp(kpn)*aivnkpt
      ENDDO

      IF (kpri .ge. 1) THEN
      WRITE (iofile,'(''generated k-points vkxyz'')')
      DO  kpn= 1,nkpt
        WRITE (iofile,'(1x,i4,1x,4(f13.10,1x),10x,''vkxyz, wghtkp'')')
     +               kpn,(vkxyz(ii,kpn),ii=1,3),wghtkp(kpn)
      ENDDO
      ENDIF

!
! ---> for nreg=1 and nfulst=0
!                 (k-points in total zone of arbitrary shape)
!       k-point generation finished
!       (k-points are actually generated in parallel-epiped
!       (r1,r2,r3) which is an equivalent unit cell)
!       The generated set does not contain all points of the
!       symmetry related stars.
!
      IF (nreg.EQ.1 .AND. nfulst.EQ.0) GOTO 1000
!
! ====================================================================
!
! --->   order generated k-points in stars by applying symmetry:
!        - determine number of different stars nkstar .le. nkpt
!        - determine order of star iostar(kpn) .le. nsym
!        - assign pointer ikpn(i,ik); i=1,iostar(ik); ik=1,nkstar
!        - determine representative vector in irrBZ for each star:
!                                   vkrep(ix,ik); ix=1,3; ik=1,nkstar
!
        CALL ordstar(
     >               iokpt,kpri,ktest,
     >               fnorm,fdist,nface,iside,
     >               nsym,ccr,rltv,mkpt,mface,mdir,
     =               nkpt,vkxyz,
     <               nkstar,iostar,ikpn,vkrep,nkrep)
!
! ---> for nreg=0 (k-points in irr part of BZ)
!
!      (a) calculate weights for k-points in irrBZ
!           - wghtkp(ik)=iostar(ik)/nkpt_old ; ik=1,nkstar
!
        DO ik = 1, nkstar
             wghtkp(ik) = wghtkp(ik)*iostar(ik)
        ENDDO
!
!      (b) final preparation of k-points for transfer to file
!           - assign nkpt= nkstar
!           - assign vkxyz(ix,ik) = vkrep(ix,ik); ix=1,3; ik=1,nkstar
!
        IF (nreg.EQ.0) THEN

          DO i1 = 1,3
            DO ik = 1,nkstar
              vkxyz(i1,ik) = vkrep(i1,ik)
            ENDDO
          ENDDO
          nkpt = nkstar
          GOTO 1000
        ENDIF
c
c ---> for nreg.eq.0: k-point generation finished
c
c ================================================================
c
c
       IF (nreg.EQ.1 .AND. nfulst.EQ.1) THEN
c
c --->   generate full stars for all representative k-points
c        - for nreg=1 and nfulst=1:
c              - determine order of full star ifstar(kpn).le.nsym
c              - assign nkpt= sum {ifstar(ik)} (ik=1,nkstar)
c              - assign vkxyz(ix,kpn) = vkstar(ix,ikpn(is,ik));
c                     ix=1,3; kpn=1,nkpt; ik=1,nstar; is=1,ifstar(ik)
c              - calculate wghtkp(kpn)=wghtkp_old(ik)/(ifstar(ik)
c                                kpn=1,nkpt; ik=1,nstar
c
c
            CALL fulstar(
     >                   iofile,iokpt,kpri,ktest,
     >                   ccr,nsym,
     >                   vkrep,nkstar,mkpt,mface,mdir,
     =                   nkpt,vkxyz,wghtkp)

            divis(4) = divis(4) * nsym

      ENDIF

 1000 CONTINUE
!      
! --> check for corner points, include them into k-point set:
!
      IF (nbound.EQ.1) THEN
        n = 1
        nc2d = 1               ! determine 2D corner points
        cp2d(:,nc2d) = cpoint(:,n)
        corn: DO n = 2, ncorn
          DO i = 1, n-1
            IF ((abs(cpoint(1,n)-cpoint(1,i)).LT.0.0001).AND.
     +          (abs(cpoint(2,n)-cpoint(2,i)).LT.0.0001)) CYCLE corn
          ENDDO
          nc2d = nc2d + 1
          cp2d(:,nc2d) = cpoint(:,n)
        ENDDO corn
        WRITE (iofile,'(''2D corner points in internal units'')')
        corn2d: DO n = 1, nc2d 
          WRITE (iofile,'(i3,3x,2(f10.7,1x))') n,cp2d(1,n),cp2d(2,n)
          DO i = 1, nkpt
            IF ((abs(cp2d(1,n)-vkxyz(1,i)).LT.0.0001).AND.
     +          (abs(cp2d(2,n)-vkxyz(2,i)).LT.0.0001)) CYCLE corn2d
          ENDDO
          nkpt = nkpt + 1
          vkxyz(:,nkpt) = cp2d(:,n)
        ENDDO corn2d
      ENDIF 
!
! --->   print out k-points and weights
!
      IF (ktest .GE. 1) THEN
        WRITE (iofile,'(''generated k-points vkxyz'')')
        WRITE (iofile,'(1x,i4,20x,
     +               ''nkpt,number of generated k-points'')') nkpt
 
       DO kpn = 1, nkpt
         WRITE (iofile,'(3(f10.7,1x),f12.10,1x,i4,3x,
     +   ''vkxyz, wghtkp'')') (vkxyz(ii,kpn),ii=1,3),wghtkp(kpn), kpn
       ENDDO

      ENDIF
      DEALLOCATE (fract,vkrep,ikpn,irrkpn,nkrep,iostar,iside)

      RETURN
      END SUBROUTINE kptmop
      END MODULE m_kptmop
