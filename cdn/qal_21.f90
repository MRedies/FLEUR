MODULE m_qal21
  !***********************************************************************
  ! Calculates qal21  needed to determine the off-diagonal parts of the
  ! DOS
  !***********************************************************************
  !
CONTAINS
  SUBROUTINE qal_21(atoms,input,noccbd,ev_list,nococonv,eigVecCoeffs,denCoeffsOffdiag,ikpt,dos)
    use m_types_nococonv
    USE m_types_setup
    USE m_types_dos
    USE m_types_cdnval, ONLY: t_eigVecCoeffs
    USE m_types_denCoeffsOffdiag
    USE m_rotdenmat
    use m_constants
    IMPLICIT NONE

    TYPE(t_input),             INTENT(IN)    :: input
    TYPE(t_nococonv),          INTENT(IN)    :: nococonv
    TYPE(t_atoms),             INTENT(IN)    :: atoms
    TYPE(t_eigVecCoeffs),      INTENT(IN)    :: eigVecCoeffs
    TYPE(t_denCoeffsOffdiag),  INTENT(IN)    :: denCoeffsOffdiag
    TYPE(t_dos),               INTENT(INOUT) :: dos

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: noccbd,ikpt

    INTEGER, INTENT (IN) :: ev_list(noccbd)

    !     .. Local Scalars ..
    INTEGER i,l,lo,lop ,natom,nn,ntyp
    INTEGER nt1,nt2,lm,n,ll1,ipol,icore,index,m
    REAL fac
    COMPLEX sumaa,sumbb,sumab,sumba

    !     .. Local Arrays ..
    COMPLEX qlo(noccbd,atoms%nlod,atoms%nlod,atoms%ntype)
    COMPLEX qaclo(noccbd,atoms%nlod,atoms%ntype),qbclo(noccbd,atoms%nlod,atoms%ntype)
    COMPLEX qcloa(noccbd,atoms%nlod,atoms%ntype),qclob(noccbd,atoms%nlod,atoms%ntype)
    COMPLEX qal21(0:3,atoms%ntype,input%neig)
    COMPLEX q_loc(2,2),q_hlp(2,2),chi(2,2)
    REAL    qmat(0:3,atoms%ntype,input%neig,4)

    !     .. Intrinsic Functions ..
    INTRINSIC conjg
    qal21=0.0
    !--->    l-decomposed density for each occupied state
    states : DO i = 1, noccbd
       nt1 = 1
       types_loop : DO n = 1 ,atoms%ntype
          nt2 = nt1 + atoms%neq(n) - 1
          ls : DO l = 0,3
             IF (i==1) THEN
             ENDIF
             sumaa = CMPLX(0.,0.) ; sumab = CMPLX(0.,0.)
             sumbb = CMPLX(0.,0.) ; sumba = CMPLX(0.,0.)
             ll1 = l* (l+1)
             ms : DO m = -l,l
                lm = ll1 + m
                atoms_loop : DO natom = nt1,nt2
                   sumaa = sumaa + eigVecCoeffs%acof(i,lm,natom,1)* CONJG(eigVecCoeffs%acof(i,lm,natom,input%jspins))
                   sumbb = sumbb + eigVecCoeffs%bcof(i,lm,natom,1)* CONJG(eigVecCoeffs%bcof(i,lm,natom,input%jspins))
                   sumba = sumba + eigVecCoeffs%acof(i,lm,natom,1) * CONJG(eigVecCoeffs%bcof(i,lm,natom,input%jspins))
                   sumab = sumab + eigVecCoeffs%bcof(i,lm,natom,1) * CONJG(eigVecCoeffs%acof(i,lm,natom,input%jspins))
                ENDDO atoms_loop
             ENDDO ms
             qal21(l,n,i) = sumaa * denCoeffsOffdiag%uu21n(l,n) + sumbb * denCoeffsOffdiag%dd21n(l,n) +&
                            sumba * denCoeffsOffdiag%du21n(l,n) + sumab * denCoeffsOffdiag%ud21n(l,n)
          ENDDO ls
          nt1 = nt1 + atoms%neq(n)
       ENDDO types_loop
    ENDDO states

    !---> initialize qlo

    qlo(:,:,:,:) = CMPLX(0.,0.)
    qaclo(:,:,:) = CMPLX(0.,0.)
    qcloa(:,:,:) = CMPLX(0.,0.)
    qclob(:,:,:) = CMPLX(0.,0.)
    qbclo(:,:,:) = CMPLX(0.,0.)

    !---> density for each local orbital and occupied state

    natom = 0
    DO ntyp = 1,atoms%ntype
       DO nn = 1,atoms%neq(ntyp)
          natom = natom + 1
          DO lo = 1,atoms%nlo(ntyp)
             l = atoms%llo(lo,ntyp)
             ll1 = l* (l+1)
             DO m = -l,l
                lm = ll1 + m
                DO i = 1, noccbd
                   qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +      &
                        eigVecCoeffs%bcof(i,lm,natom,1)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,input%jspins))
                   qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +      &
                        eigVecCoeffs%ccof(m,i,lo,natom,1)*CONJG(eigVecCoeffs%bcof(i,lm,natom,input%jspins))
                   qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +       &
                        eigVecCoeffs%acof(i,lm,natom,1)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,input%jspins))
                   qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +       &
                        eigVecCoeffs%ccof(m,i,lo,natom,1)*CONJG(eigVecCoeffs%acof(i,lm,natom,input%jspins))
                ENDDO
             ENDDO
             DO lop = 1,atoms%nlo(ntyp)
                IF (atoms%llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1, noccbd
                         qlo(i,lop,lo,ntyp) = qlo(i,lop,lo,ntyp) +  &
                              CONJG(eigVecCoeffs%ccof(m,i,lop,natom,input%jspins))*eigVecCoeffs%ccof(m,i,lo,natom,1) +&
                              CONJG(eigVecCoeffs%ccof(m,i,lo,natom,input%jspins))*eigVecCoeffs%ccof(m,i,lop,natom,1)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !---> perform brillouin zone integration and sum over bands

    DO ntyp = 1,atoms%ntype
       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          DO i = 1, noccbd
             qal21(l,ntyp,i)= qal21(l,ntyp,i)  + &
                  qaclo(i,lo,ntyp)*denCoeffsOffdiag%uulo21n(lo,ntyp) +&
                  qcloa(i,lo,ntyp)*denCoeffsOffdiag%ulou21n(lo,ntyp) +&
                  qclob(i,lo,ntyp)*denCoeffsOffdiag%ulod21n(lo,ntyp) +&
                  qbclo(i,lo,ntyp)*denCoeffsOffdiag%dulo21n(lo,ntyp)
          END DO
          DO lop = 1,atoms%nlo(ntyp)
             IF (atoms%llo(lop,ntyp).EQ.l) THEN
                DO i = 1, noccbd
                   qal21(l,ntyp,i)= qal21(l,ntyp,i)  + &
                        qlo(i,lop,lo,ntyp)*denCoeffsOffdiag%uloulop21n(lop,lo,ntyp)
                ENDDO
             ENDIF
          ENDDO
       END DO
    END DO

    DO n = 1,atoms%ntype
       fac = 1./atoms%neq(n)
       qal21(:,n,:) = qal21(:,n,:) * fac
    ENDDO
    !
    ! rotate into global frame
    !
    TYPE_loop : DO n = 1,atoms%ntype
       chi(1,1) =  EXP(-ImagUnit*nococonv%alph(n)/2)*COS(nococonv%beta(n)/2)
       chi(1,2) = -EXP(-ImagUnit*nococonv%alph(n)/2)*SIN(nococonv%beta(n)/2)
       chi(2,1) =  EXP( ImagUnit*nococonv%alph(n)/2)*SIN(nococonv%beta(n)/2)
       chi(2,2) =  EXP( ImagUnit*nococonv%alph(n)/2)*COS(nococonv%beta(n)/2)
       state : DO i = 1, noccbd
          lls : DO l = 0,3
             CALL rot_den_mat(nococonv%alph(n),nococonv%beta(n),&
                  dos%qal(l,n,ev_list(i),ikpt,1),dos%qal(l,n,ev_list(i),ikpt,2),qal21(l,n,i))
             IF (.FALSE.) THEN
                IF (n==1) WRITE(*,'(3i3,4f10.5)') l,n,i,qal21(l,n,i),dos%qal(l,n,ev_list(i),ikpt,:)
                q_loc(1,1) = dos%qal(l,n,ev_list(i),ikpt,1); q_loc(2,2) = dos%qal(l,n,ev_list(i),ikpt,2)
                q_loc(1,2) = qal21(l,n,i); q_loc(2,1) = CONJG(q_loc(1,2))
                q_hlp = MATMUL( TRANSPOSE( CONJG(chi) ) ,q_loc)
                q_loc = MATMUL(q_hlp,chi)
                qmat(l,n,i,1) = REAL(q_loc(1,1))
                qmat(l,n,i,2) = REAL(q_loc(1,2))
                qmat(l,n,i,3) = AIMAG(q_loc(1,2))
                qmat(l,n,i,4) = REAL(q_loc(2,2))
                IF (n==1) WRITE(*,'(3i3,4f10.5)') l,n,i,qmat(l,n,i,:)
             ENDIF
          ENDDO lls
       ENDDO state
    ENDDO TYPE_loop

  END SUBROUTINE qal_21
END MODULE m_qal21
