MODULE m_hssrwu
  USE m_juDFT
  !
  !*********************************************************************
  !     updates the hamiltonian and overlap matrices with the
  !     contributions from the spheres, both spherical and non-
  !     spherical, for step forward approach
  !                r. wu  1992
  !*********************************************************************
CONTAINS
  SUBROUTINE hssr_wu(atoms,sym, jsp,el,ne,usdus,lapw,input,&
       tlmplm,acof,bcof,ccof, h_r,s_r,h_c,s_c)
    !
    USE m_types
    IMPLICIT NONE


    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_usdus),INTENT(IN)       :: usdus
    TYPE(t_tlmplm),INTENT(IN)      :: tlmplm
    TYPE(t_lapw),INTENT(IN)        :: lapw
    TYPE(t_input),INTENT(IN)       :: input
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp,ne
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
    COMPLEX, INTENT (IN) :: acof(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
    COMPLEX, INTENT (IN) :: bcof(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
    COMPLEX, INTENT (IN) :: ccof(-atoms%llod:atoms%llod,input%neig,atoms%nlod,atoms%nat)

    REAL,    OPTIONAL,INTENT (INOUT) :: h_r(input%neig,input%neig),s_r(input%neig,input%neig)
    COMPLEX, OPTIONAL,INTENT (INOUT) :: h_c(input%neig,input%neig),s_c(input%neig,input%neig)

    !     ..
    !     .. Local Scalars ..
    COMPLEX dtd,dtu,hij,sij,utd,utu
    REAL invsfct
    INTEGER i,im,in,j,k,ke,l,l1,ll1,lm,lmp,lwn ,m1,n,na,nn,lmplm,m
    LOGICAL :: l_real
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: a(:,:),b(:,:),ax(:),bx(:)
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC cmplx,conjg,exp,REAL,sqrt
    !     ..

    l_real=PRESENT(h_r)

    ALLOCATE ( a(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2)),ax(input%neig) )
    ALLOCATE ( b(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2)),bx(input%neig) )
    na = 0
    DO n = 1,atoms%ntype        ! loop over atom-types
       lwn = atoms%lmax(n)
       DO nn = 1,atoms%neq(n)    ! loop over atoms
          na = na + 1
          !+inv
          IF ((sym%invsat(na).EQ.0) .OR. (sym%invsat(na).EQ.1)) THEN
             CALL timestart("hssr_wu: spherical")
             IF (sym%invsat(na).EQ.0) invsfct = 1.0
             IF (sym%invsat(na).EQ.1) invsfct = SQRT(2.0)
             DO lm = 0, atoms%lmaxd*(atoms%lmaxd+2)
                DO ke = 1, ne
                   a(ke,lm) = invsfct*acof(ke,lm,na)
                   b(ke,lm) = invsfct*bcof(ke,lm,na)
                ENDDO
             ENDDO

             DO l = 0,lwn                    ! l loop
                DO m = -l,l                  ! m loop
                   lmp = l* (l+1) + m
                   DO i = 1,ne               ! matrix update
                      DO j = 1,i - 1
                         sij = a(i,lmp)*CONJG(a(j,lmp)) +&
                              b(i,lmp)*CONJG(b(j,lmp))*usdus%ddn(l,n,jsp)
                         hij = sij * el(l,n,jsp) +&
                              0.5 * ( a(i,lmp)*CONJG(b(j,lmp)) +&
                              b(i,lmp)*CONJG(a(j,lmp)) )
                         IF (l_real) THEN
                            s_r(i,j) = s_r(i,j) + REAL(sij)
                            h_r(i,j) = h_r(i,j) + REAL(hij)
                         ELSE
                            s_c(i,j) = s_c(i,j) + sij
                            h_c(i,j) = h_c(i,j) + hij
                         ENDIF
                      ENDDO
                   ENDDO
                   DO i = 1,ne
                      sij = a(i,lmp)*CONJG(a(i,lmp)) +&
                           b(i,lmp)*CONJG(b(i,lmp))*usdus%ddn(l,n,jsp)
                      hij = sij * el(l,n,jsp) +&
                           0.5 * ( a(i,lmp)*CONJG(b(i,lmp)) +&
                           b(i,lmp)*CONJG(a(i,lmp)) )
                      IF (l_real) THEN
                         s_r(i,i) = s_r(i,i) + REAL(sij)
                         h_r(i,i) = h_r(i,i) + REAL(hij)
                      ELSE
                         s_c(i,i) = s_c(i,i) + sij
                         h_c(i,i) = h_c(i,i) + hij
                      ENDIF
                   ENDDO
                ENDDO        ! m
             ENDDO           ! l

             CALL timestop("hssr_wu: spherical")
             CALL timestart("hssr_wu: non-spherical")
             IF (atoms%lnonsph(n) >= 0 ) THEN
                DO l = 0,atoms%lnonsph(n)
                   DO m = -l,l

                      lmp = l* (l+1) + m
                      ax(:) = CMPLX(0.0,0.0)
                      bx(:) = CMPLX(0.0,0.0)

                      DO l1 = 0,atoms%lnonsph(n)         ! l', m' loop
                         DO m1 = -l1,l1
                            lm = l1* (l1+1) + m1
                            utu = CONJG(tlmplm%h_loc(lmp,lm,n,jsp,jsp))
                            dtd = CONJG(tlmplm%h_loc(lmp+tlmplm%h_loc2(n),lm+tlmplm%h_loc2(n),n,jsp,jsp))
                            utd = CONJG(tlmplm%h_loc(lmp+tlmplm%h_loc2(n),lm,n,jsp,jsp))
                            dtu = CONJG(tlmplm%h_loc(lmp,lm+tlmplm%h_loc2(n),n,jsp,jsp))
                            !--->    update ax, bx
                            DO k = 1,ne
                               ax(k) = ax(k) + utu*CONJG(a(k,lm))+&
                                    utd*CONJG(b(k,lm))
                               bx(k) = bx(k) + dtu*CONJG(a(k,lm))+&
                                    dtd*CONJG(b(k,lm))
                            ENDDO

                         ENDDO    ! m1
                      ENDDO       ! l1
                      !
                      !
                      !--->    update hamiltonian
                      IF (l_real) THEN
                         DO i = 1,ne
                            DO j = 1,i - 1
                               hij = a(i,lmp)*ax(j) + b(i,lmp)*bx(j)
                               h_r(i,j) = h_r(i,j) + REAL(hij)
                            ENDDO
                         ENDDO
                      ELSE
                         DO i = 1,ne
                            DO j = 1,i - 1
                               hij = a(i,lmp)*ax(j) + b(i,lmp)*bx(j)
                               h_c(i,j) = h_c(i,j) + hij
                            ENDDO
                         ENDDO
                      ENDIF

                      IF (l_real) THEN
                         DO i = 1,ne
                            h_r(i,i) = h_r(i,i) + REAL(a(i,lmp)*ax(i)+ b(i,lmp)*bx(i))
                         ENDDO
                      ELSE
                         DO i = 1,ne
                            h_c(i,i) = h_c(i,i) + a(i,lmp)*ax(i)+ b(i,lmp)*bx(i)
                         ENDDO
                      ENDIF

                   ENDDO ! m
                ENDDO   ! l
             ENDIF     ! atoms%lnonsph >=0
             CALL timestop("hssr_wu: non-spherical")
             !-inv
          ENDIF ! invsatom = 0 or 1
       ENDDO   ! loop over atoms
    ENDDO     ! loop over atom-types

    DEALLOCATE ( a, b, ax, bx )
  END SUBROUTINE hssr_wu
END MODULE m_hssrwu
