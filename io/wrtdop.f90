      MODULE m_wrtdop
!     ****************************************************
!     write formatted density or potential onto unit 'nu'
!     e. wimmer   march 1985
!     ****************************************************
      CONTAINS
        SUBROUTINE wrtdop(&
             &                  stars,vacuum,atoms,sphhar,&
             &                  input,sym,&
             &                  nu,&
             &                  it,fr,fpw,fz,fzxy)
          !
          USE m_constants,ONLY: namat_const
          USE m_types
          IMPLICIT NONE
          !
          !     .. Scalar Arguments ..
          TYPE(t_stars),INTENT(IN)  :: stars
          TYPE(t_vacuum),INTENT(IN) :: vacuum
          TYPE(t_atoms),INTENT(IN)  :: atoms
          TYPE(t_sphhar),INTENT(IN) :: sphhar
          TYPE(t_input),INTENT(IN)  :: input
          TYPE(t_sym),INTENT(IN)    :: sym
          INTEGER, INTENT (IN) :: nu 
          INTEGER, INTENT (IN) :: it   
          !     ..
          !     .. Array Arguments ..
          COMPLEX, INTENT (IN):: fpw(stars%n3d,input%jspins),fzxy(vacuum%nmzxyd,stars%n2d-1,2,input%jspins)
          REAL,    INTENT (IN):: fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),fz(vacuum%nmzd,2,input%jspins)
          CHARACTER(len=8):: dop,iop,name(10)
          !     .. Local Scalars ..
          INTEGER i,ivac,izn,jsp,k,lh,n,na
          !     ..
          !     .. Intrinsic Functions ..
          INTRINSIC REAL
          !     ..
          !Defaults, as the character arrays are no longer used, we should
          !write some defaults for old fleur-versions
          name    ='        '
          name(10)='ordered*'
          dop     ='in/out  '
          iop     ='char/pot'
          WRITE (nu) name
          !          WRITE (6,FMT=8000) name
8000      FORMAT (' wrtdop title:',10a8)
          WRITE (nu) iop,dop,it
          DO  jsp = 1,input%jspins
             WRITE (nu) jsp
             WRITE (nu) atoms%ntype
             na = 1
             DO  n = 1,atoms%ntype
                izn = atoms%zatom(n) + 0.01
                WRITE (nu) namat_const(izn),n,atoms%jri(n),atoms%rmt(n),atoms%dx(n)
                WRITE (nu) atoms%ntypsy(na),sphhar%nlh(atoms%ntypsy(na))
                DO  lh = 0,sphhar%nlh(atoms%ntypsy(na))
                   WRITE (nu) lh
                   WRITE (nu) (fr(i,lh,n,jsp),i=1,atoms%jri(n))
                ENDDO
                na = na + atoms%neq (n)
             ENDDO
             WRITE (nu) stars%ng3
             IF (sym%invs) THEN
                WRITE (nu) (REAL(fpw(k,jsp)),k=1,stars%ng3)
             ELSE
                WRITE (nu) (fpw(k,jsp),k=1,stars%ng3)
             END IF
             IF (input%film) THEN
                DO  ivac = 1,vacuum%nvac
                   WRITE (nu) ivac
                   WRITE (nu) vacuum%nmz,vacuum%dvac,vacuum%delz
                   WRITE (nu) (fz(i,ivac,jsp),i=1,vacuum%nmz)
                   WRITE (nu) stars%ng2,vacuum%nmzxy
                   DO  k = 2,stars%ng2
                      IF (sym%invs2) THEN
                         WRITE (nu) (REAL(fzxy(i,k-1,ivac,jsp)),i=1,vacuum%nmzxy)
                      ELSE
                         WRITE (nu) (fzxy(i,k-1,ivac,jsp),i=1,vacuum%nmzxy)
                      END IF
                   ENDDO
                ENDDO
             END IF
          ENDDO
          !
          RETURN
        END SUBROUTINE wrtdop
      END MODULE m_wrtdop
      
