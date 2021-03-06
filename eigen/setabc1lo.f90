!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_setabc1lo
!*********************************************************************
! calculates the (lower case) a, b and c coefficients for the local
! orbitals. The radial function of the local orbital is a linear
! combination of the apw radial function and its derivative and the
! extra radial funtion (a*u + b*udot + c*ulo). This function is zero
! and has zero derivative at the muffin tin boundary.
! Philipp Kurz 99/04
!*********************************************************************
      CONTAINS
      SUBROUTINE setabc1lo(atoms, ntyp,ud,usp, alo1,blo1,clo1)
      USE m_types
      IMPLICIT NONE

      TYPE(t_atoms),INTENT(IN)   :: atoms
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN)  :: ntyp,usp
!     ..
!     .. Array Arguments ..
      TYPE(t_usdus),INTENT(IN):: ud
      REAL,    INTENT (OUT) :: alo1(:,:),blo1(:,:),clo1(:,:)
!     ..
!     .. Local Scalars ..
      REAL ka,kb,ws
      INTEGER l,lo
      LOGICAL apw_at
!     ..
!     ..
! look, whether 'ntyp' is a APW atom; then set apw_at=.true.
      !
      apw_at=ANY(atoms%l_dulo(:atoms%nlo(ntyp),ntyp))

      DO lo = 1,atoms%nlo(ntyp)
         l = atoms%llo(lo,ntyp)
         IF (apw_at) THEN
           IF (atoms%l_dulo(lo,ntyp)) THEN
! udot lo
             ka=sqrt(1+(ud%us(l,ntyp,usp)/ud%uds(l,ntyp,usp))**2* ud%ddn(l,ntyp,usp))
             alo1(lo,usp)=1.00 / ka
             blo1(lo,usp)=-ud%us(l,ntyp,usp)/ (ud%uds(l,ntyp,usp) * ka )
             clo1(lo,usp)=0.00
           ELSE
! u2 lo
             alo1(lo,usp)=1.00
             blo1(lo,usp)=0.00
             clo1(lo,usp)=-ud%us(l,ntyp,usp)/ud%ulos(lo,ntyp,usp)
           ENDIF
         ELSE
           ws = ud%uds(l,ntyp,usp)*ud%dus(l,ntyp,usp) - ud%us(l,ntyp,usp)*ud%duds(l,ntyp,usp)
           ka = 1.0/ws*(ud%duds(l,ntyp,usp)*ud%ulos(lo,ntyp,usp)- ud%uds(l,ntyp,usp)*ud%dulos(lo,ntyp,usp))
           kb = 1.0/ws* (ud%us(l,ntyp,usp)*ud%dulos(lo,ntyp,usp)- ud%dus(l,ntyp,usp)*ud%ulos(lo,ntyp,usp))
           clo1(lo,usp) = 1.0/sqrt(ka**2+ (kb**2)*ud%ddn(l,ntyp,usp)+1.0+ 2.0*ka*ud%uulon(lo,ntyp,usp)+&
                2.0*kb*ud%dulon(lo,ntyp,usp))
           alo1(lo,usp) = ka*clo1(lo,usp)
           blo1(lo,usp) = kb*clo1(lo,usp)
         ENDIF
      END DO

      END SUBROUTINE setabc1lo
      END MODULE m_setabc1lo
!
!         flo = alo1(lo,usp)*us(l,ntyp) + blo1(lo,usp)*uds(l,ntyp) +
!     +         clo1(lo,usp)*ulos(lo,ntyp)
!         dflo = alo1(lo,usp)*dus(l,ntyp) + blo1(lo,usp)*duds(l,ntyp) +
!     +          clo1(lo,usp)*dulos(lo,ntyp)
!         nflo = alo1(lo,usp)**2 + (blo1(lo,usp)**2)*ddn(l,ntyp) + clo1(lo,usp)**2 +
!     +          2*alo1(lo,usp)*clo1(lo,usp)*uulon(lo,ntyp) +
!     +          2*blo1(lo,usp)*clo1(lo,usp)*dulon(lo,ntyp)
