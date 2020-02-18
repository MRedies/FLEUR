!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and avhttps://gcc.gnu.org/onlinedocs/gfortran/SQRT.htmlailable as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!------------------------------------------------------------------------------
!  This routine allows to rotate the cdn in a way that the direction of magnetization aligns with the direction of the spin quantization axis.
!  This routine also allows to reverse the rotation by using the angles stored in atoms (phi_mt_avg,theta_mt_avg) which are generated by the
!  routine magnMomFromDen.
!
! Robin Hilgers, Nov '19
MODULE m_alignSpinAxisMagn


USE m_magnMomFromDen
USE m_types
USE m_types_fleurinput
USE m_flipcdn
USE m_constants
USE m_polangle
IMPLICIT NONE

CONTAINS
SUBROUTINE rotateMagnetToSpinAxis(vacuum,sphhar,stars&
,sym,oneD,cell,noco,nococonv,input,atoms,den,l_firstIt)
   TYPE(t_input), INTENT(IN)     :: input
   TYPE(t_atoms), INTENT(IN)     :: atoms
   TYPE(t_noco), INTENT(IN)      :: noco
   TYPE(t_nococonv),INTENT(INOUT):: nococonv
   TYPE(t_stars),INTENT(IN)      :: stars
   TYPE(t_vacuum),INTENT(IN)     :: vacuum
   TYPE(t_sphhar),INTENT(IN)     :: sphhar
   TYPE(t_sym),INTENT(IN)        :: sym
   TYPE(t_oneD),INTENT(IN)       :: oneD
   TYPE(t_cell),INTENT(IN)       :: cell
   TYPE(t_potden), INTENT(INOUT) :: den

   LOGICAL                       :: l_firstIt
   REAL                          :: moments(3,atoms%ntype)
   REAL                          :: phiTemp(atoms%ntype),thetaTemp(atoms%ntype)
   REAL                          :: diffT(atoms%ntype),diffP(atoms%ntype),eps, zeros(atoms%ntype)

   INTEGER                       :: i
   eps=0.0001
   zeros(:)=0.0
   IF(l_firstIt) THEN
     CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,nococonv%alph,nococonv%beta,den)
     nococonv%alph=0.0
     nococonv%beta=0.0
   END IF

   CALL magnMomFromDen(input,atoms,noco,den,moments,thetaTemp,phiTemp)
   diffT=thetaTemp-nococonv%beta
   diffP=-(phiTemp-nococonv%alph)
   DO i=1, atoms%ntype
     IF (abs(diffT(i)).LE.eps) diffT(i)=0.0
     IF (abs(diffP(i)).LE.eps) diffP(i)=0.0
   END DO
   CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,-diffP,zeros,den)
   CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,zeros,-diffT,den)
   nococonv%beta=nococonv%beta+diffT
   nococonv%alph=nococonv%alph+diffP
   DO i=1, atoms%ntype
     IF (abs(nococonv%beta(i)).LE.eps) nococonv%beta(i)=0.0
     IF (abs(nococonv%alph(i)).LE.eps) nococonv%alph(i)=0.0
   END DO


   write(*,*) "Noco Phi"
   write(*,*) nococonv%alph
   write(*,*) "Noco Theta"
   write(*,*) nococonv%beta
END SUBROUTINE rotateMagnetToSpinAxis


SUBROUTINE rotateMagnetFromSpinAxis(noco,nococonv,vacuum,sphhar,stars&
,sym,oneD,cell,input,atoms,inDen, den)
   TYPE(t_input), INTENT(IN)  :: input
   TYPE(t_atoms), INTENT(IN)  :: atoms
   TYPE(t_noco), INTENT(IN)	  :: noco
   TYPE(t_nococonv), INTENT(INOUT)	 :: nococonv
   TYPE(t_stars),INTENT(IN)	  :: stars
   TYPE(t_vacuum),INTENT(IN)     :: vacuum
   TYPE(t_sphhar),INTENT(IN)     :: sphhar
   TYPE(t_sym),INTENT(IN)        :: sym
   TYPE(t_oneD),INTENT(IN)	 :: oneD
   TYPE(t_cell),INTENT(IN)	 :: cell
   TYPE(t_potden), INTENT(INOUT) ::  inDen
   TYPE(t_potden), OPTIONAL,INTENT(INOUT) :: den

   INTEGER                            :: i
   REAL                          :: phiTemp(atoms%ntype),thetaTemp(atoms%ntype), zeros(atoms%ntype)
   REAL                          :: moments(3,atoms%ntype)
   zeros(:)=0.0

   CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,nococonv%alph,zeros,inDen)
   CAlL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,zeros,nococonv%beta,inDen)
   IF (present(den)) THEN
     CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,nococonv%alph,zeros,den)
     CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,zeros,nococonv%beta,den)
   END IF
   nococonv%alph=0
   nococonv%beta=0


END SUBROUTINE rotateMagnetFromSpinAxis


END MODULE m_alignSpinAxisMagn
