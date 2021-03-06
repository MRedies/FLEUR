!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_setlomap
      CONTAINS
        SUBROUTINE setlomap(ntyp,l_useapw,atoms)
          !***********************************************************************
          ! sets up nlol and lo1l
          ! 
          ! nlo     : number of local orbitals for each atom type
          ! llo     : the l quantum numbers of the local orbitals
          ! nlol    : the of local orbitals with a certain l (for each atom type)
          ! lo1l    : the number of the first local orbital with that l
          ! l_dulo  : if .true., this is a local orbital formed with a $\dot u$
          !
          ! p.kurz jul. 1996
          !***********************************************************************

          USE m_juDFT
          USE m_types_atoms
          USE m_constants

          IMPLICIT NONE
          !     ..
          !     .. Scalar Arguments ..
          INTEGER, INTENT (IN) :: ntyp
          LOGICAL, INTENT (IN) :: l_useapw
          TYPE(t_atoms),INTENT(INOUT)::atoms
          !     ...
          INTEGER ilo,l
          !     ..
          WRITE (oUnit,FMT=8000) atoms%nlo(ntyp), (atoms%llo(ilo,ntyp),ilo=1,atoms%nlo(ntyp))
8000      FORMAT ('the number of local orbitals for this atom type is: ',i3,&
               &       /,'the corresponding values of l are: ',30i4)
          IF (atoms%nlo(ntyp)>atoms%nlod)  CALL juDFT_error("nlo > nlod!!!",calledby&
               &     ="setlomap")

          DO l = 0,atoms%llod
             atoms%nlol(l,ntyp) = 0
             atoms%lo1l(l,ntyp) = 0
          END DO
          DO ilo = 1,atoms%nlod
             atoms%l_dulo(ilo,ntyp) = .FALSE.
             atoms%ulo_der(ilo,ntyp)= 0
          ENDDO
          l = -1

          DO ilo = 1,atoms%nlo(ntyp)
             !+gu
             IF (atoms%llo(ilo,ntyp).LT.0) THEN
                IF (l_useapw) THEN
                   atoms%l_dulo(ilo,ntyp) = .TRUE.
                ELSE
                   atoms%l_dulo(ilo,ntyp) = .FALSE.
                   atoms%ulo_der(ilo,ntyp)= -atoms%llo(ilo,ntyp)/10+1
                   atoms%llo(ilo,ntyp) = atoms%llo(ilo,ntyp)+10*(atoms%ulo_der(ilo,ntyp)-1)
                ENDIF
                atoms%llo(ilo,ntyp) = ABS( atoms%llo(ilo,ntyp) ) - 1
                WRITE(oUnit,'(A,I2,A,I2)') 'I use',atoms%ulo_der(ilo,ntyp),&
                     &       '. derivative of l =',atoms%llo(ilo,ntyp)
             ELSE
                atoms%l_dulo(ilo,ntyp) = .FALSE.
             ENDIF
             !-gu
             IF (atoms%llo(ilo,ntyp)>atoms%llod)  CALL juDFT_error(" l > llod!!!",&
                  &        calledby="setlomap")
             IF (atoms%llo(ilo,ntyp).LT.l) THEN
                WRITE (oUnit,FMT=*)&
                     &        'setlomap: please specify the l quantum numbers ',&
                     &        'of the local orbitals is ascending order.'
                CALL juDFT_error("LO-setup",calledby="setlomap")
             END IF
             IF (atoms%llo(ilo,ntyp).GT.l) THEN
                l = atoms%llo(ilo,ntyp)
                atoms%lo1l(l,ntyp) = ilo
             END IF
             atoms%nlol(l,ntyp) = atoms%nlol(l,ntyp) + 1
          END DO

          RETURN
        END SUBROUTINE setlomap
      END
