!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_banddos
  USE m_juDFT
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_banddos
  TYPE,EXTENDS(t_fleurinput_base):: t_banddos
     LOGICAL :: dos =.FALSE.
     LOGICAL :: band =.FALSE.
     LOGICAL :: l_mcd =.FALSE.
     LOGICAL :: l_orb =.FALSE.
     LOGICAL :: l_jDOS = .FALSE.
     LOGICAL :: vacdos =.FALSE.
     INTEGER :: ndir =0
     INTEGER :: orbCompAtom=0
     INTEGER :: jDOSAtom=0
     REAL    :: e1_dos=0.5
     REAL    :: e2_dos=-0.5
     REAL    :: sig_dos=0.015
     INTEGER :: projdos !selects one atomtype and prints the projected dos if there are to many atoms
     REAL    :: e_mcd_lo =-10.0
     REAL    :: e_mcd_up= 0.0
     LOGICAL :: unfoldband =.FALSE.
     INTEGER :: s_cell_x=1
     INTEGER :: s_cell_y=1
     INTEGER :: s_cell_z=1
     REAL    :: alpha,beta,gamma !For orbital decomp. (was orbcomprot)
   CONTAINS
     PROCEDURE :: read_xml=>read_xml_banddos
     PROCEDURE :: mpi_bc=>mpi_bc_banddos
  END TYPE t_banddos
CONTAINS
  SUBROUTINE mpi_bc_banddos(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_banddos),INTENT(INOUT)::this
    integer,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    if (present(irank)) THEN
       rank=irank
    else
       rank=0
    end if
    CALL mpi_bc(this%dos ,rank,mpi_comm)
    CALL mpi_bc(this%band ,rank,mpi_comm)
    CALL mpi_bc(this%l_mcd ,rank,mpi_comm)
    CALL mpi_bc(this%l_orb ,rank,mpi_comm)
    CALL mpi_bc(this%l_jDOS,rank,mpi_comm)
    CALL mpi_bc(this%vacdos ,rank,mpi_comm)
    CALL mpi_bc(this%ndir ,rank,mpi_comm)
    CALL mpi_bc(this%orbCompAtom,rank,mpi_comm)
    CALL mpi_bc(this%jDOSAtom,rank,mpi_comm)
    CALL mpi_bc(this%e1_dos,rank,mpi_comm)
    CALL mpi_bc(this%e2_dos,rank,mpi_comm)
    CALL mpi_bc(this%sig_dos,rank,mpi_comm)
    CALL mpi_bc(this%e_mcd_lo ,rank,mpi_comm)
    CALL mpi_bc(this%e_mcd_up,rank,mpi_comm)
    CALL mpi_bc(this%unfoldband ,rank,mpi_comm)
    CALL mpi_bc(this%s_cell_x,rank,mpi_comm)
    CALL mpi_bc(this%s_cell_y,rank,mpi_comm)
    CALL mpi_bc(this%s_cell_z,rank,mpi_comm)
    CALL mpi_bc(this%alpha,rank,mpi_comm)
    CALL mpi_bc(this%beta,rank,mpi_comm)
    CALL mpi_bc(this%gamma,rank,mpi_comm)

  END SUBROUTINE mpi_bc_banddos
  SUBROUTINE read_xml_banddos(this,xml)
    USE m_types_xml
    CLASS(t_banddos),INTENT(INOUT)::this
    TYPE(t_xml),INTENT(INOUT)::xml

    CHARACTER(len=300) :: xPathA, xPathB
    INTEGER::numberNodes,iType,i,na
    LOGICAL::l_orbcomp,l_jDOS
    this%band = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@band'))
    this%dos = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@dos'))
    this%vacdos = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@vacdos'))
    this%l_mcd = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@mcd'))

    numberNodes = xml%GetNumberOfNodes('/fleurInput/output/densityOfStates')

    IF ((this%dos).AND.(numberNodes.EQ.0)) THEN
       CALL juDFT_error("dos is true but densityOfStates parameters are not set!")
    END IF

    IF (numberNodes.EQ.1) THEN
       this%ndir = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/densityOfStates/@ndir'))
       this%e2_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/densityOfStates/@minEnergy'))
       this%e1_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/densityOfStates/@maxEnergy'))
       this%sig_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/densityOfStates/@sigma'))
    END IF
    IF (this%band) THEN
       this%dos=.TRUE.
       this%ndir = -4
       WRITE(*,*) 'band="T" --> Overriding "dos" and "ndir"!'
    ENDIF

    ! Read in optional magnetic circular dichroism parameters
    numberNodes = xml%GetNumberOfNodes('/fleurInput/output/magneticCircularDichroism')

    IF ((this%l_mcd).AND.(numberNodes.EQ.0)) THEN
       CALL juDFT_error("mcd is true but magneticCircularDichroism parameters are not set!", calledby = "read_xml_banddos")
    END IF

    IF (numberNodes.EQ.1) THEN
       this%e_mcd_lo = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/magneticCircularDichroism/@energyLo'))
       this%e_mcd_up = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/magneticCircularDichroism/@energyUp'))
    END IF

    na = 0
    DO iType = 1, xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
       WRITE(xPathA,*) '/fleurInput/atomGroups/atomGroup[',iType,']'
       DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/relPos')
          na = na + 1
          WRITE(xPathB,*) TRIM(ADJUSTL(xPathA))//'/relPos[',i,']'
          l_orbcomp = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@orbcomp'))
          IF(l_orbcomp) THEN
             IF(this%l_orb) THEN
                CALL juDFT_error("Multiple orbcomp flags set.", calledby = "read_xml_banddos")
             END IF
             this%l_orb = .TRUE.
             this%orbCompAtom = na
          ENDIF
          l_jDOS = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@jDOS'))
          IF(l_jDOS) THEN
             IF(this%l_jDOS) THEN
                CALL juDFT_error("Multiple jDOS flags set.", calledby="read_xml_banddos")
             ENDIF
             this%l_jDOS = .TRUE.
             this%jDOSAtom = na
          ENDIF
       ENDDO
    ENDDO
    IF(this%l_orb.AND.this%l_jDOS) THEN
       CALL juDFT_error("Both jDOS and orbcomp flag set", calledby="read_xml_banddos")
    ENDIF

    ! Read in optional parameter for unfolding bandstructure of supercell
    numberNodes = xml%GetNumberOfNodes('/fleurInput/output/unfoldingBand')

    IF (numberNodes.EQ.1) THEN
       this%unfoldband = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@unfoldBand'))
       this%s_cell_x = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellX'))
       this%s_cell_y = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellY'))
       this%s_cell_z = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellZ'))
    END IF
  END SUBROUTINE read_xml_banddos

END MODULE m_types_banddos
