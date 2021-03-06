!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_relaxation
  USE m_judft
#ifdef CPP_MPI 
   use mpi 
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC relaxation !This is the interface. Below there are internal subroutines for bfgs, simple mixing, CG ...

CONTAINS
  SUBROUTINE relaxation(fmpi,input,atoms,cell,sym,oneD,vacuum,force_new,energies_new)
    !This routine uses the current force,energies and atomic positions to
    !generate a displacement in a relaxation step.
    !The history is taken into account by read_relax from m_relaxio
    !After generating new positions the code stops
    USE m_types
    USE m_constants
    USE m_relaxio
    USE m_mixing_history
    USE m_chkmt
    USE m_types_xml

    TYPE(t_mpi),INTENT(IN)   :: fmpi
    TYPE(t_input),INTENT(IN) :: input
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_oneD),INTENT(IN)  :: oneD
    TYPE(t_vacuum),INTENT(IN):: vacuum
    TYPE(t_cell),INTENT(IN)  :: cell
    REAL,INTENT(in)          :: force_new(:,:),energies_new !data for this iteration

    REAL,ALLOCATABLE :: pos(:,:,:),force(:,:,:),energies(:)
    REAL,ALLOCATABLE :: displace(:,:),old_displace(:,:), tempDisplace(:,:)
    REAL,ALLOCATABLE :: totalDisplace(:,:)
    REAL             :: dispAll(3,atoms%nat), overlap(0:atoms%ntype,atoms%ntype)
    REAL             :: dispLength, maxDisp, limitDisp
    INTEGER          :: iType,ierr, numDispReduce
    LOGICAL          :: l_conv

    !to calculate the current displacement
    Type(t_xml)   :: xml
    TYPE(t_atoms) :: atoms_non_displaced
    TYPE(t_atoms) :: tempAtoms

    IF (fmpi%irank==0) THEN
       call xml%init()
       ALLOCATE(pos(3,atoms%ntype,1));
       DO iType = 1, atoms%ntype
          pos(:,iType,1)=atoms%pos(:,SUM(atoms%neq(:iType-1))+1)
       END DO
       ALLOCATE(force(3,atoms%ntype,1)); force(:,:,1)=force_new
       ALLOCATE(energies(1));energies(1)=energies_new
       ALLOCATE(displace(3,atoms%ntype),old_displace(3,atoms%ntype))
       ALLOCATE(tempDisplace(3,atoms%ntype),totalDisplace(3,atoms%ntype))

       !Remove force components that are not selected for relaxation
       DO iType = 1, atoms%ntype
          IF (atoms%l_geo(iType)) THEN
             force(:,iType,1)=force(:,iType,1)*REAL(atoms%relax(:,iType))
          ELSE
             force(:,iType,1)=0.0
          ENDIF
       ENDDO

       ! add history
       CALL read_relax(pos,force,energies)

       !determine new positions
       IF (SIZE(energies)==1.OR.input%forcemix==0) THEN
          !no history present simple step
          ! choose a reasonable first guess for scaling
          ! this choice is based on a Debye temperature of 330K;
          ! modify as needed
          !alpha = (250.0/(MAXVAL(atoms%zatom)*input%xa))*((330./input%thetad)**2)
          CALL simple_step(input%forcealpha,0.25,force,displace)
       ELSE IF (input%forcemix==1) THEN
          CALL simple_cg(pos,force,displace)
       ELSE IF (input%forcemix==2) THEN
          CALL simple_bfgs(pos,force,displace)
       ELSE
          CALL juDFT_error('Unknown mixing scheme for forces', calledby='relaxation')
       END IF

       !Check for convergence of forces/displacements
       maxDisp = 0.0
       l_conv = .TRUE.
       DO iType = 1, atoms%ntype
          IF (DOT_PRODUCT(force(:,iType,SIZE(force,3)),force(:,iType,SIZE(force,3)))>input%epsforce**2) l_conv=.FALSE.
          dispLength = SQRT(DOT_PRODUCT(displace(:,iType),displace(:,iType)))
          maxDisp = MAX(maxDisp,dispLength)
          IF (dispLength>input%epsdisp) l_conv=.FALSE.
       ENDDO

       ! Limit the maximal displacement in a single force relaxation step to limitDisp = 0.2 a_0.
       limitDisp = 0.2
       IF(maxDisp.GT.limitDisp) THEN
          displace(:,:) = limitDisp*displace(:,:) / maxDisp
       END IF

       !New displacements relative to positions in inp.xml
       !CALL read_displacements(atoms,old_displace)
       call atoms_non_displaced%read_xml(xml)
       call xml%freeResources()
       call atoms_non_displaced%init(cell)
       DO iType = 1, atoms%ntype
          old_displace(:,iType) = atoms%pos(:,SUM(atoms%neq(:iType-1))+1) - &
                                  atoms_non_displaced%pos(:,SUM(atoms%neq(:iType-1))+1)
       END DO

       numDispReduce = 0
       overlap=1.0
       DO WHILE(ANY(overlap>1E-10))
          overlap = 0.0
          totalDisplace=displace+old_displace

          tempAtoms = atoms_non_displaced
          tempDisplace = MATMUL(cell%bmat,totalDisplace)/tpi_const

          CALL rotate_to_all_sites(tempDisplace,atoms_non_displaced,cell,sym,dispAll)
          tempAtoms%taual(:,:)=atoms_non_displaced%taual(:,:)+dispAll(:,:)
          tempAtoms%pos=MATMUL(cell%amat,tempAtoms%taual)
          CALL chkmt(tempAtoms,input,vacuum,cell,oneD,.TRUE.,overlap=overlap)

          IF (ANY(overlap>1E-10)) THEN
             numDispReduce = numDispReduce + 1
             IF (numDispReduce.GE.3) THEN
                CALL juDFT_warn("Strong MT spheres crash in structural relaxation")
             END IF
             displace(:,:) = 0.5 * displace(:,:)
             WRITE(oUnit,*) 'Automatically reducing atom displacements because MT spheres crash into each other!'
             WRITE(*,*) 'Automatically reducing atom displacements because MT spheres crash into each other!'
          END IF
       END DO

       !Write file
       CALL write_relax(pos,force,energies,totalDisplace)


    ENDIF
#ifdef CPP_MPI
    CALL MPI_BCAST(l_conv,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
#endif
    IF (l_conv) THEN
       CALL judft_end("Structural relaxation: Done",fmpi%irank)
    ELSE
       CALL mixing_history_reset(fmpi)
       CALL judft_end("Structural relaxation: new displacements generated",fmpi%irank)
    END IF
  END SUBROUTINE relaxation



  SUBROUTINE simple_step(alpha,maxdisp,force,displace)
    !-----------------------------------------------
    IMPLICIT NONE
    REAL,INTENT(in)  :: alpha,maxdisp
    REAL,INTENT(in)  :: force(:,:,:)
    REAL,INTENT(OUT) :: displace(:,:)

    real :: corr

    displace = alpha*force(:,:,SIZE(force,3))
    corr=maxdisp/maxval(abs(displace))
    if (corr<1.0) displace = corr*alpha*force(:,:,size(force,3))

  END SUBROUTINE simple_step

  SUBROUTINE simple_bfgs(pos,force,shift)
    !-----------------------------------------------
    !  Simple BFGS method to calculate shift out of old positions and forces
    !-----------------------------------------------
    USE m_constants
    IMPLICIT NONE
    REAL,INTENT(in)  :: pos(:,:,:),force(:,:,:)
    REAL,INTENT(OUT) :: shift(:,:)

    INTEGER         :: n,i,j,hist_length,n_force
    REAL,ALLOCATABLE:: h(:,:)
    REAL,ALLOCATABLE:: p(:),y(:),v(:)
    REAL            :: py,yy,gamma

    n_force=3*SIZE(pos,2)
    ALLOCATE(h(n_force,n_force))
    ALLOCATE(p(n_force),y(n_force),v(n_force))

    !calculate approx. Hessian
    !initialize H
    h = 0.0
    DO n = 1,n_force
       h(n,n) = 1.0
    ENDDO
    !loop over all iterations (including current)
    hist_length=SIZE(pos,3)
    DO n = 2,hist_length
       ! differences
       p(:) = RESHAPE(pos(:,:,n)-pos(:,:,n-1),(/SIZE(p)/))
       y(:) = RESHAPE(force(:,:,n-1)-force(:,:,n),(/SIZE(p)/))
       ! get necessary inner products and H|y>
       py = DOT_PRODUCT(p,y)
       v = MATMUL(y,h)
       yy = DOT_PRODUCT(y,v)
       !check that update will leave h positive definite;
       IF (py <= 0.0) THEN
          WRITE (oUnit,*) '  bfgs: <p|y> < 0'
          WRITE (oUnit,*) '  check convergence of forces'
          !Starting over with initial hessian
          h = 0.0
          DO j = 1,n_force
             h(j,j) = 1.0
          ENDDO
          CYCLE
       ELSE
          !update h
          IF (n == 2) THEN
             gamma = py/yy
          ELSE
             gamma = 1.0
          ENDIF
          DO j = 1,n_force
             DO i = 1,n_force
                h(i,j) = (h(i,j) - (v(i)*p(j)+p(i)*v(j))/py)*gamma + (1.+gamma*yy/py)*p(i)*p(j)/py
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    y(:) = RESHAPE(force(:,:,hist_length),(/SIZE(p)/))
    shift = RESHAPE(MATMUL(y,h),SHAPE(shift))
  END SUBROUTINE simple_bfgs

  SUBROUTINE simple_cg(pos,force,shift)
    !-----------------------------------------------
    IMPLICIT NONE
    REAL,INTENT(in)  :: pos(:,:,:),force(:,:,:)
    REAL,INTENT(OUT) :: shift(:,:)

    REAL                :: f1(3,SIZE(pos,2)),f2(3,SIZE(pos,2))
    REAL                :: dist(3,SIZE(pos,2))
    REAL                :: eps
    INTEGER             :: n_old, i, j

    eps = 1.0e-9

    n_old = SIZE(pos,3)-1

    dist(:,:) = pos(:,:,n_old+1)-pos(:,:,n_old)
    DO i = 1, SIZE(pos,2)
       DO j = 1, 3
          IF(ABS(dist(j,i)).LT.eps) dist(j,i) = eps ! To avoid calculation of 0.0/0.0 below.
       END DO
    END DO

    f1 = (force(:,:,n_old+1)-force(:,:,n_old))/dist
    f2 = force(:,:,n_old+1)-f1*pos(:,:,n_old+1)
    shift = -1.*f2/f1-force(:,:,n_old+1)
  END SUBROUTINE simple_cg
END MODULE m_relaxation
