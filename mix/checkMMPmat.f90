MODULE m_checkMMPmat

   !Check whether the given density matrix makes sense (only diagonal)

   USE m_types
   USE m_juDFT
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE checkMMPmat(indStart,indEnd,atoms,input,mmpmat)

      INTEGER,             INTENT(IN)  :: indStart, indEnd
      TYPE(t_atoms),       INTENT(IN)  :: atoms
      TYPE(t_input),       INTENT(IN)  :: input
      COMPLEX,             INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:,:)

      !which elements are considered to cause an error
      REAL, PARAMETER :: lowBound = -0.01
      REAL, PARAMETER :: highBound = 1.05 !Attention keep in mind jspins=1

      LOGICAL l_err
      INTEGER i_u,l,ispin,m
      REAL spindeg

      spindeg = 3-input%jspins
      l_err = .FALSE.
      DO i_u = indStart, indEnd
         l = atoms%lda_u(i_u)%l
         !Check the diagonal elements
         DO ispin = 1, input%jspins
            DO m = -l,l
               IF(REAL(mmpmat(m,m,i_u,ispin)).LT.lowBound.OR. &
                  REAL(mmpmat(m,m,i_u,ispin)).GT.highBound*spindeg) THEN
                  l_err = .TRUE.
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(l_err) THEN
         WRITE(*,*) "-----------------------------------------------------------------"
         WRITE(*,*) "Using the Quasi-Newton methods for mixing and LDA+U"
         WRITE(*,*) "from the beginning of the SCF calculaion can be unstable."
         WRITE(*,*) "You can reset the mixing_history, use straight mixing for "
         WRITE(*,*) "the first iterations or use linear mixing for the density matrix"
         WRITE(*,*) "-----------------------------------------------------------------"
         CALL juDFT_error("Invalid elements in mmpmat", calledby="checkMMPmat")

      ENDIF

   END SUBROUTINE checkMMPmat

END MODULE m_checkMMPmat