MODULE m_tetrahedronInit

   CONTAINS

   SUBROUTINE tetrahedronInit(ikpt,kpts,input,neig,eig,gCoeffs,ef,ez,nz,weights,boundInd)

      !This Subroutine sets up the weights for the brillouin zone integration using
      !the tetrahedron method for the greens function calculations
      !
      !If we are in a non-magnetic or collinear case we can decompose the greens function
      !into real and imaginary part and obtain the real part from the imaginary part 
      !via Kramers-Kronig transformation. Then we only need the weights for a dos calculations
      !
      !Otherwise we need the weights for the resolvent function

      USE m_types
      USE m_juDFT
      USE m_dosWeights
      USE m_resWeights

      IMPLICIT NONE

      INTEGER,                INTENT(IN)  :: ikpt !Current k-point
      TYPE(t_kpts),           INTENT(IN)  :: kpts 
      TYPE(t_input),          INTENT(IN)  :: input
      INTEGER,                INTENT(IN)  :: neig, nz
      REAL,                   INTENT(IN)  :: eig(neig,kpts%nkpt)
      TYPE(t_greensfCoeffs),  INTENT(IN)  :: gCoeffs
      REAL,                   INTENT(IN)  :: ef
      COMPLEX,                INTENT(IN)  :: ez(nz)

      COMPLEX,                INTENT(INOUT) :: weights(gCoeffs%ne,neig)
      INTEGER,                INTENT(INOUT) :: boundInd(neig,2) !Indices in between which the weights are non zero 
                                                                !to save computation time

      weights = 0.0
      boundInd = 0

      IF(input%film) THEN
         
      ELSE
         IF(input%l_resolvent) THEN
            CALL resWeights(ikpt,kpts,neig,eig,ez,nz,weights,boundInd)
         ELSE
            CALL dosWeights(ikpt,kpts,neig,eig,gCoeffs,ef,weights,boundInd)
         ENDIF
      ENDIF

   END SUBROUTINE tetrahedronInit

END MODULE m_tetrahedronInit