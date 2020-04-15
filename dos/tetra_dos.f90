!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_tetrados
   !----------------------------------------------------------------------
   !
   ! This subroutine evaluates the density of states (g) by the linear
   ! tetrahedron method on a energy grid (e) of 'ned' points.
   !
   ! ev()          ... eigenvalues
   ! qal()         ... partial charges
   ! ntetra)       ... number of tetrahedrons
   ! itetra(1-4,nt)... index of k-points forming tetrahedron nt
   ! voltet(nt)    ... volume of tetrahedron nt
   ! omega_bz      ... volume of irreducible part of BZ
   !
   !                                                      gb 2000
   !----------------------------------------------------------------------
   USE m_types

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE tetra_dos(qdim,neigd,ned,ntetra,nkpt,&
                        itetra,efermi,voltet,energy,nevk,&
                        ev,qal,g)

      INTEGER, INTENT(IN)    :: neigd,ned,qdim
      INTEGER, INTENT(IN)    :: ntetra,nkpt
      REAL,    INTENT(IN)    :: efermi
      INTEGER, INTENT(IN)    :: itetra(:,:) !(4,6*nkpt)
      INTEGER, INTENT(IN)    :: nevk(:)     !(nkpt)
      REAL,    INTENT(IN)    :: voltet(:)   !(6*nkpt)
      REAL,    INTENT(IN)    :: energy(:)   !(ned)
      REAL,    INTENT(IN)    :: qal(:,:,:)  !(lmax*ntype+3,neigd,nkpt)
      REAL,    INTENT(INOUT) :: ev(:,:)     !(neigd,nkpt)
      REAL,    INTENT(OUT)   :: g(:,:)      !(ned,lmax*ntype+3)

      INTEGER :: i,j,iBand,ikpt,ie,idim,itet
      REAL    :: ener,w
      REAL    :: weight(4),eval(4),ecmax(neigd),term(ned)
      REAL    :: wpar(qdim,neigd,nkpt)

      DO ikpt = 1,nkpt
         ev(nevk(ikpt)+1:neigd,ikpt) = 1.0e10
      ENDDO

      wpar = 0.0

      DO iBand = 1,neigd
         ecmax(iBand) = -1.0e25
         DO ikpt = 1,nkpt
            IF(ev(iBand,ikpt).GT.ecmax(iBand)) ecmax(iBand) = ev(iBand,ikpt)
         ENDDO
      ENDDO
      !
      !  check for energy degeneracies in tetrahedrons
      !
      DO itet = 1,ntetra
         DO iBand = 1,neigd
            DO i = 1,3
               DO j = i+1,4
                  IF (abs(ev(iBand,itetra(i,itet))-ev(iBand,itetra(j,itet))).LT.1.0e-7) THEN
                     ev(iBand,itetra(i,itet)) = ev(iBand,itetra(i,itet)) + i*(1.0e-7)*itet
                     ev(iBand,itetra(j,itet)) = ev(iBand,itetra(j,itet)) - i*(1.0e-7)*itet
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      WRITE (*,*) 'in tetra_dos'  ! do not remove  this statement

      !
      ! calculate partial weights
      !
      DO ikpt=1,nkpt
         DO iBand = 1,nevk(ikpt)
            DO itet = 1,ntetra
               IF (ALL(itetra(:,itet).ne.ikpt)) CYCLE

               eval(1:4) = ev(iBand,itetra(1:4,itet))

               IF(max(eval(1),eval(2),eval(3),eval(4)).GE.9999.9) CYCLE

               DO i=1,4
                  weight(i)=1.0
                  DO j=1,4
                     IF (i.NE.j) weight(i)=weight(i)*(eval(j)-eval(i))
                  ENDDO
                  weight(i)=6.0*voltet(itet)/weight(i)
                  DO idim=1,qdim
                     wpar(idim,iBand,itetra(i,itet)) =  wpar(idim,iBand,itetra(i,itet)) &
                                                       + 0.25*weight(i)*qal(idim,iBand,ikpt)
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !
      !---------------------------------------------------
      ! evaluate d.o.s.
      !---------------------------------------------------
      !
      g = 0.0

      DO ikpt = 1,nkpt
         DO iBand = 1,neigd

            ener = ev(iBand,ikpt)
            DO idim = 1, qdim
               w  = 0.5*wpar(idim,iBand,ikpt)
               DO ie = 1,ned
                  term(ie) = energy(ie) - ener
                  IF(energy(ie).GT.ecmax(iBand)) term(ie) = ecmax(iBand) - ener
                  IF(term(ie).LT.0.0e0)         term(ie) = 0.0e0
                  g(ie,idim) = g(ie,idim) + w * term(ie)**2
               ENDDO
            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE tetra_dos
END MODULE m_tetrados
