!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_kptgen_hybrid

   USE m_juDFT

CONTAINS

! this programm generates an aequdistant kpoint set including the
! Gamma point; it is reduced to IBZ and written in kpts (M.B.)

! odified for types D.W.

   SUBROUTINE kptgen_hybrid(input, cell, sym, kpts, l_soc)

      USE m_types
      USE m_divi

      IMPLICIT NONE

      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_kpts), INTENT(INOUT) :: kpts
! - scalars -
      LOGICAL, INTENT(IN)   ::  l_soc
! - local scalars -
      INTEGER                 ::  i, j, k, nkpt
      INTEGER                 ::  ikpt, ikpt0, nkpti
      INTEGER                 ::  nsym
! - local arrays -
      INTEGER, ALLOCATABLE   ::  rot(:, :, :), rrot(:, :, :)
      INTEGER, ALLOCATABLE   ::  invtab(:)
      INTEGER, ALLOCATABLE   ::  neqkpt(:)
      INTEGER, ALLOCATABLE   ::  pkpt(:, :, :), kptp(:), symkpt(:), iarr(:), &
                                iarr2(:)
      REAL, ALLOCATABLE      ::  rtau(:, :)
      REAL, ALLOCATABLE      ::  bk(:, :), bkhlp(:, :)
      REAL, ALLOCATABLE      ::  rarr(:)
      LOGICAL               ::  ldum

      IF (sum(kpts%nkpt3) == 0) THEN
         CALL divi(kpts%nkpt, cell%bmat, input%film, sym%nop, &
                   sym%nop2, kpts%nkpt3)
      END IF

      nkpt = kpts%nkpt3(1)*kpts%nkpt3(2)*kpts%nkpt3(3)
      ALLOCATE (bk(3, nkpt), bkhlp(3, nkpt))

      ikpt = 0
      DO i = 0, kpts%nkpt3(1) - 1
         DO j = 0, kpts%nkpt3(2) - 1
            DO k = 0, kpts%nkpt3(3) - 1
               ikpt = ikpt + 1
               bk(:, ikpt) = (/1.0*i/kpts%nkpt3(1), 1.0*j/kpts%nkpt3(2), &
               &                                      1.0*k/kpts%nkpt3(3)/)
            END DO
         END DO
      END DO

      IF (ikpt /= nkpt) call juDFT_error('failure: number of k-points')

      IF (sym%invs .OR. l_soc) THEN
         nsym = sym%nop
      ELSE
         nsym = 2*sym%nop
      END IF

      ALLOCATE (rot(3, 3, nsym), rtau(3, nsym))

      DO i = 1, sym%nop
         rot(:, :, i) = sym%mrot(:, :, i)
         rtau(:, i) = sym%tau(:, i)
      END DO

      DO i = sym%nop + 1, nsym
         rot(:, :, i) = rot(:, :, i - sym%nop)
         rtau(:, i) = rtau(:, i - sym%nop)
      END DO

      IF (any(rot(:, :, 1) - reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/)) /= 0)) &
         call juDFT_error('kptgen: First symmetry operation is not the identity.')

      ALLOCATE (rrot(3, 3, nsym), invtab(nsym))

      invtab = 0

      DO i = 1, sym%nop
         DO j = 1, sym%nop

            IF (all(matmul(rot(:, :, i), rot(:, :, j)) &
                    == reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/))) &
                .AND. all(modulo(matmul(rot(:, :, i), rtau(:, j)) + rtau(:, i), 1.0) &
                          < 1e-10)) THEN
               IF (invtab(i) /= 0) call juDFT_error('kptgen: inverse operation &
               & already defined.')
               invtab(i) = j
               rrot(:, :, i) = transpose_int(rot(:, :, j)) ! temporary fix for ifc
            END IF
         END DO
         IF (invtab(i) == 0) call judft_error( 'kptgen: inverse operation not found.')

      END DO

      DO i = sym%nop + 1, nsym
         rrot(:, :, i) = -rrot(:, :, i - sym%nop)
      END DO

      ALLOCATE (kptp(nkpt), symkpt(nkpt), rarr(3), iarr2(3), iarr(nkpt))
      ALLOCATE (pkpt(kpts%nkpt3(1) + 1, kpts%nkpt3(2) + 1, kpts%nkpt3(3) + 1))
      pkpt = 0
      DO ikpt = 1, nkpt
         iarr2 = nint(bk(:, ikpt)*kpts%nkpt3) + 1
         pkpt(iarr2(1), iarr2(2), iarr2(3)) = ikpt
      END DO

      pkpt(kpts%nkpt3(1) + 1, :, :) = pkpt(1, :, :)
      pkpt(:, kpts%nkpt3(2) + 1, :) = pkpt(:, 1, :)
      pkpt(:, :, kpts%nkpt3(3) + 1) = pkpt(:, :, 1)

      IF (any(pkpt == 0)) THEN
         CALL juDFT_error('kptgen: Definition of pkpt-pointer failed.', &
                          calledby='kptgen_hybrid')
      END IF
      iarr = 1
      ldum = .FALSE.
      DO i = 1, nkpt
         IF (iarr(i) == 0) CYCLE
         kptp(i) = i
         symkpt(i) = 1
         DO k = 2, nsym
            rarr = matmul(rrot(:, :, k), bk(:, i))*kpts%nkpt3
            iarr2 = nint(rarr)
            IF (any(abs(iarr2 - rarr) > 1e-10)) THEN
               WRITE (6, '(A,I3,A)') 'kptgen: Symmetry operation', k, &
                  ' incompatible with k-point set.'
               ldum = .TRUE.
            END IF
            iarr2 = modulo(iarr2, kpts%nkpt3) + 1
            IF (any(iarr2 > kpts%nkpt3)) &
               call judft_error( 'kptgen: pointer indices exceed pointer dimensions.')
            j = pkpt(iarr2(1), iarr2(2), iarr2(3))
            IF (j == 0) call judft_error( 'kptgen: k-point index is zero (bug?)')
            IF (iarr(j) == 0 .OR. j == i) CYCLE
            iarr(j) = 0
            kptp(j) = i
            symkpt(j) = k
         END DO
      END DO
      IF (ldum) &
         call judft_error( 'kptgen: Some symmetry operations are incompatible '&
                          //'with k-point set.')
      i = 0
      DO ikpt = 1, nkpt
         IF (iarr(ikpt) == 1) THEN
            i = i + 1
            iarr(ikpt) = i
         END IF
      END DO
      nkpti = i
      DO ikpt = 1, nkpt
         IF (iarr(ikpt) == 0) THEN
            i = i + 1
            iarr(ikpt) = i
         END IF
      END DO
      bk(:, iarr) = bk
      kptp = iarr(kptp)
      kptp(iarr) = kptp
      symkpt(iarr) = symkpt
      DO i = 1, kpts%nkpt3(1) + 1
         DO j = 1, kpts%nkpt3(2) + 1
            DO k = 1, kpts%nkpt3(3) + 1
               pkpt(i, j, k) = iarr(pkpt(i, j, k))
            END DO
         END DO
      END DO

      ALLOCATE (neqkpt(nkpti))
      neqkpt = 0
      DO ikpt0 = 1, nkpti
         DO ikpt = 1, nkpt
            IF (kptp(ikpt) == ikpt0) neqkpt(ikpt0) = neqkpt(ikpt0) + 1
         END DO
      END DO

!     Do not do any IO, but store in kpts
      kpts%nkpt = nkpti
      if (allocated(kpts%bk)) deallocate (kpts%bk)
      if (allocated(kpts%wtkpt)) deallocate (kpts%wtkpt)
      ALLOCATE (kpts%bk(3, kpts%nkpt), kpts%wtkpt(kpts%nkpt))

      DO ikpt = 1, nkpti
         kpts%bk(:, ikpt) = bk(:, ikpt)
         kpts%wtkpt(ikpt) = neqkpt(ikpt)
      END DO
      kpts%posScale = 1.0

   CONTAINS

! Returns least common multiple of the integers iarr(1:n).
      FUNCTION kgv(iarr, n)
         IMPLICIT NONE
         INTEGER              :: kgv
         INTEGER, INTENT(IN)  :: n, iarr(n)
         LOGICAL              :: lprim(2:maxval(iarr))
         INTEGER, ALLOCATABLE :: prim(:), expo(:)
         INTEGER              :: nprim, marr
         INTEGER              :: i, j, ia, k
! Determine prime numbers
         marr = maxval(iarr)
         lprim = .TRUE.
         DO i = 2, marr
            j = 2
            DO WHILE (i*j <= marr)
               lprim(i*j) = .FALSE.
               j = j + 1
            END DO
         END DO
         nprim = count(lprim)
         ALLOCATE (prim(nprim), expo(nprim))
         j = 0
         DO i = 2, marr
            IF (lprim(i)) THEN
               j = j + 1
               prim(j) = i
            END IF
         END DO
! Determine least common multiple
         expo = 0
         DO i = 1, n
            ia = iarr(i)
            IF (ia == 0) CYCLE
            DO j = 1, nprim
               k = 0
               DO WHILE (ia/prim(j)*prim(j) == ia)
                  k = k + 1
                  ia = ia/prim(j)
               END DO
               expo(j) = max(expo(j), k)
            END DO
         END DO
         kgv = 1
         DO j = 1, nprim
            kgv = kgv*prim(j)**expo(j)
         END DO
         DEALLOCATE (prim, expo)
      END FUNCTION kgv

!     ifc seems to have problems transposing integer arrays. this is a fix.
      FUNCTION transpose_int(a)
         IMPLICIT NONE
         integer :: transpose_int(3, 3), a(3, 3)
         integer :: i, j
         DO i = 1, 3
            DO j = 1, 3
               transpose_int(i, j) = a(j, i)
            END DO
         END DO
      END FUNCTION transpose_int
   END SUBROUTINE kptgen_hybrid

END MODULE m_kptgen_hybrid