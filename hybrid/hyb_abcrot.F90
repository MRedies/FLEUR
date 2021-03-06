MODULE m_hyb_abcrot
CONTAINS
   SUBROUTINE hyb_abcrot(hybinp, atoms, neig, sym,&
                    acof, bcof, ccof)
!     ***************************************************************
!     * This routine transforms a/b/cof which are given wrt rotated *
!     * MT functions (according to invsat/ngopr) into a/b/cof wrt   *
!     * unrotated MT functions. Needed for GW calculations.         *
!     *                                                             *
!     * Christoph Friedrich Mar/2005                                *
!     ***************************************************************
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_atoms), INTENT(IN)  :: atoms
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: neig
!     ..
!     .. Array Arguments ..

      COMPLEX, INTENT(INOUT) :: acof(:, 0:, :) !(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd)
      COMPLEX, INTENT(INOUT) :: bcof(:, 0:, :) !(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd)
      COMPLEX, INTENT(INOUT) :: ccof(-atoms%llod:, :, :, :)!(-llod:llod,input%neig,atoms%nlod,atoms%natd)
!     ..
!     .. Local Scalars ..
      INTEGER itype, ineq, iatom, iop, ilo, i, l, ifac
!     ..
!     .. Local Arrays ..
!***** COMPLEX, ALLOCATABLE :: d_wgn(:,:,:,:) !put into module m_savewigner
!

      IF (.NOT. ALLOCATED(hybinp%d_wgn2)) THEN    !calculate sym%d_wgn only once
         PRINT *, "calculate wigner-matrix"
         call judft_error('WIGNER MATRIX should be available in hybinp part')
         !IF (.NOT.oneD%odi%d1) THEN
         !  allocate(sym%d_wgn(-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,sym%nop))
         !  CALL d_wigner(sym%nop,sym%mrot,cell%bmat,atoms%lmaxd,sym%d_wgn)
         !ELSE
         !  allocate(sym%d_wgn(-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,oneD%ods%nop))
         !  CALL d_wigner(oneD%ods%nop,oneD%ods%mrot,cell%bmat,atoms%lmaxd,sym%d_wgn)
         !ENDIF
      ENDIF

      iatom = 0
      DO itype = 1, atoms%ntype
         DO ineq = 1, atoms%neq(itype)
            iatom = iatom + 1
            iop = sym%ngopr(iatom)
!                                    l                        l    l
! inversion of spherical harmonics: Y (pi-theta,pi+phi) = (-1)  * Y (theta,phi)
!                                    m                             m
            ifac = 1
            IF (sym%invsat(iatom) == 2) THEN
               iop = sym%ngopr(sym%invsatnr(iatom))

               ifac = -1
            ENDIF
            DO l = 1, atoms%lmax(itype)
!  replaced d_wgn by conjg(d_wgn),FF October 2006
               DO i = 1, neig
                  acof(i, l**2:l*(l + 2), iatom) = ifac**l*matmul(conjg(hybinp%d_wgn2(-l:l, -l:l, l, iop)),&
                                                                  acof(i, l**2:l*(l + 2), iatom))
                  bcof(i, l**2:l*(l + 2), iatom) = ifac**l*matmul(conjg(hybinp%d_wgn2(-l:l, -l:l, l, iop)),&
                                                                  bcof(i, l**2:l*(l + 2), iatom))
               ENDDO
            ENDDO
            DO ilo = 1, atoms%nlo(itype)
               l = atoms%llo(ilo, itype)
               IF (l > 0) THEN
                  DO i = 1, neig
                     ccof(-l:l, i, ilo, iatom) = ifac**l*matmul(conjg(hybinp%d_wgn2(-l:l, -l:l, l, iop)), &
                                                                ccof(-l:l, i, ilo, iatom))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE hyb_abcrot
END MODULE m_hyb_abcrot
