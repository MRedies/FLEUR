!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_nocoInputCheck

   CONTAINS

   SUBROUTINE nocoInputCheck(atoms,input,sym,vacuum,noco)

      USE m_juDFT
      USE m_types_atoms
      USE m_types_input
      USE m_types_sym
      USE m_types_vacuum
      USE m_types_noco
      USE m_constants

      IMPLICIT NONE

      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_input),  INTENT(IN)    :: input
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      TYPE(t_noco),   INTENT(IN)    :: noco

      INTEGER itype
      LOGICAL l_relax_any

!---> make sure second variation is switched off
      IF (input%secvar) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'cannot be used with the second variation!!'
         CALL juDFT_error("Second variation cannot be used!!!" ,calledby="nocoInputCheck")
      END IF

!---> make sure histogram method is used
      IF (input%bz_integration==1) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'cannot be used with the Gaussian smearing for '
         WRITE (oUnit,*) 'the Brillouin zone integration!!'
         WRITE (oUnit,*) 'Please use the histogram method.'
         CALL juDFT_error("Only histogram Brillouin zone integration can be used!!!",calledby ="nocoInputCheck")
      END IF

!---> make sure force is switched off
      IF (input%l_f) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'does not support force calculations.'
         CALL juDFT_error("force calculations not supported!!!",calledby="nocoInputCheck")
      END IF

!---> make sure nstm equals zero
      IF (vacuum%nstm.NE.0) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'does not support STM calculations(nstm .NE. 0).'
         CALL juDFT_error("nstm /= 0 not supported!",calledby ="nocoInputCheck")
      END IF

!---> make sure starcoeff is switched off
!      IF (starcoeff) THEN
!         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
!         WRITE (oUnit,*) 'does not support starcoefficients output.'
!     CALL juDFT_error("starcoefficients output (for STM) cannot be !!!"
!     generated
!      ENDIF

!---> make sure coretails are switched off
      IF (input%ctail) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'cannot be used with the coretail option!! '
         CALL juDFT_error("Coretail option cannot be used!!!",calledby="nocoInputCheck")
      END IF

!---> make sure that moments are not relaxed and constrained
      l_relax_any = .FALSE.
      DO itype = 1,atoms%ntype
         l_relax_any = l_relax_any.OR.noco%l_relax(itype)
      END DO
      IF (l_relax_any.AND.noco%l_constr) THEN
         WRITE (oUnit,*)'The relaxation of the moment is switched on for at'
         WRITE (oUnit,*)'least one atom. At the same time the constrained'
         WRITE (oUnit,*)'moment option has been switched on!!!'
!          CALL juDFT_error("relaxation of moments and constraint are sw
      ENDIF
      if (l_relax_any.or.noco%l_constr) CALL judft_warn("Constraint moments and relaxations are untested in this version!")
!---> make sure that perp. component of mag. is calculated if needed
      IF ( (l_relax_any .or. noco%l_constr) .and. (.not. noco%l_mperp) ) THEN
         WRITE (oUnit,*)'The relaxation of the moment is switched on for at'
         WRITE (oUnit,*)'least one atom or the constrained moment option is'
         WRITE (oUnit,*)'switched on. In either case, you need to set'
         WRITE (oUnit,*)'l_mperp=T !!'
         CALL juDFT_error("Stop: Set l_mperp = T to relax or constrain the moments!!",calledby ="nocoInputCheck")
      ENDIF
!---> make sure l_constr is switched off in the case of spin spirals
      IF (noco%l_constr .and. noco%l_ss) THEN
         WRITE (oUnit,*)'The constraint moment option is not implemeted'
         WRITE (oUnit,*)'for spin spirals.'
         CALL juDFT_error("Stop: constraint not implemented for spin spirals!!",calledby ="nocoInputCheck")
      ENDIF

      IF(noco%l_mtnocoPot.AND.atoms%n_hia+atoms%n_u>0.AND.sym%nop.NE.1) THEN
         CALL juDFT_warn("LDA+U and FullyFullyNoco with symmetries is not correctly implemented at the moment",calledby="nocoInputCheck")
      ENDIF
      
    IF(noco%l_mtnocoPot.AND.atoms%n_hia+atoms%n_u>0.AND.(.NOT.noco%l_alignMT)) THEN
         CALL juDFT_warn("LDA+U and FullyFullyNoco should only be used together with the l_RelaxAlpha/Beta=T setting to achieve reasonable results.",calledby="nocoInputCheck")
      ENDIF
      
          !Warning on strange choice of switches
    IF (noco%l_mtNocoPot.AND..NOT.noco%l_mperp) THEN
    	CALL juDFT_error("l_mperp='F' and l_mtNocoPot='T' makes no sense.",calledby='nocoInputCheck')
    END IF
  
    IF(noco%l_alignmt.AND..NOT.(noco%l_RelaxAlpha.OR.noco%l_RelaxBeta)) CALL juDFT_warn("No relaxation is performed if neither l_RelaxBeta nor l_RelaxAlpha is true.",calledby="nocoInputCheck")
    IF(.NOT.noco%l_alignmt.AND.(noco%l_RelaxAlpha.OR.noco%l_RelaxBeta)) CALL juDFT_warn("No relaxation is performed with l_RelaxMT=F even if either l_RelaxBeta or l_RelaxAlpha is True.",calledby="nocoInputCheck")

   END SUBROUTINE nocoInputCheck

END MODULE m_nocoInputCheck
