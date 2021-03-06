MODULE m_fleurinput_postprocess
  USE m_types_fleurinput
CONTAINS
  SUBROUTINE fleurinput_postprocess(Cell,Sym,Atoms,Input,Noco,Vacuum,&
    Banddos,Oned,Xcpot,Kpts,gfinp)
    USE m_juDFT
    USE m_types_fleurinput
    use m_make_sym
    USE m_chkmt
    use m_make_xcpot
    use m_lapwdim
    use m_checks
    USE m_relaxio
    USE m_types_nococonv

    TYPE(t_cell),INTENT(INOUT)  ::cell
    TYPE(t_sym),INTENT(INOUT)   ::sym
    TYPE(t_atoms),INTENT(INOUT) ::atoms
    TYPE(t_input),INTENT(INOUT) ::input
    TYPE(t_noco),INTENT(IN)     ::noco
    TYPE(t_vacuum),INTENT(INOUT)::vacuum
    TYPE(t_banddos),INTENT(IN)  ::banddos
    TYPE(t_oneD),INTENT(INOUT)  ::oneD
    CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT)::xcpot
    TYPE(t_kpts),INTENT(IN)     ::kpts
    TYPE(t_gfinp),INTENT(IN)    ::gfinp

    call cell%init(DOT_PRODUCT(atoms%volmts(:),atoms%neq(:)))
    call atoms%init(cell)
    CALL sym%init(cell,input%film)
    call vacuum%init(sym)

    CALL make_sym(sym,cell,atoms,noco,oneD,input,gfinp)
    call make_xcpot(xcpot,atoms,input)
    call oneD%init(atoms)

    call check_input_switches(banddos,vacuum,noco,atoms,input,sym)
    ! Check muffin tin radii, only checking, dont use new parameters
    CALL chkmt(atoms,input,vacuum,cell,oneD,.TRUE.)
    !adjust positions by displacements
    CALL apply_displacements(cell,input,vacuum,oneD,sym,noco,atoms,gfinp)


  END SUBROUTINE fleurinput_postprocess
END MODULE m_fleurinput_postprocess
