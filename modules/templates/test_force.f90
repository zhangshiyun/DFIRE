program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_mode
    implicit none

    procedure(abstract_fun_force), pointer :: calc_fun_force_h => calc_fun_force

    ! vars
    call testvar ()

    !calc_fun_force_h => calc_fun_force

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

    ! fire
    call init_confire( confire, con )
    call mini_fire_cv( confire, force_type = calc_fun_force_h )
    con = confire

contains

    subroutine testvar ()
        implicit none

        sets%natom = 256
        sets%phi = 0.86d0
        sets%seed = 202

    end subroutine testvar 

end program main
