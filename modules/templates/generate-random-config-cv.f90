program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_fire
    implicit none

    ! vars
    call testvar

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

    ! fire
    call init_confire( confire, con )
    call mini_fire_cv( confire )
    con = confire

    call save_config_to( con, "con.dat")

contains

    subroutine testvar
        implicit none

        sets%natom = 256
        sets%phi = 0.86d0
        sets%seed = 202
    end subroutine testvar

end program main
