program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_fire
    use mo_static
    implicit none

    ! vars
    call testvar

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )
    call init_fourier( fourier, con, tcutoff = set_fourier%cutoff, tbin = set_fourier%bin)


    do step=1, 100

        sets%seed = step
        call gen_rand_config( con, sets%seed )

        ! fire
        call init_confire( confire, con )
        call mini_fire_cv( confire )
        con = confire

        set_fourier%calc_flag = .true.
        call calc_fourier( fourier, con, set_fourier%calc_flag )

        print*, step

    end do


    call outp_fourier( fourier, set_fourier%file_sq )

contains

    subroutine testvar
        implicit none

        sets%natom = 256
        sets%phi = 0.86d0
        sets%seed = 202
    end subroutine testvar

end program main
