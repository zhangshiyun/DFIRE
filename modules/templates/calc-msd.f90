program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_dynamic
    implicit none

    ! vars
    call testvar

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

    ! list
    call init_list( nb, con )
    call make_list( nb, con )

    ! md
    call init_md( 1.d-3 )
    call pre_nvt( con )

    ! msd
    call init_msd( con, msd1, 500, 200, 1.d-2, 7.2d0 )
    call init_vcorr( con, vcorr, 10, 200 )

    ! main
    do step=1, 10000
        if ( check_list( nb, con ) ) call make_list( nb, con )
        call md_nvt( con, nb )
    end do

    do step=1, 200000

        if ( check_list( nb, con ) ) call make_list( nb, con )

        call md_nvt( con, nb )

        if ( mod( step, 10000 ) == 0 ) write(*,*) step
        call calc_msd( con, msd1 )
        call calc_vcorr( con, vcorr )

    end do

    ! call save_config_to( con, "./con.dat" )

    call endof_msd( msd1 )
    call endof_vcorr( vcorr )

    open(11, file="msd.dat")
        do step=1, msd1%nhist
           write(11,*) 1.d-2 * (step-1) * msd1%ndt, msd1%msd(step), msd1%fkt(step), msd1%alpha2(step)
        end do
    close(11)

    open(11, file="vcorr.dat")
        do step=1, vcorr%nhist
            write(11,*) 1.d-2 * (step-1) * vcorr%ndt, vcorr%vcorr(step)
        end do
    close(11)
contains

    subroutine testvar
        implicit none

        sets%natom = 64
        sets%phi = 0.84d0
        sets%seed = 202
    end subroutine testvar

end program main
