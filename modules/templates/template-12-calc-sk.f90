program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_fire
    use mo_structure
    implicit none

    type(sk_t) :: sk1

    ! vars
    call testvar()

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

    ! initiate sk1
    call init_sk(con, sk1)

    do step=1, 2000

        print*, "#", step

        sets%seed = step
        call gen_rand_config( con, sets%seed )

        ! fire
        call init_confire( confire, con )
        call mini_fire_cv( confire )
        con = confire

        call calc_sk(con, sk1)

    end do


    call endofsk(sk1)

    do i=1, sk1%nmax
        print*, sk1%ks(free+1, i), sk1%sks(i)
    end do


contains

    subroutine testvar
        implicit none

        sets%natom = 1024
        sets%phi = 0.9d0
        sets%seed = 202
    end subroutine testvar

    subroutine readvar()
        implicit none

        read(*,*) sets%natom

    end subroutine

end program main
