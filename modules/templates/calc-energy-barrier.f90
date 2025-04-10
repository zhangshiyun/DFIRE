program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_mode
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

    ! matrix
    call init_mode( mode, con )

    call make_dymatrix( mode, con )
    mode%dymatrix0 = mode%dymatrix
    call make_trimatrix( mode, con )

    call mode%solve

    mode%modez = - mode%dymatrix(:,15)

    do step=1, 1000
        call clac_modez( mode, temp1 )
        write(*,*) maxval(abs(mode%modez/mode%dymatrix(:,15)))
    end do

contains

    subroutine testvar
        implicit none

        sets%natom = 256
        sets%phi = 0.86d0
        sets%seed = 202
    end subroutine testvar

end program main
