program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_disorder
    use mo_dynamic
    implicit none

    ! vars
    call testvar

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_lattice_triangle( con )

    ! net
    call init_network( net, con )
    call make_network( net, con )

    ! fire
    call init_confire( confire, con, tnet=net )
    call check_system_force( confire, tnet=net )

    ! disorder
    call init_disorder( radisorder, con, sets%seed )

    ! main
    do step=10, 48

        eta = 0.01d0 * step
        testp = 1.d-6

        confire = con

        call make_disorder( confire, radisorder, eta )
        call remake_network( net, confire )

        ! B
        confire2 = confire
        confire2%la = confire2%la * ( 1.d0 - testp )
        confire2%ra = confire2%ra * ( 1.d0 - testp )

        call mini_fire_cv( confire2, tnet=net )
        call calc_Bi( confire2, net, testp )

        ! Gs
        confire2 = confire
        confire2%strain = testp
        call mini_fire_cv( confire2, tnet=net )
        call calc_gis( confire2, net, testp )

        ! Gxy
        confire2 = confire
        confire2%la(1) = confire2%la(1) * ( 1.d0 - testp )
        confire2%la(2) = confire2%la(2) * ( 1.d0 + testp )
        call mini_fire_cv( confire2, tnet=net )
        call calc_gixy( confire2, net, testp )
        write(*,*) net%mb, net%mgs, net%mgxy

        ! corr


    end do
contains

    subroutine testvar
        implicit none

        sets%natom = 256
        sets%phi = 0.92d0
        sets%seed = 202
    end subroutine testvar

end program main
