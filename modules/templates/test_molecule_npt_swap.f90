program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_molecule
    use mo_extra_molecule
    use mo_molecule_npt
    use mo_molecule_npt_swap
    implicit none
    
    type( tpmolecule_npt_swap ) :: test_npt_swap
    type( tpset_npt_swap )      :: set_npt_swap
    type( tpout_molecule )      :: out_npt_swap    
    
    ! vars
    call testvar

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

    ! list
    call init_list( nb, con )
    call make_list( nb, con )

    !npt   
    set_npt_swap%dt        = 1.d-2
    set_npt_swap%temper    = 3.d-3
    set_npt_swap%pre       = 3.d-1
    set_npt_swap%resc_prob = 3.d-4   
    set_npt_swap%swap_prob = 1.d-2  

    call test_npt_swap%initial( con, nb, calc_extra_force, set_npt_swap )  

    do step = 1, 10000000

        call test_npt_swap%body( con, nb, calc_extra_force, out_npt_swap )   

        if( mod(step, 100000) .eq. 0 ) then
            print*, step, out_npt_swap%T, out_npt_swap%Ea, out_npt_swap%press, out_npt_swap%phi
        end if   

    end do 

    call test_npt_swap%clean  

    do i = 1, con%natom
        con%ra(:,i) = con%ra(:,i) - anint( con%ra(:,i) / con%la ) * con%la
        con%r(i)    = con%r(i) * 0.7
    end do

    call save_config_to( con, "con.dat")

contains

    subroutine testvar
        implicit none

        sets%natom = 1024
        sets%phi   = 1.54d0
        sets%seed  = 202

    end subroutine testvar

end program main
