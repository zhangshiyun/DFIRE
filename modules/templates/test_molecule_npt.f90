program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_molecule
    use mo_extra_molecule
    use mo_molecule_npt
    implicit none
    
    type( tpmolecule_npt ) :: test_npt
    type( tpset_npt ) :: set_npt
    type( tpout_molecule ) :: out_npt    
    
    ! vars
    call testvar

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

    ! list
    call init_list( nb, con )
    call make_list( nb, con )

    !npt   
    set_npt%dt        = 1.d-2
    set_npt%temper    = 2.d-3
    set_npt%pre       = 2.d-2
    set_npt%resc_prob = 1.d-4   

    call test_npt%initial( con, nb, calc_extra_force, set_npt )  

    do step=1, 1000000

        call test_npt%body( con, nb, calc_extra_force, out_npt )   

        if( mod(step, 10000) .eq. 0 ) then
            set_npt%temper = set_npt%temper - 1.d-5 ! simulate the cooling (from 2e-3 to 1e-3)
            call test_npt%initial( con, nb, calc_extra_force, set_npt ) 
        end if

        if( mod(step, 7777) .eq. 0 ) then
            print*, step, out_npt%T, out_npt%press, out_npt%phi
        end if   

    end do

    call test_npt%clean    

contains

    subroutine testvar
        implicit none

        sets%natom = 128
        sets%phi = 0.84d0
        sets%seed = 202

    end subroutine testvar

end program main
