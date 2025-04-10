program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_fire
    use mo_field
    use mo_rs_corr
    implicit none

    type( tpfield_sca ) :: field_dens
    type( tpset_corr )  :: set_ls_gr,  set_ll_gr
    type( tpcorr_sca )  :: calc_ls_gr, calc_ll_gr
    type( tpout_corr )  :: out_ls_gr,  out_ll_gr
          
! vars
    call testvar

! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

!initialize the density field
    field_dens%num     =  con%natom
    field_dens%ra      => con%ra
    field_dens%la      => con%la
    call field_dens%init_field() 
    field_dens%ary     = 1.d0
    field_dens%unit_d  = 1.d0

!the settings of correlation
    set_ls_gr%unit_d   = (1.4 + 1.0)*0.5
    set_ll_gr%unit_d   = 1.4
    set_ls_gr%rbin     = con%la(1) * 1e-4
    set_ll_gr%rbin     = con%la(1) * 1e-4
    set_ls_gr%rmin     = 0.8
    set_ll_gr%rmin     = 0.8
    set_ls_gr%rmax     = con%la(1) * 0.2
    set_ll_gr%rmax     = con%la(1) * 0.2
    set_ls_gr%idi_num  = con%natom / 2
    set_ls_gr%idj_num  = con%natom / 2
    set_ll_gr%idi_num  = con%natom / 2
    set_ll_gr%idj_num  = con%natom / 2
    set_ls_gr%filename = 'gr_ls.dat'
    set_ll_gr%filename = 'gr_ll.dat'
    call set_ls_gr%init_setcorr()
    call set_ll_gr%init_setcorr()
    do i = 1, set_ls_gr%idi_num
        set_ls_gr%idi(i) = i    
    end do
    do j = 1, set_ls_gr%idj_num
        set_ls_gr%idj(j) = con%natom / 2 + j    
    end do
    do i = 1, set_ll_gr%idi_num
        set_ll_gr%idi(i) = con%natom / 2 + i    
    end do
    do j = 1, set_ll_gr%idj_num
        set_ll_gr%idj(j) = con%natom / 2 + j    
    end do

!initialize the computing of correlation    
    call calc_ls_gr%initial( set_ls_gr ) 
    call calc_ll_gr%initial( set_ll_gr ) 
    
    do step=1, 1000
        sets%seed = step
        call gen_rand_config( con, sets%seed )
! fire
        call init_confire( confire, con )
        call mini_fire_cv( confire )
        con = confire
! comput gr
        field_dens%ra      => con%ra
        field_dens%la      => con%la
        set_ls_gr%calc_flag = .true.
        set_ll_gr%calc_flag = .true.
        call calc_ls_gr%body( field_dens, field_dens, set_ls_gr%calc_flag )
        call calc_ll_gr%body( field_dens, field_dens, set_ll_gr%calc_flag )

        print*, step
    end do

!do i = 1, con%natom
!    print*, con%ra(i,1), field_dens%ra(i,1)
!enddo

!output the results
    call calc_ls_gr%output( out_ls_gr, set_ls_gr%filename )
    call calc_ll_gr%output( out_ll_gr, set_ll_gr%filename )
!visit the results
    step = size(out_ls_gr%absr)
    do i = 1, step
        print*, i, out_ls_gr%absr(i), out_ls_gr%corr(i)
    end do

!clean everything
    call field_dens%clean_field() 
    call calc_ls_gr%clean()
    call calc_ll_gr%clean()
    call set_ls_gr%clean_setcorr()
    call set_ll_gr%clean_setcorr()
    call out_ls_gr%clean_outcorr()
    call out_ll_gr%clean_outcorr()

contains

    subroutine testvar
        implicit none

        sets%natom = 256
        sets%phi = 0.86d0
        sets%seed = 202
    end subroutine testvar

end program main
