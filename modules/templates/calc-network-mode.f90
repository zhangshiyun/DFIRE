program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_disorder
    use mo_dynamic
    use mo_mode
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
    do step=2, 48, 2

        eta = 0.01d0 * step
        testp = 1.d-6

        confire = con

        call make_disorder( confire, radisorder, eta )
        call remake_network( net, confire )

        call init_mode( mode, confire, opxyflag=.true., opshearflag=.true.)
        call make_dymatrix( mode, confire, net, opflag=3 )
        

        call init_mode( mode1, confire )
        call make_dymatrix( mode1, confire, net )
        call mode1%solve

        temp1 = 0.d0
        do i=3, mode1%mdim
            temp1 = temp1 + 1.d0 / mode1%egdymatrix(i) * dot_product( mode1%dymatrix(:,i), mode%dymatrix(1:mode1%mdim,mode1%mdim+1)) &
                                                     & * dot_product( mode1%dymatrix(:,i), mode%dymatrix(1:mode1%mdim,mode1%mdim+2) )
           !write(*,"(i8, 3es26.16)") i, mode1%egdymatrix(i), dot_product( mode1%dymatrix(:,i), mode%dymatrix(1:mode1%mdim,mode1%mdim+1) ), &
           !                                                & dot_product( mode1%dymatrix(:,i), mode%dymatrix(1:mode1%mdim,mode1%mdim+2) )
           !write(*,'(i8,2es16.6)') i, mode%dymatrix(mode%mdim-1,i), mode%dymatrix(mode%mdim,i)
           !write(*,'(i8,4es16.6)') i, mode%dymatrix(mode%mdim-2:mode%mdim,i)!, mode1%dymatrix(mode1%mdim,i)
           !write(*,'(4es26.16)') eta, mode%egdymatrix(mode%mdim-2:mode%mdim)
           !write(*,'(5es16.6)') confire%ra(:,i), confire%r(i), mode%dymatrix(free*(i-1)+1:free*i,5)
           !write(*,*) i, mode1%dymatrix(i,mode%mdim-1) / mode%dymatrix(i,mode%mdim-1)
        end do

        write(*,*) eta, mode%dymatrix(mode1%mdim+1,mode1%mdim+2), temp1
       !write(*,*)''

    end do

contains

    subroutine testvar
        implicit none

        sets%natom = 256
        sets%phi = 0.92d0
        sets%seed = 202
    end subroutine testvar

end program main
