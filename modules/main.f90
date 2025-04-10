program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_mode
    implicit none

    real(8) :: istrain(150),istress(150)
    real(8) :: iEa(150),iPress(150), iG(150)
    real(8) :: bistrain(150),bistress(150)
    real(8) :: biEa(150),biPress(150), biG(150)
    real(8), allocatable, dimension(:)  :: u_tmp
    real(8), allocatable, dimension(:)  :: nonaffine
    real(8), allocatable, dimension(:,:)  :: nonaffine2
    real(8), allocatable, dimension(:,:)  :: bnonaffine
    real(8), allocatable, dimension(:,:)  :: bnonaffine2
    integer, allocatable, dimension(:,:)  :: bonds
    integer, allocatable, dimension(:,:)  :: bonds1
    integer, allocatable, dimension(:,:)  :: bonds2
    real(8) :: r_ij(2)
    real(8) :: r_ij2(2)
    real(8) :: gggtmp
    real(8) :: affine_term
    real(8) :: dstrain, dstrain_ini
    real(8) :: ptmp
    real(8) :: sdisk, theta
    real(8) :: g_tmp, g_tmp1, g_tmp2, g0, g_min, g_max
    real(8) :: bg_tmp
    integer :: round1
    integer :: round10
    integer :: round2
    integer :: countn, iseed, seed_index
    integer :: round20, n_theta
    integer :: coolsteps
    character(250) :: fconfig1
    character(250) :: configname, configname1, configname2
    character(250) :: nonaffinefile
    character(250) :: nonaffinefile2
    character(250) :: bnonaffinefile
    character(250) :: bnonaffinefile2
    character(250) :: ssfile1
    character(250) :: ssfile2
    character(250) :: z_file, g_file_1, g_file_2, g_file_0
    character(250) :: g_file_11, g_file_22, g_file_00
    character(250) :: shear_config_first, read_config_first
    character(6) :: steptemp
    character(80) :: tmodedim
    character(250) :: filesave
    character(250) :: mode_save
    character(250) :: mode_first
    character(250) :: egvalue_file
    character(250) :: egvalue_first
    character(80) :: format_tmode
    character(250) :: psith_file

    ! vars
    !call testvar
    call readvar

    allocate(u_tmp(free*sets%natom))

    allocate(nonaffine(sets%natom))

    allocate(bonds(sets%natom,sets%natom))
    allocate(bonds1(sets%natom,sets%natom))
    allocate(bonds2(sets%natom,sets%natom))

    bonds1 = 0
    bonds2 = 0

    iseed = ( iseed - 1 ) * 100

    open(101,file=egvalue_file)
    open(102,file=g_file_1)

    call init_system( con, sets%natom)

    do seed_index = 1, 100
        iseed = iseed + 1

        write(steptemp, '(i6)') iseed

        configname1 = trim(adjustl(read_config_first)) // '/2d_kalj_1.2_fast_' // trim(adjustl(steptemp)) //  '.dat'

        ! system
        call read_config( con, configname1, sets%natom)

        call init_confire(confire, con)

        call check_system_force(con)

        print *, maxval( abs(con%fa))

        con0 = con

        ! calc G from H
        call init_mode(mode, con)
        call make_dymatrix(mode,con, 0.d0)
        call mode%solve
        call mode%inv
        call calc_pw(mode)

        print*, 'shear angle...'
        do n_theta = 1, 100
            theta = pi * 0.5 / 100 * (n_theta - 1)
            call init_mode(mode0, con,opshearflag=.true.)
            call make_dymatrix(mode0, con, theta,  opflag=2)

            affine_term = mode0%dymatrix(free*sets%natom+1,free*sets%natom+1) / (con%la(1) * con%la(2))
            u_tmp = matmul(mode0%dymatrix(1:free*sets%natom,free*sets%natom+1),mode%invmatrix)
            !g_tmp = (affine_term - dot_product(u_tmp,mode0%dymatrix(free*sets%natom+1,1:free*sets%natom))) / (con%la(1) * con%la(2))
            gggtmp = dot_product(u_tmp,mode0%dymatrix(free*sets%natom+1,1:free*sets%natom)) / (con%la(1) * con%la(2))

            g_tmp = affine_term - gggtmp

            if (n_theta == 1) then
                g0    = g_tmp
                g_min = g_tmp
                g_max = g_tmp
            end if
            
            if (g_tmp < g_min) then
                g_min = g_tmp
            end if

            if (g_tmp > g_max) then
                g_max = g_tmp
            end if

        end do

        print*, 'extend...'
        call init_mode(mode1, con, opshearflag=.true., opxyflag=.true.)
        call make_dymatrix(mode1, con, 0.d0,  opflag=3)
        call mode1%solve

        do i=1,2*con%natom
            if (abs(mode%egdymatrix(i)) < 1.d-20 ) then
                mode%egdymatrix(i) = 1.d-20
            end if
            write(101,'(2es26.16)') mode%egdymatrix(i), mode%pw(i)
        end do

        write(102,'(4es26.16)') g0, g_min, g_max, minval(mode1%egdymatrix)

    end do

    close(101)
    close(102)


contains

    subroutine testvar
        implicit none

        sets%natom = 1000
        sets%seed = 2000
        sets%phi = 0.85d0

        !configname1 = '/home/sz433/project/instability/smaller_data_1024/1024_0.85_shear_fconfig_1005_increase_6.dat'
        !configname2 = '/home/sz433/project/instability/smaller_data_1024/1024_0.85_shear_fconfig_1005_increase_7.dat'
        iseed = 1
        egvalue_file = 'test_egs.dat'
        g_file_1 = 'test_g.dat'
        read_config_first = '/data/shyzhang/jamming_dos_shear_angle/LJ/lammps_config/config'


    end subroutine testvar

    subroutine readvar
        implicit none

        read(*,*) sets%natom
        read(*,*) iseed

        read(*,'(A)') g_file_1
        read(*,'(A)') egvalue_file
        read(*,'(A)') read_config_first

    !    read(*,'(A)') shear_config_first
    !    read(*,'(A)') z_file
    !    read(*,'(A)') mode_file
    !    read(*,'(A)') egvalue_file
    !    read(*,'(A)') psith_file


    end subroutine readvar


end program main

