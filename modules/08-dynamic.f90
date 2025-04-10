module mo_dynamic
    use mo_syst
    use mo_config
    implicit none

    type msd_t
        integer :: natom                                 ! Number of particles
        integer :: ndt                                   ! sample frequency
        integer :: nhist                                 ! length of histdata
        logical :: calc_flag = .false., msd_flag, fkt_flag
        real(8) :: dt, kx
        integer :: cumstep = -1
        integer :: histidx = 0
        integer :: cumcount
        real(8), allocatable, dimension(:)   :: msd      ! msd
        real(8), allocatable, dimension(:)   :: alpha2   ! non-gaussian parameter
        real(8), allocatable, dimension(:)   :: fkt      ! intemediate scatter function
        real(8), allocatable, dimension(:,:) :: histdata ! array used to stor history config
    contains
        procedure :: pmsd => print_msd
        procedure :: pfkt => print_fkt
    end type

    type vcorr_t
        integer :: natom
        integer :: ndt
        integer :: nhist
        integer :: cumstep = -1
        integer :: histidx = 0
        logical :: calc_flag = .false.
        integer :: cumcount
        real(8), allocatable, dimension(:)   :: vcorr
        real(8), allocatable, dimension(:,:) :: histdata
    end type

    type(msd_t) :: msd1, msd2, msd3
    type(vcorr_t) :: vcorr1
contains

    subroutine init_msd( con, msd, ndt, nhist, dt, kx )
        implicit none

        ! para list
        type(con_t), intent(in)    :: con
        type(msd_t), intent(inout) :: msd
        integer,     intent(in)    :: ndt, nhist
        real(8),     intent(in)    :: dt, kx

        ! readin variables
        msd%natom = con%natom
        msd%ndt   = ndt
        msd%nhist = nhist

        ! allocate memory
        if ( dt > 0.d0 ) then
            allocate( msd%msd(nhist)    )
            allocate( msd%alpha2(nhist) )
            msd%msd_flag = .true.
            msd%dt       = dt
            msd%msd      = 0.d0
            msd%alpha2   = 0.d0
        else
            msd%msd_flag = .false.
        end if
        if ( kx > 0.d0 ) then
            allocate( msd%fkt(nhist) )
            msd%fkt_flag = .true.
            msd%kx       = kx
            msd%fkt      = 0.d0
        else
            msd%fkt_flag = .false.
        end if

        if ( msd%msd_flag .or. msd%fkt_flag ) then
            allocate( msd%histdata(msd%natom,nhist) )
            msd%histdata = 0.d0
        end if

        ! initiate variables
        msd%cumstep   = -1
        msd%cumcount  = 0
        msd%calc_flag = .false.
    end subroutine

    subroutine calc_msd( con, tmsd )
        implicit none

        ! para list
        type(con_t), intent(in)    :: con
        type(msd_t), intent(inout) :: tmsd

        ! local
        integer :: idx
        integer :: i
        real(8), dimension(con%natom) :: thi, td

        associate(                       &
            ra        => con%ra,         &
            msd       => tmsd%msd,       &
            alpha2    => tmsd%alpha2,    &
            fkt       => tmsd%fkt,       &
            calc_flag => tmsd%calc_flag, &
            msd_flag  => tmsd%msd_flag,  &
            fkt_flag  => tmsd%fkt_flag,  &
            histidx   => tmsd%histidx,   &
            cumstep   => tmsd%cumstep,   &
            ndt       => tmsd%ndt,       &
            nhist     => tmsd%nhist,     &
            kx        => tmsd%kx,        &
            cumcount  => tmsd%cumcount,  &
            histdata  => tmsd%histdata   &
            )

            ! take one snapshot every ndt step
            cumstep = cumstep + 1
            if ( cumstep == ndt+1 ) cumstep = 0
            if ( cumstep /= 0 ) return

            ! main v

            histidx = histidx + 1

            ! if histidx less then nhist, then just store configure, not calc msd
            if ( .not. calc_flag .and. histidx == nhist ) calc_flag = .true.

            ! if histidx largeer than dimension of histdata, then reset it to one
            if ( histidx == nhist + 1 ) histidx = 1

            ! store config in "history data"
            histdata(:, histidx) = ra(1, :)
            thi = ra(1, :)

            if ( calc_flag ) then

                cumcount = cumcount + 1

                do i=1, nhist

                    ! if histidx >= i; idx = histidx - i + 1
                    ! if histidx <  i; idx = histidx - i + 1 + nhist
                    ! 1  2  3  4  5 !
                    !    hi         !
                    !    i          !  idx = 1
                    !          i    !  idx = 4
                    ! i             !  idx = 5
                    idx = histidx-i+1
                    if ( idx <= 0 ) idx = idx + nhist

                    td = thi - histdata(:, i)
                    if ( msd_flag ) then
                        msd(idx)    = msd(idx)    + sum( td**2 )
                        alpha2(idx) = alpha2(idx) + sum( td**4 )
                    end if
                    if ( fkt_flag ) then
                        fkt(idx) = fkt(idx) + sum( cos(td*kx) )
                    end if

                end do

            end if

        end associate
    end subroutine

    pure function print_msd( this, i ) result( re )
        implicit none

        ! para list
        class(msd_t), intent(in) :: this
        integer,      intent(in) :: i

        ! result
        real(8) :: re

        if ( this%msd_flag .and. this%calc_flag ) then
            re = this%msd(i) / this%cumcount / this%natom
        else
            re = - 1.d0
        end if
    end function

    pure function print_fkt( this, i ) result( re )
        implicit none

        ! para list
        class(msd_t), intent(in) :: this
        integer,      intent(in) :: i

        ! result
        real(8) :: re

        if ( this%fkt_flag .and. this%calc_flag ) then
            re = this%fkt(i) / this%cumcount / this%natom
        else
            re = - 1.d0
        end if
    end function

    subroutine endof_msd( tmsd )
        implicit none

        ! para list
        type(msd_t), intent(inout) :: tmsd

        associate(                    &
            natom    => tmsd%natom,   &
            cumcount => tmsd%cumcount &
            )

            if ( tmsd%msd_flag .and. tmsd%calc_flag ) then
                tmsd%msd    = tmsd%msd    / ( cumcount * natom )
                tmsd%alpha2 = tmsd%alpha2 / ( cumcount * natom )
                tmsd%alpha2 = tmsd%alpha2 / tmsd%msd**2 / 3.d0 - 1.d0
            end if

            if ( tmsd%msd_flag .and. tmsd%calc_flag ) then
                tmsd%fkt = tmsd%fkt / ( cumcount * natom )
            end if

        end associate
    end subroutine

    subroutine init_vcorr( con, vcorr, ndt, nhist )
        implicit none

        ! para list
        type(con_t),   intent(in)    :: con
        type(vcorr_t), intent(inout) :: vcorr
        integer,       intent(in)    :: ndt
        integer,       intent(in)    :: nhist

        ! readin variables
        vcorr%natom = con%natom
        vcorr%ndt   = ndt
        vcorr%nhist = nhist

        ! allocate memory
        allocate( vcorr%vcorr(nhist) )
        allocate( vcorr%histdata(con%natom, nhist) )
        vcorr%vcorr    = 0.d0
        vcorr%histdata = 0.d0

        ! initiate variables
        vcorr%cumcount  = 0
        vcorr%cumstep   = -1
        vcorr%calc_flag = .false.
    end subroutine

    subroutine calc_vcorr( con, tvcorr )
        implicit none

        ! para list
        type(con_t),   intent(in)    :: con
        type(vcorr_t), intent(inout) :: tvcorr

        ! local
        integer :: idx
        integer :: i
        real(8), dimension(con%natom) :: thi

        associate(                         &
            va        => con%va,           &
            vcorr     => tvcorr%vcorr,     &
            calc_flag => tvcorr%calc_flag, &
            histidx   => tvcorr%histidx,   &
            cumstep   => tvcorr%cumstep,   &
            ndt       => tvcorr%ndt,       &
            nhist     => tvcorr%nhist,     &
            cumcount  => tvcorr%cumcount,  &
            histdata  => tvcorr%histdata   &
            )

            cumstep = cumstep + 1
            if ( cumstep >= ndt ) cumstep = 0
            if ( cumstep /= 0 ) return

            ! main v

            histidx = histidx + 1

            ! if histidx less then nhist, then just store configure, not calc vcorr
            if ( .not. calc_flag .and. histidx == nhist ) calc_flag = .true.

            ! if histidx greater then max dim of histdata, reset to one
            if ( histidx == nhist + 1 ) histidx = 1

            ! store config to "history data"
            histdata(:, histidx) = va(1, :)
            thi = va(1, :)

            if ( calc_flag ) then

                cumcount = cumcount + 1

                do i=1, nhist

                    ! if histidx >= i; idx = histidx - i + 1
                    ! if histidx <  i; idx = histidx - i + 1 + nhist
                    ! 1  2  3  4  5 !
                    !    hi         !
                    !    i          !  idx = 1
                    !          i    !  idx = 4
                    ! i             !  idx = 5
                    idx = histidx-i+1
                    if ( idx <= 0 ) idx = idx + nhist

                    vcorr(idx) = vcorr(idx) + dot_product( histdata(:, i), thi )

                end do
            end if

        end associate
    end subroutine

    subroutine endof_vcorr( tvcorr )
        implicit none

        ! para list
        type(vcorr_t), intent(inout) :: tvcorr

        tvcorr%vcorr = tvcorr%vcorr / tvcorr%cumcount
        tvcorr%vcorr = tvcorr%vcorr / tvcorr%vcorr(1)
    end subroutine

end module
