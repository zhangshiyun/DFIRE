module mo_static
    use mo_syst
    use mo_config
    use mo_math
    implicit none

    type tpset_fourier
        character(250)  ::  file_sq   = "short_sq.dat" ! file for the structure factor s(|q|)
        logical         ::  calc_flag = .true.         ! to calculate or not, the condition is assigned in program main
        ! control parameters for the output type, one could also adjust them in program main
        real(8)         ::  cutoff  = 3.d0             ! the maximal qsc in the calculation and the output
        real(8)         ::  bin     = 0.5d0            ! the bin of qsc for average
    end type

    type(tpset_fourier) :: set_fourier

    type tpfourier
        ! output: dimensionless qualities, more qualities in q space can be added here
        real(8), allocatable, dimension(:) :: qsc
        real(8), allocatable, dimension(:) :: sqsc   ! the structure factor s(|q|)
        ! output: control parameter
        integer, allocatable, dimension(:) :: qord   ! the ordered index of qsc
        logical, allocatable, dimension(:) :: q_flag ! for cutoff
        ! poilot process: control parameters, do not attach values here!
        integer :: ave_times                         ! to record the number of effective call for calculation
        integer :: ave_num                           ! ave_num = ave_times * natom
        ! poilot process: wave vector, being calculated in the init_fourier
        logical :: bin_flag = .false.
        logical :: cut_flag = .false.
        real(8) :: unit_d, cutoff, bin
        real(8) :: kve(free)                         ! kve(i) = 2.0*pi / la(i), dimensional quality, la(i) is the box length
        integer :: nve                               ! nve = int(dble(natom)**(1.d0/dble(free))), natom is the number of particles
        integer :: dimq                              ! the len of qsc/sqsc/qord (it is according to the set of cutoff, if cutoff presents)
        contains
            ! the initial procedure
            ! procedure :: init_fourier     => init_fourier
            ! calculation
            ! procedure :: calc_fourier     => calc_fourier
            procedure :: calc_unit_length => calc_unit_length
            ! the output procedure for the structure factor s(|q|)
            ! procedure :: outp_fourier     => outp_fourier
    end type

    type(tpfourier) :: fourier

contains

    ! free == 2

    subroutine init_fourier( tfourier, tcon, tcutoff, tbin )
        !
        ! initiate fourier, allocate memory of fourier
        !
        implicit none

        ! para list
        type(tpfourier),     intent(inout) :: tfourier
        type(con_t),         intent(in)    :: tcon
        real(8), intent(in), optional      :: tcutoff
        real(8), intent(in), optional      :: tbin

        ! local
        real(8), allocatable, dimension(:) ::  qtest
        real(8)  :: temp_qsc
        integer  :: dimqtest
        integer  :: nn(free)
        integer  :: i, j, countn, counti

        if(free .ne. 2) then
            write(*,*) "the module is only for 2d, modify it, if you need 3d version"
            stop
        end if

        associate(                            &
            kve       => tfourier%kve,        &
            nve       => tfourier%nve,        &
            dimq      => tfourier%dimq,       &
            unit_d    => tfourier%unit_d,     &
            cutoff    => tfourier%cutoff,     &
            cut_flag  => tfourier%cut_flag,   &
            bin       => tfourier%bin,        &
            bin_flag  => tfourier%bin_flag,   &
            ave_times => tfourier%ave_times,  &
            ave_num   => tfourier%ave_num,    &
            natom     => tcon%natom,          &
            la        => tcon%la              &
            )

        kve = 2.d0*pi / la
        nve = int( dble(natom)**( 1.d0/dble(free) ) )
        unit_d = tfourier%calc_unit_length( tcon )

        dimqtest = 1
        do i = 1, free
            dimqtest = dimqtest * ( nve-i+1 )
        end do

        dimqtest = dimqtest / int(gamma( dble(free)+1.d0 ))

        if ( present( tbin ) ) then
            bin_flag = .true.
            bin = tbin
        end if

        if ( present( tcutoff ) ) then
            cut_flag = .true.
            cutoff = tcutoff
            allocate( qtest(dimqtest) )
            countn = 0
            nn(1) = 0
            do while (nn(1) < nve)
                nn(1) = nn(1) + 1
                nn(2) = nn(1)
                do while (nn(2) < nve)
                    nn(2) = nn(2) + 1
                    countn  = countn + 1
                    qtest(countn) = dsqrt( sum( ( dble(nn)/la )**2 ) ) * 2.d0*pi * unit_d
                end do
            end do

            call qsort(qtest)

            loop1: do i = 1, dimqtest
                if(qtest(i) >= cutoff) then
                    dimq = i - 1
                    exit loop1
                end if
            end do loop1
            deallocate( qtest )
        else
            dimq = dimqtest
        end if

        if(dimq < 1) then
            write(*,*) "the cutoff is too larger, please reenter"
            stop
        end if


        allocate( tfourier%qsc (dimq) )
        allocate( tfourier%sqsc(dimq) )
        allocate( tfourier%qord(dimq) )
        allocate( tfourier%q_flag (dimqtest) )

        associate(                            &
            qsc       => tfourier%qsc,        &
            sqsc      => tfourier%sqsc,       &
            qord      => tfourier%qord,       &
            q_flag    => tfourier%q_flag      &
        )

        qsc    = 0.d0
        sqsc   = 0.d0
        qord   = 0
        q_flag = .false.

        if( present( tcutoff ) ) then
            countn = 0
            counti = 0
            nn(1) = 0
            loop1: do while (nn(1) < nve)
                nn(1) = nn(1) + 1
                nn(2) = nn(1)
                loop2: do while (nn(2) < nve)
                    nn(2) = nn(2) + 1
                    counti  = counti + 1
                    temp_qsc = dsqrt( sum( ( dble(nn)/la )**2 ) ) * 2.d0*pi * unit_d
                    if(temp_qsc >= cutoff) cycle
                    countn  = countn + 1
                    if (countn > dimq) exit loop1
                    qsc(countn) = temp_qsc
                    q_flag(counti) = .true.
                end do loop2
            end do loop1
        else
            countn = 0
            nn(1) = 0
            do while (nn(1) < nve)
                nn(1) = nn(1) + 1
                nn(2) = nn(1)
                do while (nn(2) < nve)
                    nn(2) = nn(2) + 1
                    countn  = countn + 1
                    qsc(countn) = dsqrt( sum( ( dble(nn)/la )**2 ) ) * 2.d0*pi * unit_d
                end do
            end do
        end if

        ave_num   = natom
        ave_times = 0
        qord      = sortperm(dimq, qsc)


        ! do i = 1, dimq
        !     write(*,*) i, qsc(i), qord(i)
        ! end do

        end associate
        end associate
    end subroutine

    subroutine calc_fourier( tfourier, tcon, tcalc_flag )
        !
        ! if tcalc_flag is true, calculate fourier
        !
        implicit none

        ! para list
        type(tpfourier), intent(inout) :: tfourier
        type(con_t),     intent(in)    :: tcon
        logical,         intent(in)    :: tcalc_flag

        ! local
        integer :: nn(free)
        real(8) :: kk(free)
        integer :: i, countn, counti
        real(8) :: cosqra, sinqra, qra, temp
        real(8) :: rsq, isq

        if(free .ne. 2) then
            write(*,*) "the module is only for 2d, modify it, if you need 3d version"
            stop
        end if

        if(tcalc_flag) then

            associate(                           &
                q_flag    => tfourier%q_flag,    &
                sqsc      => tfourier%sqsc,      &
                kve       => tfourier%kve,       &
                nve       => tfourier%nve,       &
                dimq      => tfourier%dimq,      &
                unit_d    => tfourier%unit_d,    &
                cutoff    => tfourier%cutoff,    &
                cut_flag  => tfourier%cut_flag,  &
                ave_times => tfourier%ave_times, &
                natom     => tcon%natom,         &
                la        => tcon%la,            &
                ra        => tcon%ra             &
                )

            if(cut_flag) then
                counti = 0
                countn = 0
                nn(1) = 0
                loop1: do while (nn(1) < nve)
                    nn(1) = nn(1) + 1
                    nn(2) = nn(1)
                    loop2: do while (nn(2) < nve)
                        nn(2) = nn(2) + 1
                        counti = counti + 1
                        if(q_flag(counti)) then
                            countn  = countn + 1
                            if(countn > dimq) exit loop1
                            rsq = 0.d0
                            isq = 0.d0
                            kk = dble(nn) * kve
                            do i = 1, natom
                                qra = sum(kk * ra(:,i))
                                cosqra = cos(qra)
                                sinqra = sin(qra)
                                rsq = rsq + cosqra
                                isq = isq + sinqra
                            end do
                            sqsc(countn) = sqsc(countn) + rsq**2 + isq**2
                        end if
                    end do loop2
                end do loop1
            else
                countn = 0
                nn(1) = 0
                do while (nn(1) < nve)
                    nn(1) = nn(1) + 1
                    nn(2) = nn(1)
                    do while (nn(2) < nve)
                        nn(2) = nn(2) + 1
                        countn  = countn + 1
                        rsq = 0.d0
                        isq = 0.d0
                        kk = dble(nn) * kve
                        do i = 1, natom
                            qra = sum(kk * ra(:,i))
                            cosqra = cos(qra)
                            sinqra = sin(qra)
                            rsq = rsq + cosqra
                            isq = isq + sinqra
                        end do
                        sqsc(countn) = sqsc(countn) + rsq**2 + isq**2
                    end do
                end do
            end if

            ave_times = ave_times + 1

            end associate
        end if
    end subroutine

    subroutine outp_fourier( tfourier, tfile_sq )
        !
        ! output s(|q|), if tbin presents, average s(|q|) according to tbin.
        !
        implicit none

        ! para list
        type(tpfourier),     intent(inout) :: tfourier
        character(250),      intent(inout) :: tfile_sq

        ! local
        real(8), allocatable, dimension(:)  :: out_qsc
        real(8), allocatable, dimension(:)  :: out_sqsc ! the structure factor s(|q|)
        integer  :: i, counti

        real(8)  :: maxqsc, temp_sqsc, temp_qsc


        associate(                            &
            qsc       => tfourier%qsc,        &
            sqsc      => tfourier%sqsc,       &
            qord      => tfourier%qord,       &
            nve       => tfourier%nve,        &
            dimq      => tfourier%dimq,       &
            bin       => tfourier%bin,        &
            bin_flag  => tfourier%bin_flag,   &
            ave_times => tfourier%ave_times,  &
            ave_num   => tfourier%ave_num     &
            )

        if(ave_times < 1) then
            print *, "fail to output, no calculation of fourier is done"
            return
        end if

        ave_num = ave_num * ave_times
        sqsc = sqsc / dble(ave_num)
        allocate( out_qsc(dimq) )
        allocate( out_sqsc(dimq) )
        out_qsc = 0.d0
        out_sqsc = 0.d0

        do i = 1, dimq
             out_qsc (i) = qsc( qord(i) )
             out_sqsc(i) = sqsc( qord(i) )
        end do

        if ( bin_flag ) then
            maxqsc = out_qsc(1) * 0.9d0
            counti = 0
            temp_qsc = 0.d0
            temp_sqsc = 0.d0
            open(2, file = tfile_sq)
            do i = 1, dimq
                if(out_qsc(i) > maxqsc) then
                    if(i > 1) then
                        temp_qsc  = temp_qsc  / dble(counti)
                        temp_sqsc = temp_sqsc / dble(counti)
                        write(2,*) temp_qsc, temp_sqsc
                    end if
                    counti    = 0
                    temp_qsc  = 0.d0
                    temp_sqsc = 0.d0
                    maxqsc = maxqsc + bin
                end if
                temp_qsc  = temp_qsc  + out_qsc(i)
                temp_sqsc = temp_sqsc + out_sqsc(i)
                counti    = counti + 1
            end do
            close(2)
        else
            open(1, file = tfile_sq)
            do i = 1, dimq
                write(1,*) out_qsc(i), out_sqsc(i)
            end do
            close(1)
        end if

        deallocate( out_qsc )
        deallocate( out_sqsc )

        end associate
    end subroutine

    pure function calc_unit_length( this, tcon ) result(unit_d)
        !
        ! unit_d is the average diameter, the unit of length
        !
        implicit none

        ! para list
        class(tpfourier), intent(in) :: this
        type(con_t),      intent(in) :: tcon

        ! result
        real(8) :: unit_d

        ! local
        real(8) :: sdisk

        associate(                  &
            r      => tcon%r,       &
            natom  => tcon%natom    &
            )

        sdisk  = sum((2.d0*r)**free)
        sdisk  = sdisk / dble(natom)
        unit_d = sdisk ** (1.d0 / dble(free))

        end associate
    end function

end module
