module mo_md
    use mo_syst
    use mo_config
    use mo_list
    ! use mo_network
    use mo_force
    implicit none

    type mdargs_t
        real(8) :: temper
        real(8) :: press
    end type

    type(mdargs_t) :: mdargs

    real(8), private, parameter :: dt   = 1.d-2
    real(8), private, parameter :: dt2  = 1.d-4
    real(8), private, parameter :: hdt  = 5.d-3
    real(8), private, parameter :: hdt2 = 5.d-5

    integer, private, save      :: cumstep = 0
    integer, private, parameter :: nscale = 33

contains

    subroutine init_md( ttemper, oppress )
        implicit none

        ! para list
        real(8), intent(in) :: ttemper
        real(8), intent(in), optional :: oppress

        mdargs%temper = ttemper

        if ( present( oppress ) ) then
            mdargs%press  = oppress
        end if
    end subroutine

    subroutine pre_nvt( tcon )
        implicit none

        ! para list
        type(con_t), intent(inout) :: tcon

        associate(                 &
            natom => tcon%natom,   &
            va    => tcon%va,      &
            Tk    => mdargs%temper &
            )

            call random_number(va)
            call scale_temper( tcon, Tk )

        end associate
    end subroutine

    subroutine md_nvt( tcon, tnb )
        implicit none

        ! para list
        type(con_t),    intent(inout)        :: tcon
        type(list_t),   intent(in), optional :: tnb

        ! local
        real(8) :: chipxi

        associate(               &
            ra  => tcon%ra,      &
            va  => tcon%va,      &
            fa  => tcon%fa,      &
            Tk  => mdargs%temper &
            )

            ! velocity verlet : move 1
            va = va + fa * hdt
            ra = ra + va * dt

            ! velocity verlet : force
            call calc_force( tcon, tnb )

            ! constraint methond to control Temper
            chipxi = sum( fa*va ) / sum(va**2)
            fa = fa - chipxi * va

            ! velocity verlet : move 2
            va = va + fa * hdt

            cumstep = cumstep + 1
            if ( cumstep == nscale ) then
                cumstep = 0
                call scale_temper( tcon, Tk )
            end if

        end associate
    end subroutine

    subroutine scale_temper( tcon, tTk )
        implicit none

        ! para list
        type(con_t), intent(inout) :: tcon
        real(8), intent(in) :: tTk

        ! local
        integer :: i
        real(8) :: temp

        associate(               &
            natom => tcon%natom, &
            va    => tcon%va     &
            )

            do i=1, free
                temp    = sum( va(i,:) ) / natom
                va(i,:) = va(i,:) - temp
            end do

            temp = sum( va**2 ) / ( natom * free )
            va = sqrt( tTk / temp ) * va

        end associate
    end subroutine

end module
