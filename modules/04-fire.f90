module mo_fire
    !
    !  FIRE method of local minimization, see more at paper https://doi.org/10.1103/PhysRevLett.97.170201
    !
    use mo_syst
    use mo_config
    use mo_list
    use mo_network
    use mo_force
    implicit none

    !-- fire var
    real(8), private, parameter :: fmax    = 1.d-11
    integer, private, parameter :: stepmax = 5e8
    integer, private :: step
    !--
    real(8), private, parameter :: dt0   = 1.d-5
    real(8), private, parameter :: dtmax = 1.d-3
    real(8), private, parameter :: beta0 = 1.d-1
    real(8), private, parameter :: finc  = 1.1d0
    real(8), private, parameter :: fdec  = 0.5d0
    real(8), private, parameter :: fbeta = 0.99d0
    integer, private, parameter :: nmin  = 5
    !--
    real(8), private :: dt, dt2, dt22, beta, power, fn, vn
    integer, private :: count


    type(con_t)  :: confire, confire1, confire2, bconfire
    type(list_t) :: nbfire

    procedure(abstract_force), pointer :: calc_force_h => null()

contains

    subroutine init_confire( tconfire, tcon, tnet )
        implicit none

        ! para list
        type(con_t),     intent(inout)        :: tconfire
        type(con_t),     intent(in)           :: tcon
        type(network_t), intent(in), optional :: tnet

        ! copy con to confire
        tconfire = tcon

        ! if net don't exist, allocate list
        if ( .not. present( tnet ) ) then
            call init_list( nbfire, tconfire )
        end if

        calc_force_h => calc_force_lj
    end subroutine

    subroutine check_system_force( tcon, tnet )
        implicit none

        ! para list
        type(con_t),     intent(inout)           :: tcon
        type(network_t), intent(inout), optional :: tnet

        if ( .not. present( tnet ) ) then
            call make_list( nbfire, tcon )
            call calc_force_h( tcon, nbfire )
        else
           !call make_network( tnet, tcon )
            call calc_force_spring( tcon, tnet )
        end if
    end subroutine

    ! main 1
    subroutine mini_fire_cv( tcon, tnet, force_type)
        implicit none

        ! para list
        type(con_t),     intent(inout)           :: tcon
        type(network_t), intent(in),    optional :: tnet

        logical, external, optional :: force_type
        logical :: get_force

        ! local
        logical :: nonnetwork_flag
        real(8) :: dt, beta, temp
        real(8) :: onembeta, betavndfn
        integer :: cumn

        ! network system or not
        nonnetwork_flag = .true.
        if ( present( tnet ) ) nonnetwork_flag = .false.

        associate(              &
            ra    => tcon%ra,   &
            va    => tcon%va,   &
            fa    => tcon%fa,   &
            natom => tcon%natom &
            )

            ! initial sets
            fa   = 0.d0 ; va   = 0.d0  ;
            dt   = dt0  ; beta = beta0 ;
            cumn = 0

            ! calc fortran before iteration
            if ( nonnetwork_flag ) then
                call make_list( nbfire, tcon )
                if (present( force_type )) then
                    get_force = force_type( tcon, nbfire )
                else
                    call calc_force_h( tcon, nbfire )
                end if
            else
                call calc_force_spring( tcon, tnet )
            end if

            ! main
            do step=1, stepmax

                dt2  = 0.5d0 * dt
                dt22 = 0.5d0 * dt**2

                ! velocity verlet method / move 1
                ra = ra + va * dt + fa * dt22
                va = va + fa * dt2

                ! check list
                !if ( nonnetwork_flag ) then
                !    if( check_list( nbfire, tcon ) ) then
                !        call make_list( nbfire, tcon )
                !    end if
                !end if
                call make_list(nbfire, tcon)

                ! velocity verlet method / force
                if ( nonnetwork_flag ) then
                    if (present( force_type )) then
                        get_force = force_type( tcon, nbfire )
                    else
                        call calc_force_h( tcon, nbfire )
                    end if
                else
                    call calc_force_spring( tcon, tnet )
                end if

                ! velocity verlet method / move 2
                va = va + fa * dt2

                ! fire
                cumn  = cumn + 1
                power = sum( fa * va )

                fn = sqrt(sum(fa**2))
                vn = sqrt(sum(va**2))

                onembeta = 1.d0 - beta
                betavndfn = beta*vn/fn
                va  = onembeta * va + betavndfn * fa
                !va = ( 1.d0 - beta ) * va + beta * (vn/fn) * fa

                if ( power > 0.d0 .and. cumn > nmin ) then
                    dt   = min( dt*finc, dtmax )
                    beta = beta * fbeta
                end if

                if ( power < 0.d0 ) then
                    cumn = 0
                    dt   = dt * fdec
                    beta = beta0
                    va   = 0.d0
                end if

                temp = maxval( abs(fa) )
                if ( temp < fmax ) exit

                if ( step == stepmax ) then
                    write(*,*) "subroutine reached step maximum and existed with no force balance"
                    stop
                end if

                if( step == 1 ) write(*,*) 'step', '    dt         ', '              fmax   '

                if( mod(step,100) == 0 ) then
                    write(*,'(i6,2e25.15)') step, dt, temp
                end if

            end do

        end associate
    end subroutine

    ! main 2
    subroutine mini_fire_cp( tcon, tnet, opboxp_set, opxyzp_set, opxp_set, opyp_set, opzp_set, opstress_set, oppin_flag )
        implicit none

        ! para list
        type(con_t),     intent(inout)        :: tcon
        type(network_t), intent(in), optional :: tnet
        real(8),         intent(in), optional :: opboxp_set    ! target press
        real(8),         intent(in), optional :: opxyzp_set
        real(8),         intent(in), optional :: opxp_set, opyp_set, opzp_set
        real(8),         intent(in), optional :: opstress_set
        logical,         intent(in), optional :: oppin_flag

        ! local
        real(8) :: dstrain
        logical :: nonnetwork_flag
        real(8) :: boxp_set,  xyzp_set(free)
        logical :: boxp_flag, xyzp_flag(free)
        logical :: cp_flag
        real(8) :: stress_set
        logical :: cs_flag
        logical :: pin_flag

        real(8) :: lainv(free)
        real(8) :: dt, beta, temp
        real(8) :: onembeta, betavndfn
        integer :: cumn
        integer :: i, ipin

        ! 1. network or not
        nonnetwork_flag = .true.
        if ( present( tnet ) ) nonnetwork_flag = .false.

        ! 2 stress
        cs_flag = .false.
        if ( present(opstress_set) ) then
            stress_set = opstress_set
            cs_flag    = .true.
        end if

        ! 3.0 model of constant pressure
        ! 3.1 constant box press, change box xy length simutaneously
        boxp_flag = .false.
        if ( present( opboxp_set ) ) then
            if ( present(opxp_set) .or. present(opyp_set) .or. present(opzp_set) .or. &
                &present(opxyzp_set) ) then
                write(*,*) "Error set in cp"
                stop
            end if
            boxp_set  = opboxp_set
            boxp_flag = .true.
        end if
        ! 3.2 constant box press, change box xy length seperately.
        xyzp_flag = .false.
        if ( present(opxyzp_set) ) then
            if ( present(opxp_set) .or. present(opyp_set) .or. present(opzp_set) ) then
                write(*,*) "Error set in xyzp_set"
                stop
            end if
            xyzp_set  = opxyzp_set
            xyzp_flag = .true.
        else
            if( present( opxp_set ) ) then
                xyzp_set(1)  = opxp_set
                xyzp_flag(1) = .true.
            end if
            if ( present( opyp_set ) ) then
                xyzp_set(2)  = opyp_set
                xyzp_flag(2) = .true.
            end if
            if ( free == 3 .and. present( opzp_set ) ) then
                xyzp_set(free)  = opzp_set
                xyzp_flag(free) = .true.
            end if
        end if

        ! cp or not
        cp_flag = .false.
        if ( boxp_flag .or. any(xyzp_flag) ) cp_flag = .true.

        ! pin
        pin_flag = .false.
        if( present(oppin_flag) ) pin_flag = oppin_flag

        associate(                     &
            natom    => tcon%natom,    &
            ra       => tcon%ra,       &
            va       => tcon%va,       &
            fa       => tcon%fa,       &
            pinflag  => tcon%pinflag,  &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            strain   => tcon%strain,   &
            strainv  => tcon%strainv,  &
            strainf  => tcon%strainf,  &
            stress   => tcon%stress,   &
            la       => tcon%la,       &
            lav      => tcon%lav,      &
            laf      => tcon%laf       &
            )

            lainv = 1.d0 / la

            ! initial sets
            fa      = 0.d0 ; va      = 0.d0
            lav     = 0.d0 ; laf     = 0.d0
            strainv = 0.d0 ; strainf = 0.d0
            dt      = dt0  ; beta    = beta0
            cumn = 0

            ! pre
            if ( nonnetwork_flag ) then
                call make_list( nbfire, tcon )
                call calc_force_h( tcon, nbfire )
            else
                call calc_force_spring( tcon, tnet )
            end if

            !! pin
            if ( pin_flag ) then
                do ipin=1, natom
                    if ( pinflag(ipin) == 1 ) fa(:,ipin) = 0.d0
                end do
            end if
            !! end pin

            if ( cp_flag ) then
                if ( boxp_flag ) laf = press - boxp_set
                do i=1, free
                    if ( xyzp_flag(i) ) laf(i) = pressxyz(i) - xyzp_set(i)
                end do
            end if
            if ( cs_flag ) strainf = stress_set - stress

            do step=1, stepmax

                ! step length
                dt2  = dt    * 0.5d0
                dt22 = dt**2 * 0.5d0

                ! velocity verlet method / move 1 / config
                ra = ra + va * dt + fa * dt22
                va = va + fa * dt2

                ! velocity verlet method / move 1 / box
                if ( cp_flag ) then
                    !v affine deformation
                    if ( boxp_flag ) then
                        la(1)  = la(1)  + lav(1) * dt + laf(1) * dt22
                        lav(1) = lav(1) + laf(1) * dt2
                        !
                        la(2:free) = la(2:free) * ( la(1)*lainv(1) )
                        ra = ra * ( la(1)*lainv(1) )
                    else
                        do i=1, free
                            if ( xyzp_flag(i) ) then
                                la(i)   = la(i)  + lav(i) * dt + laf(i) * dt22
                                lav(i)  = lav(i) + laf(i) * dt2
                                ra(i,:) = ra(i,:) * ( lainv(i)*la(i) )
                            end if
                        end do
                    end if
                    !^
                    lainv = 1.d0 / la
                end if

                ! velocity verlet method / move 1 / shear
                if ( cs_flag ) then
                    dstrain = strainv * dt + strainf * dt22
                    strain  = strain  + dstrain
                    strainv = strainv + strainf * dt2

                    ra(1,:) = ra(1,:) + dstrain * ra(2,:)
                end if

                ! check list
                if ( nonnetwork_flag .and. check_list( nbfire, tcon ) ) then
                    call make_list( nbfire, tcon )
                end if

                ! velocity verlet method / force
                if ( nonnetwork_flag ) then
                    call calc_force_h( tcon, nbfire )
                else
                    call calc_force_spring( tcon, tnet )
                end if

                !! pin
                if ( pin_flag ) then
                    do ipin=1, natom
                        if ( pinflag(ipin) == 1 ) fa(:,ipin) = 0.d0
                    end do
                end if
                !! end pin

                if ( cp_flag ) then
                    if ( boxp_flag ) laf = press - boxp_set
                    do i=1, free
                        if ( xyzp_flag(i) ) laf(i) = pressxyz(i) - xyzp_set(i)
                    end do
                end if
                if ( cs_flag ) strainf = stress_set - stress

                ! velocity verlet method / move 2
                va = va + fa * dt2
                if ( cp_flag ) lav = lav + laf * dt2
                if ( cs_flag ) strainv = strainv + strainf * dt2

                ! fire
                cumn  = cumn + 1
                power = sum( fa * va ) + sum( laf*lav ) + strainv * strainf

                fn = sqrt( sum( fa**2 ) + sum(laf**2) + strainf**2 )
                vn = sqrt( sum( va**2 ) + sum(lav**2) + strainv**2 )

                onembeta  = 1.d0 - beta
                betavndfn = beta*vn/fn

                va      = onembeta * va      + betavndfn * fa
                lav     = onembeta * lav     + betavndfn * laf
                strainv = onembeta * strainv + betavndfn * strainf

                if ( power > 0.d0 .and. cumn > nmin ) then
                    dt   = min( dt*finc, dtmax )
                    beta = beta * fbeta
                end if

                if ( power < 0.d0 ) then
                    cumn    = 0
                    dt      = dt * fdec
                    beta    = beta0
                    va      = 0.d0
                    lav     = 0.d0
                    strainv = 0.d0
                end if

                temp = maxval( abs(fa) )
                if ( cp_flag ) then
                    if ( boxp_flag ) temp = max( temp, abs( press - boxp_set  ) )
                    do i=1, free
                        if ( xyzp_flag(i) ) temp = max( temp, abs( pressxyz(i) - xyzp_set(i) ) )
                    end do
                end if
                if ( cs_flag ) temp = max( temp, abs( stress - stress_set ) )

                if ( temp < fmax ) then
                    !write(*,'(5es16.6)') 1.0, tcon%press, tcon%pressx, tcon%pressy, tcon%stress
                    exit
                end if

                if ( step == stepmax ) then
                    write(*,*) "subroutine reached step maximum and existed with no force balance"
                    stop
                end if

                !if( step == 1 ) write(*,*) 'step', '    dt         ', '              fmax   '

                !if( mod(step,100) == 0 ) then
                !    if ( free == 2 ) then
                !        write(*,'(i6,5e16.6)') step, dt, temp, press, pressxyz(:)
                !    elseif ( free == 3 ) then
                !        write(*,'(i6,6e16.6)') step, dt, temp, press, pressxyz(:)
                !    end if
                !end if

            end do

            ! write(*,*) tcon%Ea

        end associate
    end subroutine

end module
