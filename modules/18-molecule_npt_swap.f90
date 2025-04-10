module  mo_molecule_npt_swap
    use mo_syst
    use mo_config
    use mo_list
    use mo_molecule
    use mo_extra_molecule
    use mo_molecule_npt
    use mo_math, only: swap
    implicit none
!Assumption: 1) Integration method: prediction-correction; 2) Ensemble simulation: constraint method; 3) Soft particles; 4) Using neighbor list 5) The length unit is the smallest particle diameter
!procedure(abstract_extra_force), pointer :: tp_abstr_force => null() !the force type in this module shall be abstract_extra_force

    type, extends(tpmolecule_npt), public :: tpmolecule_npt_swap
        real(8)  :: swap_prob
    contains
        !override
        procedure, pass  :: initial    =>  initial_npt_swap
        procedure, pass  :: body       =>  body_npt_swap
        procedure, pass  :: rescale_p  =>  rescale_p_npt_swap
        !extend
        procedure, pass  :: sweep_swap =>  sweep_swap
        procedure, pass  :: pair_swap  =>  pair_swap

        !inherit  ::  variables, clean_npt, predic_npt, correc_npt, constri_npt, rescale_t_npt
    end type tpmolecule_npt_swap



contains

    subroutine initial_npt_swap( this, tcon, tnb, tp_abstr_force, tset_molecule )
        implicit none

        ! para list
        class(tpmolecule_npt_swap), intent(inout) :: this
        type(con_t),                intent(inout) :: tcon
        type(list_t),               intent(inout) :: tnb
        class(tpset_molecule),      intent(in)    :: tset_molecule
        logical,                    external      :: tp_abstr_force

        !local
        integer :: i
        logical :: get_force

        select type(tset_molecule)
        type is (tpset_npt_swap)

        associate(                    &
            unit_v   =>  this%unit_v, &
            temper   =>  this%temper, &
            pre      =>  this%pre,    &
            dt       =>  this%dt,     &

            c1       =>  this%c1,     &
            c2       =>  this%c2,     &
            c3       =>  this%c3,     &
            coeff0   =>  this%coeff0, &
            coeff2   =>  this%coeff2, &
            coeff3   =>  this%coeff3, &

            la_list  =>  this%la_list,&
            la1      =>  this%la1,    &
            la2      =>  this%la2,    &
            la3      =>  this%la3,    &

            lagr_t   => this%lagr_t,  &
            lagr_p   => this%lagr_p,  &

            natom    =>  tcon%natom,  &
            radius   =>  tcon%r,      &
            ra       =>  tcon%ra,     &
            va       =>  tcon%va,     &
            fa       =>  tcon%fa,     &
            la       =>  tcon%la      &
            )

        if ( .not. allocated(this%ra1) )   allocate( this%ra1(free,natom) )
        if ( .not. allocated(this%ra2) )   allocate( this%ra2(free,natom) )
        if ( .not. allocated(this%ra3) )   allocate( this%ra3(free,natom) )
        if ( .not. allocated(this%va1) )   allocate( this%va1(free,natom) )
        if ( .not. allocated(this%va2) )   allocate( this%va2(free,natom) )
        if ( .not. allocated(this%va3) )   allocate( this%va3(free,natom) )

        this%swap_prob = tset_molecule%swap_prob
        this%resc_prob = tset_molecule%resc_prob

        temper = tset_molecule%temper
        pre    = tset_molecule%pre
        dt     = tset_molecule%dt

        c1 = dt
        c2 = c1 * dt / 2.d0
        c3 = c2 * dt / 3.d0
        coeff0 = this%gear0 * c1
        coeff2 = this%gear2 * c1 / c2
        coeff3 = this%gear3 * c1 / c3

        !unit_v = sum( (2.d0*radius)**dble(free) ) / dble(natom)

        unit_v = 1.d0

        call full_make_list( tnb, tcon)
        la_list = la
        this%nlcut0 = tnb%nlcut

        call random_number(va)
        call rescale_t_npt( this, tcon )
        call rescale_p_npt_swap( this, tcon, tnb, tp_abstr_force )

        get_force = tp_abstr_force( tcon, tnb, 0 )

        call constri_npt( this, tcon )

        associate(                    &
            ra1      =>  this%ra1,    &
            ra2      =>  this%ra2,    &
            ra3      =>  this%ra3,    &
            va1      =>  this%va1,    &
            va2      =>  this%va2,    &
            va3      =>  this%va3     &
            )

        ra1 = va + lagr_p * ra
        ra2 = 0.d0
        ra3 = 0.d0

        va1 = fa - lagr_t * va
        va2 = 0.d0
        va3 = 0.d0

        la1 = lagr_p * la
        la2 = 0.d0
        la3 = 0.d0

        end associate
        end associate

        end select
    end subroutine


    subroutine body_npt_swap( this, tcon, tnb, tp_abstr_force, tout_molecule )
        implicit none

         ! para list
        class(tpmolecule_npt_swap), intent(inout) :: this
        type(con_t),                intent(inout) :: tcon
        type(list_t),               intent(inout) :: tnb
        logical,                    external      :: tp_abstr_force
        type(tpout_molecule),       intent(inout), optional :: tout_molecule

        ! local
        integer :: i
        logical :: get_force
        real(8) :: prob, temp

        call predic_npt( this, tcon )

        call sweep_swap( this, tcon, tnb, tp_abstr_force )

        temp = dble(free) * (tcon%la(1) - this%la_list(1))
        if(temp < 0) then
            tnb%nlcut = this%nlcut0 + temp
        endif

        if ( tnb%nlcut < 1.d-2 .or. check_list( tnb, tcon )) then
            call full_make_list( tnb, tcon )
            this%la_list = tcon%la
        end if

        get_force = tp_abstr_force( tcon, tnb, 0 )

        call constri_npt( this, tcon )

        call correc_npt( this, tcon )

        if ( present( tout_molecule ) ) then

            associate(                    &
                unit_v   => this%unit_v,  &
                natom    => tcon%natom,   &
                radius   => tcon%r,       &
                la       => tcon%la,      &
                vili     => ex_var%vili   &
            )

            tout_molecule%phi   = sqrt( pi**free ) / gamma( dble(free)/2.d0+1.d0 ) * sum( radius**free ) / product( la )
            tout_molecule%T     = 2.d0 * tcon%Ek / dble(free)
            tout_molecule%press = ( tout_molecule%T * dble(natom) + vili ) * unit_v / product( la )
            tout_molecule%Ek    = tcon%Ek
            tout_molecule%Ea    = tcon%Ea
            tout_molecule%Ev    = tcon%Ek + tcon%Ea

            end associate
        end if

        call random_number( prob )
        if( prob < this%resc_prob ) then
            call rescale_t_npt( this, tcon )
            call rescale_p_npt_swap( this, tcon, tnb, tp_abstr_force )
        end if

    end subroutine


    subroutine sweep_swap( this, tcon, tnb, tp_abstr_force )
        implicit none

        ! para list
        class(tpmolecule_npt_swap), intent(inout) :: this
        type(con_t),                intent(inout) :: tcon
        type(list_t),               intent(inout) :: tnb
        logical,                    external      :: tp_abstr_force

        !local
        integer :: n, ti, tj, tn
        real(8) :: var, prob

        associate(                      &
            natom     => tcon%natom,    &
            radius    => tcon%r,        &
            swap_prob => this%swap_prob &
            )

        do n = 1, natom

            call random_number(prob)

            if(prob < swap_prob) then

                call random_number(var)
                ti = floor( var * dble(natom) ) + 1
                call random_number(var)
                tj = floor( var * dble(natom) ) + 1

                if( ti .eq. tj ) cycle

                if( pair_swap( this, tcon, tnb, tp_abstr_force, ti, tj ) ) then
                    if( radius(ti) > radius(tj) ) then
                        tn = ti
                    else
                        tn = tj
                    end if
                    call one_make_list( tnb, tcon, tn )
                end if

            end if

        enddo

        end associate
    end subroutine


    logical function pair_swap( this, tcon, tnb, tp_abstr_force, ti, tj )
        implicit none

        ! para list
        class(tpmolecule_npt_swap), intent(inout) :: this
        type(con_t),                intent(inout) :: tcon
        type(list_t),               intent(inout) :: tnb
        logical,                    external      :: tp_abstr_force
        integer,                    intent(in)    :: ti, tj

        !local
        real(8) :: prob, temp
        real(8) :: ea0, ea1
        logical :: get_potential

        associate(                  &
            temper => this%temper,  &
            radius => tcon%r,       &
            dea    => ex_var%dea    &
            )

        pair_swap = .false.

        temp = abs( radius(ti) - radius(tj) )
        if(temp > tnb%nlcut) then
            go to 100
        end if

        ea0 = 0.d0
        get_potential = tp_abstr_force( tcon, tnb, ti )
        ea0 = ea0 + dea
        get_potential = tp_abstr_force( tcon, tnb, tj )
        ea0 = ea0 + dea

        call swap( radius(ti), radius(tj) )

        ea1 = 0.d0
        get_potential = tp_abstr_force( tcon, tnb, ti )
        ea1 = ea1 + dea
        get_potential = tp_abstr_force( tcon, tnb, tj )
        ea1 = ea1 + dea

        if( ea1 < ea0 ) then
            pair_swap = .true.
            go to 100
        end if

        prob = exp( (ea0 - ea1)/temper )
        call random_number(temp)
        if( temp < prob ) then
            pair_swap = .true.
            go to 100
        else
            call swap( radius(ti), radius(tj) )
        end if

100     end associate

        return
    end function


    subroutine rescale_p_npt_swap( this, tcon, tnb, tp_abstr_force )
         implicit none

        ! para list
        class(tpmolecule_npt_swap), intent(inout) :: this
        type(con_t),           intent(inout) :: tcon
        type(list_t),          intent(inout) :: tnb
        logical,               external      :: tp_abstr_force

        !local
        real(8) :: pvar, ttemper, temp, factor
        logical :: get_force

        associate(                     &
            unit_v   => this%unit_v,   &
            pre      => this%pre,      &
            la_list  => this%la_list,  &
            natom    => tcon%natom,    &
            ra       => tcon%ra,       &
            la       => tcon%la,       &
            chixi    => ex_var%chixi,  &
            vili     => ex_var%vili    &
            )

        ttemper = sum( tcon%va**2.d0 ) / dble(free)

        call full_make_list( tnb, tcon )
        la_list = la
        get_force = tp_abstr_force( tcon, tnb, 0 )
        pvar = ( ttemper + vili ) * unit_v / product( la )

        do while( abs( (pre-pvar)/pre ) > 1.d-2)

            factor = 1.d0 + (pvar - pre) / (pvar + chixi/dble(free**2)/product(la) ) / dble(free)

            la = la * factor
            ra = ra * factor

            temp = dble(free) * (tcon%la(1) - this%la_list(1))
            if(temp < 0) then
                tnb%nlcut = this%nlcut0 + temp
            endif
            if ( tnb%nlcut < 1.d-2 ) then
                call full_make_list( tnb, tcon )
                this%la_list = tcon%la
            end if

            get_force = tp_abstr_force( tcon, tnb, 0 )
            pvar = ( ttemper + vili ) * unit_v / product( la )

        end do

        end associate
    end subroutine

end module


