module  mo_molecule_npt
    use mo_syst
    use mo_config
    use mo_list
    use mo_molecule
    use mo_extra_molecule
    implicit none
!Assumption: 1) Integration method: prediction-correction; 2) Ensemble simulation: constraint method; 3) Soft particles; 4) Using neighbor list 5) The length unit is the averaged particle diameter
!procedure(abstract_extra_force), pointer :: tp_abstr_force => null() !the force type in this module shall be abstract_extra_force

    type, extends(Base_molecule), public :: tpmolecule_npt
        real(8) :: gear0 = 3.d0 / 8.d0
        real(8) :: gear2 = 3.d0 / 4.d0
        real(8) :: gear3 = 1.d0 / 6.d0

        real(8) :: dt, temper, pre, resc_prob
        real(8) :: c1, c2, c3
        real(8) :: coeff0, coeff2, coeff3

        real(8) :: lagr_t, lagr_p
        real(8) :: unit_v, nlcut0
        real(8) :: la1(free), la2(free), la3(free), la_list(free)
        real(8), allocatable, dimension(:,:) :: ra1, ra2, ra3
        real(8), allocatable, dimension(:,:) :: va1, va2, va3
    contains
        procedure, pass  :: initial  =>  initial_npt
        procedure, pass  :: body     =>  body_npt
        procedure, pass  :: clean    =>  clean_npt

        procedure, pass :: predic    =>  predic_npt
        procedure, pass :: correc    =>  correc_npt
        procedure, pass :: constri   =>  constri_npt
        procedure, pass :: rescale_t =>  rescale_t_npt
        procedure, pass :: rescale_p =>  rescale_p_npt
    end type tpmolecule_npt

    !private ::  predic_npt
    !private ::  correc_npt
    !private ::  constri_npt
    !private ::  rescale_t_npt
    !private ::  rescale_p_npt

contains

    subroutine initial_npt( this, tcon, tnb, tp_abstr_force, tset_molecule )
        implicit none

        ! para list
        class(tpmolecule_npt), intent(inout) :: this
        type(con_t),           intent(inout) :: tcon
        type(list_t),          intent(inout) :: tnb
        class(tpset_molecule), intent(in)    :: tset_molecule
        logical,               external      :: tp_abstr_force

        !local
        integer :: i
        logical :: get_force

        select type(tset_molecule)
        type is (tpset_npt)

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

            lagr_t   =>  this%lagr_t, &
            lagr_p   =>  this%lagr_p, &

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

        !do i = 1, natom
        !    call one_make_list (tnb, tcon, i)
        !end do
        !call full_make_list( tnb, tcon)
        call make_list( tnb, tcon)
        la_list = la
        this%nlcut0 = tnb%nlcut

        unit_v = sum( (2.d0*radius)**dble(free) ) / dble(natom)

        call random_number(va)
        call rescale_t_npt( this, tcon )
        call rescale_p_npt( this, tcon, tnb, tp_abstr_force )

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


    subroutine body_npt( this, tcon, tnb, tp_abstr_force, tout_molecule )
        implicit none

         ! para list
        class(tpmolecule_npt), intent(inout) :: this
        type(con_t),           intent(inout) :: tcon
        type(list_t),          intent(inout) :: tnb
        logical,               external      :: tp_abstr_force
        type(tpout_molecule),  intent(inout), optional :: tout_molecule

        ! local
        logical :: get_force
        real(8) :: prob, temp

        call predic_npt( this, tcon )

        temp = dble(free) * (tcon%la(1) - this%la_list(1))
        if(temp < 0) then
            tnb%nlcut = this%nlcut0 + temp
        endif

        if ( tnb%nlcut < 1.d-2 .or. check_list( tnb, tcon )) then
            call make_list( tnb, tcon )
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
            call rescale_p_npt( this, tcon, tnb, tp_abstr_force )
        end if

    end subroutine


    subroutine correc_npt( this, tcon )
        implicit none

        ! para list
        class(tpmolecule_npt), intent(inout) :: this
        type(con_t),           intent(inout) :: tcon

        !local
        integer :: i
        real(8) :: corr_la(free), corr_ra(free), corr_va(free), temp(free)

        associate(                   &
            lagr_t  =>  this%lagr_t, &
            lagr_p  =>  this%lagr_p, &
            coeff0  =>  this%coeff0, &
            coeff2  =>  this%coeff2, &
            coeff3  =>  this%coeff3, &

            natom   =>  tcon%natom,  &
            la      =>  tcon%la,     &
            ra      =>  tcon%ra,     &
            va      =>  tcon%va,     &
            fa      =>  tcon%fa,     &
            c1      =>  this%c1,     &
            c2      =>  this%c2,     &
            c3      =>  this%c3,     &
            ra1     =>  this%ra1,    &
            ra2     =>  this%ra2,    &
            ra3     =>  this%ra3,    &
            va1     =>  this%va1,    &
            va2     =>  this%va2,    &
            va3     =>  this%va3,    &
            la1     =>  this%la1,    &
            la2     =>  this%la2,    &
            la3     =>  this%la3     &
        )

        do i = 1, natom
            temp     = ra1(:,i)
            ra1(:,i) = va (:,i) + lagr_p * ra(:,i)
            corr_ra  = ra1(:,i) - temp
            ra (:,i) = ra (:,i) + coeff0 * corr_ra
            ra2(:,i) = ra2(:,i) + coeff2 * corr_ra
            ra3(:,i) = ra3(:,i) + coeff3 * corr_ra

            temp     = va1(:,i)
            va1(:,i) = fa (:,i) - lagr_t * va(:,i)
            corr_va  = va1(:,i) - temp
            va (:,i) = va (:,i) + coeff0 * corr_va
            va2(:,i) = va2(:,i) + coeff2 * corr_va
            va3(:,i) = va3(:,i) + coeff3 * corr_va
        end do

        temp    = la1
        la1     = lagr_p * la
        corr_la = la1 - temp
        la      = la  + coeff0 * corr_la
        la2     = la2 + coeff2 * corr_la
        la3     = la3 + coeff3 * corr_la

        end associate
    end subroutine


    subroutine predic_npt( this, tcon )
        implicit none

        ! para list
        class(tpmolecule_npt), intent(inout) :: this
        type(con_t),           intent(inout) :: tcon

        associate(                    &
            la      =>  tcon%la,      &
            ra      =>  tcon%ra,      &
            va      =>  tcon%va,      &
            c1      =>  this%c1,      &
            c2      =>  this%c2,      &
            c3      =>  this%c3,      &
            ra1     =>  this%ra1,     &
            ra2     =>  this%ra2,     &
            ra3     =>  this%ra3,     &
            va1     =>  this%va1,     &
            va2     =>  this%va2,     &
            va3     =>  this%va3,     &
            la1     =>  this%la1,     &
            la2     =>  this%la2,     &
            la3     =>  this%la3      &
        )

        ra  = ra  + c1 * ra1 + c2 * ra2 + c3 * ra3
        ra1 = ra1 + c1 * ra2 + c2 * ra3
        ra2 = ra2 + c1 * ra3
        va  = va  + c1 * va1 + c2 * va2 + c3 * va3
        va1 = va1 + c1 * va2 + c2 * va3
        va2 = va2 + c1 * va3
        la  = la  + c1 * la1 + c2 * la2 + c3 * la3
        la1 = la1 + c1 * la2 + c2 * la3
        la2 = la2 + c1 * la3

        end associate
    end subroutine



    subroutine constri_npt( this, tcon )
         implicit none

        ! para list
        class(tpmolecule_npt), intent(inout) :: this
        type(con_t),           intent(inout) :: tcon

        associate(                     &
            lagr_t   => this%lagr_t,   &
            lagr_p   => this%lagr_p,   &
            natom    => tcon%natom,    &
            va       => tcon%va,       &
            fa       => tcon%fa,       &
            Ek       => tcon%Ek,       &
            vili     => ex_var%vili,   &
            chixi    => ex_var%chixi,  &
            forchi   => ex_var%forchi  &
            )

        Ek = sum( va**2.d0 )
        lagr_t = sum( va * fa) / Ek
        lagr_p = - forchi / ( dble(free**2) * ( vili + Ek/dble(free) ) + chixi )
        Ek = 0.5 * Ek  / dble(natom)

        end associate
    end subroutine



    subroutine rescale_p_npt( this, tcon, tnb, tp_abstr_force )
         implicit none

        ! para list
        class(tpmolecule_npt), intent(inout) :: this
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

        call make_list( tnb, tcon )
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
                call make_list( tnb, tcon )
                this%la_list = tcon%la
            end if

            get_force = tp_abstr_force( tcon, tnb, 0 )
            pvar = ( ttemper + vili ) * unit_v / product( la )

        end do

        end associate
    end subroutine



 subroutine rescale_t_npt( this, tcon )
         implicit none

        ! para list
        class(tpmolecule_npt), intent(inout) :: this
        type(con_t),           intent(inout) :: tcon

        !local
        integer :: i
        real(8) :: tvar

        associate(                     &
            temper   => this%temper,   &
            natom    => tcon%natom,    &
            va       => tcon%va        &
            )

         do i=1, free
            tvar    = sum( va(i,:) ) / natom
            va(i,:) = va(i,:) - tvar
         end do

        tvar = sum( va**2 ) / ( natom * free )
        va   = sqrt( temper / tvar ) * va

        end associate
    end subroutine


    subroutine clean_npt( this )
        implicit none

        ! para list
        class(tpmolecule_npt),  intent(inout) :: this

        if ( allocated(this%ra1) )     deallocate( this%ra1 )
        if ( allocated(this%ra2) )     deallocate( this%ra2 )
        if ( allocated(this%ra3) )     deallocate( this%ra3 )
        if ( allocated(this%va1) )     deallocate( this%va1 )
        if ( allocated(this%va2) )     deallocate( this%va2 )
        if ( allocated(this%va3) )     deallocate( this%va3 )

    end subroutine


end module


