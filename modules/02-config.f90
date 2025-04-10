module mo_config
    !
    !  base struct and method of configuration generating, storing...
    !
    use mo_syst
    use mo_math, only: init_rand
    implicit none

    type con_t
        ! Atom
        !! number of particles
        integer :: natom
        !! configuration, velocity, force
        real(8), allocatable, dimension(:,:) :: ra, va, fa
        !! radius
        real(8), allocatable, dimension(:)   :: r, iiea
        !! v for polydisperse system
        real(8), allocatable, dimension(:)   :: radius_dispersity
        !! v for pin particle
        integer, allocatable, dimension(:)   :: pinflag

        ! box
        !! volume fraction
        real(8) :: phi
        real(8) :: dstrain
        !! lengthes of box, la = (lx, ly, [lz])
        real(8) :: la(free)
        !! v for vary lengthes of box
        real(8) :: lav(free)
        real(8) :: laf(free)
        !! strain, for shear
        real(8) :: strain
        real(8) :: strainv, strainf

        ! sets
        real(8) :: T
        ! properties
        !! energy
        real(8) :: Ea, Ek, Ev
        !! stress tensor
        real(8) :: stress, press
        !!! pressures of x, y and z
        real(8) :: pressxyz(free)
    contains
        ! con%dra(i,j) => (xij, yij, zij)
        procedure :: dra         => calc_dra
        ! con%vec(i, [xk,yk,zk]) => round( xj-xk, yj-yk, zj-zk )
        procedure :: vec         => calc_vec
        ! con%dra(i,j) => ||rij||
        procedure :: len         => calc_len
        ! con%calc_phi() => phi
        procedure :: calc_phi    => calc_phi
        ! compress and shear
        procedure :: concompress => concompress
        procedure :: conshear    => conshear
    end type

    type(con_t), target :: con, con0, contemp, contemp1, contemp2, con_g, con1, con2
    type(con_t), target :: bcon, bcon0, bcontemp, bcontemp1, bcontemp2, bcon_g

contains

    subroutine init_system( tcon, tnatom, tphi )
        !
        ! initiate system, allocate memory of system
        !
        implicit none

        ! para list
        type(con_t),       intent(inout) :: tcon
        integer,           intent(in)    :: tnatom
        real(8), optional, intent(in)    :: tphi

        tcon%natom = tnatom
        if ( present( tphi ) ) tcon%phi = tphi

        allocate( tcon%ra(free,tnatom), tcon%r(tnatom), tcon%va(free,tnatom), tcon%fa(free,tnatom) )
        allocate( tcon%iiea(tnatom))
    end subroutine

    subroutine add_radius_dispersity( tcon, deta, opseed )
        implicit none

        ! para list
        type(con_t),       intent(inout) :: tcon
        real(8),           intent(in)    :: deta
        integer, optional, intent(in)    :: opseed

        if ( present(opseed) ) then
            call init_rand(opseed)
            if ( .not. allocated(tcon%radius_dispersity) ) then
                allocate( tcon%radius_dispersity(tcon%natom) )
            end if
            call random_number(tcon%radius_dispersity)
            tcon%radius_dispersity = tcon%radius_dispersity - 0.5d0
        else
            if ( .not. allocated(tcon%radius_dispersity) ) then
                write(*,*) "Error, You should give seed at least one time"
                stop
            end if
        end if

        tcon%r = tcon%r + tcon%radius_dispersity * deta
    end subroutine

    subroutine gen_rand_config( tcon, tseed, tphi )
        !
        !  generate random configuration with given seed and phi
        !
        implicit none

        ! para list
        type(con_t), intent(inout)        :: tcon
        integer,     intent(inout)        :: tseed
        real(8),     intent(in), optional :: tphi

        ! local
        integer :: i, j

        ! initialized rand
        call init_rand(tseed)
        tseed = 0

        if ( present( tphi ) ) tcon%phi = tphi

        associate(                &
            natom  => tcon%natom, &
            ra     => tcon%ra,    &
            fa     => tcon%fa,    &
            va     => tcon%va,    &
            r      => tcon%r,     &
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            ! radius
            r(1:natom/2)       = 0.5d0
            r(natom/2+1:natom) = 0.5d0 * ratio

            ! box length
            la     = calc_box_length(tcon)
            strain = 0.d0

            ! config
            call random_number(ra)
            do i=1, natom
                do j=1, free
                    ra(j,i) = ( ra(j,i) - 0.5d0 ) * la(j)
                end do
            end do

            ! f v
            va = 0.d0
            fa = 0.d0

        end associate
    end subroutine

    subroutine gen_lattice_triangle( tcon, tphi )
        !
        !  generate 2D triangle lattice
        !  phi_c = pi / 2 / sqrt(3) .= 0.9068996821171089
        !
        implicit none

        ! para list
        type(con_t), intent(inout)        :: tcon
        real(8),     intent(in), optional :: tphi

        ! local
        integer :: nxy, i, ii, jj
        real(8) :: a, b
        real(8), parameter :: xyoffset = 1.d-2

        if ( present( tphi ) ) tcon%phi = tphi

        associate(                &
            natom  => tcon%natom, &
            ra     => tcon%ra,    &
            r      => tcon%r,     &
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            ! box length
            la     = calc_box_trianglelattice( tcon )
            strain = 0.d0

            ! cell numbers and unit
            nxy = nint( sqrt( dble(natom) ) )
            a   = la(1) / nxy
            b   = a * sqrt(3.d0) / 2.d0

            ! radius
            r = 0.5d0

            ! config
            i=0
            do ii=1, nxy
                do jj=1, nxy
                    i = i + 1
                    ra(1,i) = a * mod( (i-1), nxy ) + a/2.d0 * ( (i-1)/nxy ) + xyoffset
                    ra(2,i) = b * ( (i-1)/nxy ) + xyoffset
                end do
            end do

            do i=1, natom
                ra(:,i) = ra(:,i) - anint( ra(:,i) / la ) * la
            end do

        end associate
    end subroutine

    subroutine gen_lattice_sc( tcon, tphi )
        !
        !  generate 2d simple-cube lattice
        !
        implicit none

        ! para list
        type(con_t), intent(inout)        :: tcon
        real(8),     intent(in), optional :: tphi

        ! local
        integer :: nxy, i, ii, jj
        real(8) :: volume, sdisk
        real(8) :: a, b
        real(8), parameter :: xyzoffset = 1.d-2

        if ( present( tphi ) ) tcon%phi = tphi

        associate(                &
            natom  => tcon%natom, &
            ra     => tcon%ra,    &
            r      => tcon%r,     &
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            ! sc close packing
            ! box: 1
            ! sphere: pi/4
            ! phi: pi/4 = 0.7853981633974483 <= 0.78540

            ! radius
            r = 0.5d0

            ! box length
            sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(r**free)
            volume = sdisk / tcon%phi
            la = volume**(1.d0/2.d0)
            strain = 0.d0

            ! cell numbers and unit
            nxy = nint( sqrt( dble(natom) ) )
            a   = la(1) / nxy

            ! config
            i=0
            do ii=1, nxy
                do jj=1, nxy
                    i = i + 1
                    ra(1,i) = (ii-1) * a
                    ra(2,i) = (jj-1) * a
                end do
            end do

            do i=1, natom
                ra(:,i) = ra(:,i) + xyzoffset
            end do

            do i=1, natom
                ra(:,i) = ra(:,i) - anint( ra(:,i) / la ) * la
            end do

        end associate
    end subroutine

    subroutine gen_lattice_fcc( tcon, tphi )
        !
        !  3D fcc
        !
        implicit none

        ! para list
        type(con_t), intent(inout)        :: tcon
        real(8),     intent(in), optional :: tphi

        ! local
        integer :: nxyz, i, ii, jj, kk
        real(8) :: a, sdisk, volume
        real(8), parameter :: xyzoffset = 1.d-2

        if ( present( tphi ) ) tcon%phi = tphi

        associate(                &
            natom  => tcon%natom, &
            ra     => tcon%ra,    &
            r      => tcon%r,     &
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            ! fcc close packing
            ! box: 2 x^2 = 4 => x^3 = 2^1.5
            ! sphere: 4*4/3*pi*1/8 = 2pi/3
            ! phi: pi/3/2^0.5 = 0.7404804896930609 <= 0.7404805

            ! radius
            r = 0.5d0

            ! box length
            sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(r**free)
            volume = sdisk / tcon%phi
            la = volume**(1.d0/3.d0)
            strain = 0.d0

            ! primitive cell
            nxyz = nint( dble(natom/4)**(1.d0/3.d0) )
            if ( 4*nxyz**3 /= natom ) then
                stop "wrong natom"
            end if
            a = la(1) / nxyz

            ! config
            i = 0
            do ii=1, nxyz
                do jj=1, nxyz
                    do kk=1, nxyz
                        i = i + 1
                        ra(1,i) = (ii-1) * a
                        ra(2,i) = (jj-1) * a
                        ra(3,i) = (kk-1) * a
                        i = i + 1
                        ra(1,i) = (ii-1+0.5d0) * a
                        ra(2,i) = (jj-1+0.5d0) * a
                        ra(3,i) = (kk-1      ) * a
                        i = i + 1
                        ra(1,i) = (ii-1+0.5d0) * a
                        ra(2,i) = (jj-1      ) * a
                        ra(3,i) = (kk-1+0.5d0) * a
                        i = i + 1
                        ra(1,i) = (ii-1      ) * a
                        ra(2,i) = (jj-1+0.5d0) * a
                        ra(3,i) = (kk-1+0.5d0) * a
                    end do
                end do
            end do

            do i=1, natom
                ra(:,i) = ra(:,i) + xyzoffset
            end do
        end associate

        call trim_config( tcon )
    end subroutine

    subroutine gen_lattice_bcc( tcon, tphi )
        !
        !  3D bcc
        !
        implicit none

        ! para list
        type(con_t), intent(inout)        :: tcon
        real(8),     intent(in), optional :: tphi

        ! local
        integer :: nxyz, i, ii, jj, kk
        real(8) :: a, sdisk, volume
        real(8), parameter :: xyzoffset = 1.d-2

        if ( present( tphi ) ) tcon%phi = tphi

        associate(                &
            natom  => tcon%natom, &
            ra     => tcon%ra,    &
            r      => tcon%r,     &
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            ! bcc close packing
            ! box: 3*x^2 = 4 => x^3 = (4/3)^1.5
            ! shperes: 2*4/3*pi/8 = pi / 3
            ! phi: pi / 3 / (4/3)^1.5 = 0.6801747615878317 <= 0.680175

            ! radius
            r = 0.5d0

            ! box length
            sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(r**free)
            volume = sdisk / tcon%phi
            la = volume**(1.d0/3.d0)
            strain = 0.d0

            ! primitive cell
            nxyz = nint( dble(natom/2)**(1.d0/3.d0) )
            if ( 2*nxyz**3 /= natom ) then
                write(*,*) "wrong natom"; stop
            end if
            a = la(1) / nxyz

            ! config
            i = 0
            do ii=1, nxyz
                do jj=1, nxyz
                    do kk=1, nxyz
                        i = i + 1
                        ra(1,i) = (ii-1) * a
                        ra(2,i) = (jj-1) * a
                        ra(3,i) = (kk-1) * a
                        i = i + 1
                        ra(1,i) = (ii-1+0.5) * a
                        ra(2,i) = (jj-1+0.5) * a
                        ra(3,i) = (kk-1+0.5) * a
                    end do
                end do
            end do

            do i=1, natom
                ra(:,i) = ra(:,i) + xyzoffset
            end do
        end associate

        call trim_config( tcon )
    end subroutine

    subroutine read_config( tcon, tfilename, tnatom, tphi )
        !
        !  read configuration from text file
        !
        implicit none

        ! para list
        type(con_t),  intent(inout)        :: tcon
        character(*), intent(in)           :: tfilename
        integer,      intent(in)           :: tnatom
        real(8),      intent(in), optional :: tphi

        ! local
        integer :: i
        real(8) :: factor

        ! allocate array of tcon
        !if ( present( tphi ) ) then
        !    call init_system( tcon, tnatom, tphi )
        !else
        !    call init_system( tcon, tnatom )
        !end if

        ! read config
        open(901,file=tfilename)
            read(901, *) tcon%la, tcon%strain
            do i=1, tnatom
                read(901,*) tcon%ra(:,i), tcon%r(i)
            end do
        close(901)

        !!!!
        !tcon%r = tcon%r * 0.5

        ! phi
        tcon%phi = pi * sum(tcon%r**2) / product(tcon%la)
    end subroutine

    pure function calc_box_length(tcon) result(l)
        implicit none

        ! para list
        type(con_t), intent(in) :: tcon

        ! result
        real(8) :: l

        ! local
        real(8) :: phi
        real(8) :: sdisk, volume

        phi = tcon%phi

        ! V_n(r) = pi^(n/2) / Gamma( n/2 + 1 ) * r^n
        sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(tcon%r**free)

        ! box length
        volume = sdisk / phi
        l      = volume ** ( 1.d0 / dble(free) )
    end function

    function calc_box_trianglelattice( tcon ) result( tla )
        implicit none

        ! para list
        type(con_t), intent(in)  :: tcon

        ! result
        real(8), dimension(free) :: tla

        ! local
        integer :: nxy

        associate(               &
            natom => tcon%natom, &
            phi   => tcon%phi    &
            )

            nxy = nint( sqrt( dble(natom) ) )
            if ( nxy**2 /= natom .and. mod(nxy,2) == 0 ) then
                stop "wrong natom"
            end if

            tla(1) = sqrt( natom * pi / sqrt(12.d0) / phi )
            tla(2) = tla(1) * sqrt(3.d0) / 2.d0

        end associate
    end function

    pure function calc_dra( this, ti, tj ) result(dra)
        !
        !  dra = rj - ri = rij = [xij, yij, zij]
        !
        implicit none

        ! para list
        class(con_t), intent(in) :: this
        integer,      intent(in) :: ti
        integer,      intent(in) :: tj

        ! result
        real(8), dimension(free) :: dra

        ! local
        real(8) :: rai(free), raj(free)
        integer :: k, cory, iround(free)

        associate(                &
            ra     => this%ra,    &
            la     => this%la,    &
            strain => this%strain &
            )

            rai = ra(:,ti)
            raj = ra(:,tj)

            dra = raj - rai

            cory = nint( dra(free) / la(free) )
            dra(1) = dra(1) - strain * la(free) * cory

            do k=1, free-1
                iround(k) = nint( dra(k) / la(k) )
            end do
            iround(free) = cory

            do k=1, free
                dra(k) = dra(k) - iround(k) * la(k)
            end do

        end associate
    end function

    pure function calc_vec( this, ti, traj ) result(dra)
        !
        !  calculate vector distance of particle i and point[x,y,z]
        !
        implicit none

        ! para list
        class(con_t), intent(in) :: this
        integer,      intent(in) :: ti
        real(8), dimension(free), intent(in) :: traj

        ! result
        real(8), dimension(free) :: dra

        ! local
        real(8) :: rai(free)
        integer :: k, cory, iround(free)

        associate(                &
            ra     => this%ra,    &
            la     => this%la,    &
            strain => this%strain &
            )

            rai = ra(:,ti)

            dra = traj - rai

            cory = nint( dra(free) / la(free) )
            dra(1) = dra(1) - strain * la(free) * cory

            do k=1, free-1
                iround(k) = nint( dra(k) / la(k) )
            end do
            iround(free) = cory

            do k=1, free
                dra(k) = dra(k) - iround(k) * la(k)
            end do

        end associate
    end function

    pure function find_closest( this, traj, n ) result(list)
        !
        !  select nth closest particles to point[x,y,z]
        !
        use mo_math, only: swap
        implicit none

        ! para list
        class(con_t), intent(in) :: this
        real(8), dimension(free), intent(in) :: traj
        integer, intent(in) :: n

        ! result
        integer :: list(n)

        ! local
        real(8) :: closest(n), dis
        integer :: i, j, itemp
        real(8) :: rtemp

        closest = 1.d8
        list = 0

        do i=1, this%natom
            dis = norm2(this%vec(i, traj))
            if ( dis < closest(n) ) then
                closest(n) = dis
                list(n) = i
                do j=n, 2, -1
                    if ( closest(j) < closest(j-1) ) then
                        itemp = list(j)
                        list(j) = list(j-1)
                        list(j-1) = itemp
                        rtemp = closest(j)
                        closest(j) = closest(j-1)
                        closest(j-1) = rtemp
                    else
                        exit
                    end if
                end do
            end if
        end do
    end function

    pure function calc_len( this, ti, tj ) result(tl)
        !
        !  l = ||dra||
        !
        implicit none

        ! para list
        class(con_t), intent(in) :: this
        integer,      intent(in) :: ti, tj

        ! result
        real(8) :: tl

        ! local
        real(8) :: dra(free)

        dra = calc_dra( this, ti, tj )
        tl = norm2(dra)
    end function

    subroutine trim_config( tcon, opsumxyz )
        !
        !  make sure all the partiles locate in the periodical cell
        !
        implicit none

        ! para list
        type(con_t), intent(inout)        :: tcon
        logical,     intent(in), optional :: opsumxyz  ! set center of mass to zero

        ! local
        real(8) :: lainv(free), temp
        integer :: iround(free), cory
        integer :: i, k

        associate(                &
            natom  => tcon%natom, &
            ra     => tcon%ra,    &
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            lainv = 1.d0 / la

            do i=1, natom

                cory = nint( ra(free,i) * lainv(free) )
                ra(1,i) = ra(1,i) - strain * la(free) * cory

                do k=1, free-1
                    iround(k) = nint( ra(k,i) * lainv(k) )
                end do
                iround(free) = cory

                do k=1, free
                    ra(k,i) = ra(k,i) - iround(k) * la(k)
                end do

            end do

            if ( present( opsumxyz ) ) then
                if ( opsumxyz .eqv. .true. ) then
                    do i=1, free
                        temp = sum(ra(i,:)) / natom
                        ra(i,:) = ra(i,:) - temp
                    end do
                end if
            end if

        end associate
    end subroutine

    pure function calc_phi(this) result(re)
        !
        !  calculate volume fraction
        !  phi = volume of particles / volume of box
        !
        implicit none

        ! para list
        class(con_t), intent(in) :: this

        ! result
        real(8) :: re

        ! local
        real(8) :: sdisk, volume

        volume = product( this%la(1:free) )

        !if ( free == 2 ) then
        !    sdisk = pi * sum( this%r(1:this%natom)**2 )
        !elseif ( free == 3 ) then
        !    sdisk = 4.d0/3.d0 * pi * sum( this%r(1:this%natom)**3 )
        !end if
        sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(this%r**free)

        re = sdisk / volume
    end function

    subroutine gen_pin( tcon, tn )
        implicit none

        ! para list
        type(con_t), intent(inout) :: tcon
        integer,     intent(in)    :: tn

        ! local
        integer :: i, testi
        integer :: flag
        real(8) :: temp, dij

        if ( allocated(tcon%pinflag) .and. size(tcon%pinflag) /= tcon%natom ) then
            deallocate( tcon%pinflag )
        end if

        if ( .not. allocated( tcon%pinflag ) ) then
            allocate( tcon%pinflag( tcon%natom ) )
        end if

        associate(                   &
            pinflag => tcon%pinflag, &
            r       => tcon%r,       &
            natom   => tcon%natom    &
            )

            pinflag = 0
            pinflag(1) = 1

            testi = 1
            do while ( sum(pinflag) < tn )

                testi = testi + 1

                flag = 0
                do i=1, testi
                    if ( pinflag(i) == 1 ) then
                        temp = calc_len( tcon, i, testi )
                        dij = r(testi) + r(i)
                        if ( temp < 3.d0*dij ) flag = flag + 1
                    end if
                end do

                if ( flag == 0 ) then
                    pinflag( testi ) = 1
                end if

            end do

        end associate
    end subroutine

    subroutine concompress( this, xyz, de, opaffine)
        implicit none

        ! para list
        class(con_t), intent(inout) :: this
        integer, intent(in) :: xyz
        real(8), intent(in) :: de
        logical, optional   :: opaffine

        ! local
        logical :: affine_flag

        this%la(xyz) = this%la(xyz) * ( 1.d0 - de )

        affine_flag = .true.
        if ( present(opaffine) .and. ( opaffine .eqv. .false. ) ) affine_flag = .false.

        if ( affine_flag ) then
            this%ra(xyz,:) = this%ra(xyz,:) * ( 1.d0 - de )
        end if
    end subroutine

    subroutine conshear( this, de, opy, opx, opaffine)
        implicit none

        ! para list
        class(con_t), intent(inout) :: this
        real(8), intent(in) :: de
        logical, optional   :: opaffine
        integer, optional   :: opx, opy
        integer :: x, y

        ! local
        logical :: affine_flag

        x = 1; y = free
        if ( present(opx) ) x = opx
        if ( present(opy) ) y = opy

        this%strain = this%strain + de

        affine_flag = .true.
        if ( present(opaffine) .and. ( opaffine .eqv. .false. ) ) affine_flag = .false.

        if ( affine_flag ) then
            this%ra(x,:) = this%ra(x,:) + de * this%ra(y,:)
        end if
    end subroutine

end module

! ToDo : use hdf5 insteadly
subroutine save_config_to( tcon, tfilename )
    use mo_config
    implicit none

    ! para list
    type(con_t),  intent(in) :: tcon
    character(*), intent(in) :: tfilename

    ! local
    integer :: i

    associate(                  &
        natom   => tcon%natom,  &
        ra      => tcon%ra,     &
        r       => tcon%r,      &
        la      => tcon%la,     &
        strain  => tcon%strain  &
        )

        open(901,file=tfilename)
            write(901,'(3es26.16)') dble(natom), tcon%phi, 0.d0
            write(901,'(3es26.16)') la, strain
            do i=1, natom
                if ( free == 3 ) then
                    write(901,'(4es26.16)') ra(:,i), r(i)
                elseif( free == 2 ) then
                    write(901,'(3es26.16)') ra(:,i), r(i)
                end if
            end do
        close(901)

    end associate
end subroutine

subroutine save_config_debug( tcon, tfilename )
    use mo_config
    implicit none

    ! para list
    type(con_t),  intent(in) :: tcon
    character(*), intent(in) :: tfilename

    ! local
    integer :: i

    associate(                &
        natom  => tcon%natom, &
        ra     => tcon%ra,    &
        va     => tcon%va,    &
        fa     => tcon%fa,    &
        r      => tcon%r,     &
        la     => tcon%la,    &
        strain => tcon%strain &
        )

        open(901,file=tfilename)
            write(901,'(3es26.16)') dble(natom), tcon%phi, 0.d0
            write(901,'(3es26.16)') la, strain
            do i=1, natom
                write(901,'(7es26.16)') ra(:,i), r(i), va(:,i), fa(:,i)
            end do
        close(901)

    end associate
end subroutine

! save config with wall
subroutine save_config_copy( tcon, tfilename )
    use mo_config
    implicit none

    ! para list
    type(con_t),  intent(in) :: tcon
    character(*), intent(in) :: tfilename

    ! local
    integer :: i, ii, jj
    real(8), dimension(free) :: xyoffset

    associate(                &
        natom  => tcon%natom, &
        ra     => tcon%ra,    &
        r      => tcon%r,     &
        la     => tcon%la,    &
        strain => tcon%strain &
        )

        open(901,file=tfilename,status="new")

            write(901,'(3es26.16)') dble(natom), tcon%phi, 0.d0
            write(901,'(3es26.16)') 3*la(1:free), strain

            do ii=-1, 1
                do jj=-1, 1

                    xyoffset(1) = ii * la(1)
                    xyoffset(2) = jj * la(2)
                    do i=1, natom
                        write(901,'(3es26.16)') ra(:,i)+xyoffset, r(i)
                    end do

                end do
            end do

        close(901)

    end associate
end subroutine
