module mo_list
    !
    !  verlet list, see more at page 147 of <computer simulation of liquids>
    !
    use mo_syst
    use mo_config
    implicit none

    ! global constants.
    !! skin
    real(8), private, parameter :: set_nlcut = 3.0d0
    !! for sake of memory saving, we consider listmax neighbor of one particle at most
    !! enlarge this if you study 3D system or high volume fraction system
    integer, private, parameter :: listmax = 80

    ! neighbor list of one particle
    type listone_t
        ! neighbor number
        integer    :: nbsum
        ! neighbor index
        integer    :: nblist(listmax)
        ! precalculated relative relation, used for distance calculation in perodical cell
        integer(1) :: iround(free,listmax)
        integer(1) :: cory(listmax)
        !
        real(4)    :: con0(free)
    end type

    ! neighbor list struct of system
    type list_t
        integer                                    :: natom
        type(listone_t), allocatable, dimension(:) :: list
        ! contact number
        integer, allocatable, dimension(:)         :: nbi
        ! tag particle is rattler or not
        integer, allocatable, dimension(:)         :: rattlerflag
        real(8)    :: nlcut
    end type

    type(list_t) :: nb, nb1, nb2


    ! voronoi cell
    type voro_one_t
        integer :: nbsum
        integer :: nblist(listmax)
        integer :: vid1(listmax), vid2(listmax)
    end type

    type voro_t
        integer :: natom
        type(voro_one_t), allocatable, dimension(:) :: list
        real(8), allocatable, dimension(:,:) :: center, vertex
        real(8) :: la(free), strain
    contains
        procedure :: init      => init_voro
        procedure :: decompose => calc_voro
    end type

    type(voro_t) :: voro

contains

    subroutine init_list(tnb, tcon)
        !
        !  allocate memory
        !
        implicit none

        ! para list
        type(list_t), intent(inout) :: tnb
        type(con_t),  intent(in)    :: tcon

        ! local
        integer :: tnatom

        tnatom    = tcon%natom
        tnb%natom = tnatom
        tnb%nlcut = set_nlcut

        if ( allocated(tnb%list) ) then
            if ( size(tnb%list) /= tnatom ) then
                deallocate( tnb%list )
                allocate( tnb%list(tnatom) )
            end if
        else
            allocate( tnb%list(tnatom) )
        end if
    end subroutine

    subroutine make_list(tnb, tcon)
        !
        !  make list of system
        !
        implicit none

        ! para list
        type(con_t),  intent(in)    :: tcon
        type(list_t), intent(inout) :: tnb

        ! local
        real(8) :: lainv(free), dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, itemp

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            r      => tcon%r,      &
            la     => tcon%la,     &
            strain => tcon%strain, &
            list   => tnb%list,    &
            nlcut  => tnb%nlcut    &
            )

            nlcut = set_nlcut
            lainv = 1.d0 / la

            ! set nbsum to zero
            list(:)%nbsum = 0

            do i=1, natom

                list(i)%con0 = ra(:,i)
                rai          = ra(:,i)
                ri           = r(i)

                do j=i+1, natom

                    raj = ra(:,j)
                    rj  = r(j)

                    dra = raj - rai

                    cory   = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij  = ri + rj

                    if ( rij2 > ( nlcut )**2 ) cycle

                    if ( list(i)%nbsum < listmax ) then
                        itemp                   = list(i)%nbsum
                        itemp                   = itemp + 1
                        list(i)%nbsum           = itemp
                        list(i)%nblist(itemp)   = j
                        list(i)%iround(:,itemp) = iround
                        list(i)%cory(itemp)     = cory
                    end if

                end do
            end do

        end associate
    end subroutine

    subroutine make_bond(tcon, tnb, bonds)

        implicit none

        type(con_t), intent(in)     :: tcon
        type(list_t), intent(in)     :: tnb
        integer, dimension(tcon%natom,tcon%natom), intent(inout)      :: bonds


        ! local
        real(8) :: lainv(free), dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, itemp, jj

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            radius => tcon%r,      &
            la     => tcon%la,     &
            strain => tcon%strain, &
            list   => tnb%list,    &
            nlcut  => tnb%nlcut    &
            )

        do i=1, natom

            rai = ra(:,i)
            ri  = radius(i)

            do jj=1, list(i)%nbsum

                j = list(i)%nblist(jj)
                iround = list(i)%iround(:,jj)
                cory = list(i)%cory(jj)

                raj = ra(:,j)
                rj  = radius(j)

                dra = raj - rai
                dra(1) = dra(1) - cory * strain * la(free)

                do k=1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )
                dij = ri + rj

                if ( rij2 > dij**2 ) cycle

                bonds(i,j) = 1
                bonds(j,i) = 1

            end do
        end do
        end associate
    end subroutine

    subroutine one_make_list( tnb, tcon, tn )
        !
        !  make list for one particle
        !
        implicit none

        ! para list
        type(con_t),  intent(in)    :: tcon
        type(list_t), intent(inout) :: tnb
        integer,      intent(in)    :: tn

        ! local
        real(8) :: lainv(free), dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, itemp

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            r      => tcon%r,      &
            la     => tcon%la,     &
            strain => tcon%strain, &
            list   => tnb%list,    &
            nlcut  => tnb%nlcut    &
            )

            lainv = 1.d0 / la
            i = tn

            list(i)%nbsum = 0
            list(i)%con0  = ra(:,i)
            rai           = ra(:,i)
            ri            = r(i)

            do j = 1, natom

                raj = ra(:,j)
                rj  = r(j)

                dra = raj - rai
                cory   = nint( dra(free) * lainv(free) )
                dra(1) = dra(1) - strain * la(free) * cory

                do k = 1, free-1
                    iround(k) = nint( dra(k) * lainv(k) )
                end do
                iround(free) = cory

                do k = 1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )
                dij  = ri + rj

                if ( rij2 > ( dij+nlcut )**2 ) cycle

                if ( list(i)%nbsum < listmax ) then
                    itemp                   = list(i)%nbsum
                    itemp                   = itemp + 1
                    list(i)%nbsum           = itemp
                    list(i)%nblist(itemp)   = j
                    list(i)%iround(:,itemp) = iround
                    list(i)%cory(itemp)     = cory
                end if

                if ( list(j)%nbsum < listmax ) then
                    itemp = list(j)%nbsum
                    if ( all( list(j)%nblist(1:itemp) /= i ) ) then
                        itemp                   = itemp + 1
                        list(j)%nbsum           = itemp
                        list(j)%nblist(itemp)   = i
                        list(j)%iround(:,itemp) = - iround
                        list(j)%cory(itemp)     = - cory
                    end if
                end if

            end do

        end associate
    end subroutine

    subroutine full_make_list( tnb, tcon )
        !
        !  make list for all the particles
        !
        implicit none

        ! para list
        type(con_t),  intent(in)    :: tcon
        type(list_t), intent(inout) :: tnb

        ! local
        real(8) :: lainv(free), dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, itemp

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            r      => tcon%r,      &
            la     => tcon%la,     &
            strain => tcon%strain, &
            list   => tnb%list,    &
            nlcut  => tnb%nlcut    &
            )

            nlcut = set_nlcut
            lainv = 1.d0 / la

            ! set nbsum to zero
            list(:)%nbsum = 0

            do i=1, natom

                list(i)%con0 = ra(:,i)
                rai          = ra(:,i)
                ri           = r(i)

                do j=i+1, natom

                    raj = ra(:,j)
                    rj  = r(j)

                    dra = raj - rai

                    cory   = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij  = ri + rj

                    if ( rij2 > ( dij+nlcut )**2 ) cycle

                    if ( list(i)%nbsum < listmax ) then
                        itemp                   = list(i)%nbsum
                        itemp                   = itemp + 1
                        list(i)%nbsum           = itemp
                        list(i)%nblist(itemp)   = j
                        list(i)%iround(:,itemp) = iround
                        list(i)%cory(itemp)     = cory
                    end if
                    if ( list(j)%nbsum < listmax ) then
                        itemp                   = list(j)%nbsum
                        itemp                   = itemp + 1
                        list(j)%nbsum           = itemp
                        list(j)%nblist(itemp)   = i
                        list(j)%iround(:,itemp) = - iround
                        list(j)%cory(itemp)     = - cory
                    end if

                end do
            end do

        end associate
    end subroutine

    subroutine calc_z( tnb, tcon )
        !
        !  calculate coordination number
        !
        implicit none

        ! para list
        type(con_t),  intent(in)    :: tcon
        type(list_t), intent(inout) :: tnb

        ! local
        real(8) :: dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, jj

        if ( allocated( tnb%nbi ) .and. size(tnb%nbi) /= tcon%natom ) then
            deallocate( tnb%nbi )
        end if

        if ( .not. allocated( tnb%nbi ) ) then
            allocate( tnb%nbi(tcon%natom) )
        end if

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            r      => tcon%r,      &
            la     => tcon%la,     &
            strain => tcon%strain, &
            list   => tnb%list,    &
            nbi    => tnb%nbi      &
            )

            ! set nbsum to zero
            nbi = 0

            do i=1, natom

                rai          = ra(:,i)
                ri           = r(i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)

                    raj = ra(:,j)
                    rj  = r(j)

                    dra = raj - rai

                    cory = list(i)%cory(jj)
                    iround = list(i)%iround(:,jj)

                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij  = ri + rj

                    if ( rij2 > ( dij )**2 ) cycle

                    nbi(i) = nbi(i) + 1
                    nbi(j) = nbi(j) + 1

                end do
            end do

        end associate
    end subroutine

    function check_list( tnb, tcon ) result(flag)
        !
        !  determine remake list or not
        !
        implicit none

        ! para list
        type(list_t), intent(in) :: tnb
        type(con_t),  intent(in) :: tcon

        ! result
        logical :: flag

        ! local
        real(8) :: maxdis, dra(free), dr2
        integer :: i

        associate(               &
            natom => tcon%natom, &
            ra    => tcon%ra,    &
            nlcut => tnb%nlcut   &
            )

            maxdis = 0.d0
            do i=1, tcon%natom
                dra = tcon%ra(:,i) - tnb%list(i)%con0
                dr2 = sum( dra**2 )
                if ( maxdis < dr2 ) maxdis = dr2
            end do

        flag = .false.
        if ( maxdis > 0.25 * nlcut**2 ) flag = .true.

        end associate
    end function

    ! ToDo
    ! subroutine check_rattler
    function calc_rattler( tcon, tnblist ) result(flag)
        implicit none

        ! para list
        type(con_t),  intent(in)           :: tcon
        type(list_t), intent(in), optional :: tnblist

        ! result
        integer, dimension(tcon%natom) :: flag

        ! local
        type(list_t) :: lclist
        integer      :: i

        if ( present(tnblist) ) then
            lclist = tnblist
        else
            call init_list( lclist, tcon )
            call make_list( lclist, tcon )
        end if

        call calc_z( lclist, tcon )

        flag = 0
        do i=1, lclist%natom
            if ( lclist%nbi(i) == 0 ) then
                flag(i) = 1
            elseif ( lclist%nbi(i) <= (free-1) ) then
                write(*,*) "There exist unstable particle(s)"
                stop
            end if
        end do
    end function

    ! voro
    subroutine init_voro( this, tcon )
        implicit none

        ! para list
        class(voro_t), intent(inout) :: this
        type(con_t),  intent(in)     :: tcon

        this%natom  = tcon%natom
        this%center = tcon%ra
        this%la     = tcon%la
        this%strain = tcon%strain

        if( .not. allocated(this%list) ) then
            allocate( this%list(this%natom) )
        else if( size(this%list) /= tcon%natom ) then
            deallocate( this%list )
            allocate( this%list(tcon%natom) )
        end if
    end subroutine

    subroutine calc_voro( this )
        implicit none

        ! para list
        class(voro_t), intent(inout) :: this

        ! local
        integer, parameter :: maxcan = 200
        integer, dimension(maxcan) :: verts
        real(8), dimension(maxcan) :: px, py, ps
        integer, dimension(maxcan) :: tag
        real(8), parameter :: rcut = 3.5

        real(8) :: rai(free), raij(free), rijsq
        integer :: i, j, k, cory, ncan

        associate( natom  => this%natom,  &
                   list   => this%list,   &
                   con    => this%center, &
                   la     => this%la,     &
                   strain => this%strain  &
                   )

        list(:)%nbsum = 0

        ! main loop
        do i=1, natom
            rai = con(:,i)
            k   = 0
            do j=1, natom
                if ( i == j ) cycle

                raij = con(:,j) - rai
                cory = nint( raij(2) / la(2) )
                raij(1) = raij(1) - cory * strain * la(2)
                raij(:) = raij(:) - anint( raij(:) / la ) * la

                rijsq = sum(raij**2)

                if ( rijsq < rcut**2 ) then
                    if ( k == maxcan ) then
                        write(*,*) k
                        cycle
                    end if
                    k      = k + 1
                    px(k)  = raij(1)
                    py(k)  = raij(2)
                    ps(k)  = rijsq
                    tag(k) = j
                end if
            end do

            ncan = k

            call sortps( px, py, ps, tag, ncan )

            call calc_voro_single( px, py, ps, maxcan, ncan, verts )

            do k=1, ncan
                if ( verts(k) /= 0 ) then
                    list(i)%nbsum = list(i)%nbsum + 1
                    list(i)%nblist(list(i)%nbsum) = tag(k)
                end if
            end do

        end do

        end associate
    end subroutine

    subroutine sortps( px, py, ps, tag, ncan )
        use mo_math, only: swap
        implicit none

        ! para list
        integer :: ncan
        real(8), dimension(ncan) :: px, py, ps
        integer, dimension(ncan) :: tag

        ! local
        integer :: i, imin

        do i=1, ncan-1

            imin = minloc( ps(i:ncan), 1 ) + i - 1
            if ( i == imin ) cycle

            call swap( px(i),  px(imin)  )
            call swap( py(i),  py(imin)  )
            call swap( ps(i),  ps(imin)  )
            call swap( tag(i), tag(imin) )

        end do
    end subroutine

    subroutine calc_voro_single( px, py, ps, maxcan, ncan, verts )
        implicit none

        integer :: maxcan, ncan, nv, ne
        real(8), dimension(maxcan) :: px, py, ps
        integer, dimension(maxcan) :: verts
        integer, parameter :: maxv = 200
        real(8), dimension(maxv) :: vx, vy
        integer, dimension(maxv) :: iv, jv

        logical :: flag
        real(8) :: ai,bi,ci, aj,bj,cj, det, detinv
        real(8) :: vxij, vyij
        real(8), parameter :: tol=1.d-8

        integer :: i, j, l, v

        !-- check
        if ( ncan < 3 ) then
            write(*,'('' less than 3 points given to work '',i5)') ncan
            stop
        end if

        v = 0
        do i=1, ncan-1

            ai =  px(i)
            bi =  py(i)
            ci = -ps(i)

            do j=i+1, ncan

                aj =  px(j)
                bj =  py(j)
                cj = -ps(j)

                det = ai*bj - aj*bi

                if( abs(det) >= tol ) then

                    detinv = 1.d0 / det
                    vxij = ( bi * cj - bj * ci ) * detinv
                    vyij = ( aj * ci - ai * cj ) * detinv
                    flag = .true.
                    l  = 1
                    do while( flag .and. l < ncan )
                        if ( l /= i .and. l /= j ) then
                            flag = ( ( px(l) * vxij + py(l) * vyij ) .le. ps(l) )
                        end if
                        l = l + 1
                    end do

                    if ( flag ) then
                        v = v + 1
                        if ( v > maxv ) stop 'too many vertices'
                        iv(v)  = i
                        jv(v)  = j
                        vx(v) = 0.5 * vxij
                        vy(v) = 0.5 * vyij
                    end if
                end if
            end do
        end do

        nv = v                      ! total vertex found
        if ( nv < 3 ) then
            write(*,'('' less than 3 vertices found in work '',i5)') nv
            stop
        end if

        verts(1:ncan) = 0

        do i = 1, nv
            verts(iv(i)) = verts(iv(i)) + 1
            verts(jv(i)) = verts(jv(i)) + 1
        end do

        flag = .true.
        ne   = 0

        do i = 1, ncan
            if ( verts(i) > 0 ) then
                ne = ne + 1
                if ( verts(i) /= 2 ) then    ! it is supposed that every neighbor contribute 2 vertex
                    flag = .false.
                end if
            end if
        end do

        if ( flag .eqv. .false. ) then
            write (*,'('' **** vertex error: degeneracy ? **** '')')
        end if

        if ( ne /= nv ) then
            write(*,'('' **** edge   error: degeneracy ? ****  '')')
        end if
    end subroutine

end module
