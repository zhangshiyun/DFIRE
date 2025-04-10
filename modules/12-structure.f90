module mo_structure
    use mo_syst
    use mo_config
    use mo_list
    use mo_math, only: sortperm, unitv
    implicit none

    type sk_t
        integer :: nmax, natom
        real(8), allocatable, dimension(:,:) :: ks  ! save kx, ky, k
        real(8), allocatable, dimension(:)   :: sks ! save sk
        real(8), dimension(free) :: kbase, la
        integer, dimension(free) :: nxy
        logical :: init_flag = .false.
        logical :: cumflag
        integer :: cumcount  ! average counter
    end type

contains

    subroutine init_sk(con, sk)
        implicit none

        ! para list
        type(con_t), intent(in)    :: con   ! need natom and length of box
        type(sk_t),  intent(inout) :: sk

        ! local
        real(8), allocatable, dimension(:,:) :: kxy1, kxy2  ! save candidate ks
        integer, allocatable, dimension(:)   :: perm        ! order of kxy
        integer :: i, j, icum1, icum2

        if ( sk%init_flag ) return

        ! zero counter
        sk%cumcount = 0

        sk%natom = con%natom
        sk%la    = con%la

        ! kbase = 2pi / Lxy
        sk%kbase = 2.d0 * pi / sk%la
        sk%nxy   = 2*nint(sk%la)

        allocate( kxy1(free+1, product(sk%nxy)) )
        allocate( kxy2(free+1, product(sk%nxy)) )

        ! for isotropic sys
        if ( all( sk%la == sk%la(1) ) ) then
            icum1 = 0
            do i=0, sk%nxy(1)
                do j=i, sk%nxy(2)
                    if ( i==0 .and. j==0 ) cycle
                    icum1 = icum1 + 1
                    kxy1(1, icum1) = i*sk%kbase(1)
                    kxy1(2, icum1) = j*sk%kbase(2)
                    kxy1(3, icum1) = sqrt( (sk%kbase(1)*i)**2 + (sk%kbase(2)*j)**2 )
                end do
            end do
        else
            print*, 51
        end if

        ! select ks
        ! sort ks, get a order list
        perm = sortperm(icum1, kxy1(free+1,1:icum1))
        ! save first k
        icum2 = 1
        kxy2(:, icum2) = kxy1(:, perm(1))
        ! filter redundant k, save to kxy2
        do i=2, icum1
            if ( kxy2(free+1, icum2) /= kxy1(free+1,perm(i)) ) then
                icum2 = icum2 + 1
                kxy2(:, icum2) = kxy1(:, perm(i))
            end if
        end do

        ! save k, allocate sk array
        sk%ks  = kxy2

        ! decide length of sk
        sk%nmax = icum2
        allocate(sk%sks(icum2))
        sk%sks = 0.d0

        ! deallocate temp array
        deallocate(kxy1, kxy2, perm)
    end subroutine

    subroutine calc_sk(con, sk)
        implicit none

        ! para list
        type(con_t), intent(in)    :: con
        type(sk_t),  intent(inout) :: sk

        ! local
        real(8) :: k(free)
        integer :: i, j
        real(8) :: kr, rhiok_c, rhiok_s

        sk%cumcount = sk%cumcount + 1

        do i=1, sk%nmax
            k = sk%ks(1:free, i)
            rhiok_c = 0.d0
            rhiok_s = 0.d0
            do j=1, con%natom
                kr = sum( k * con%ra(:,j))
                rhiok_c = rhiok_c + cos(kr)
                rhiok_s = rhiok_s + sin(kr)
            end do
            sk%sks(i) = sk%sks(i) + rhiok_s**2 + rhiok_c**2
        end do
    end subroutine

    subroutine endofsk(sk)
        implicit none

        ! para list
        type(sk_t), intent(inout) :: sk

        sk%sks = sk%sks / sk%cumcount / sk%natom
    end subroutine

    function calc_psi(tcon, n) result(psi)
        implicit none

        ! para list
        type(con_t), intent(in) :: tcon
        integer,     intent(in) :: n

        ! result
        complex(16), dimension(tcon%natom) :: psi

        ! lcoal
        type(voro_t) :: voroni
        real(8)      :: nij(free)
        complex(16)  :: cpxtemp
        integer      :: i, j, k

        call init_voro( voroni, tcon )
        call calc_voro( voroni )

        psi(:) = cmplx(0.0, 0.0)

        do i=1, tcon%natom
            do k=1, voroni%list(i)%nbsum
                j       = voroni%list(i)%nblist(k)
                nij     = unitv( tcon%dra(i,j) )
                cpxtemp = cmplx(nij(1),nij(2))
                psi(i)  = psi(i) + cpxtemp**n
            end do
        end do
        psi(:)  = psi(:) / voroni%list(:)%nbsum
    end function

end module
