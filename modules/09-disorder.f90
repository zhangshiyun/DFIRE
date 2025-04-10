module mo_disorder
    use mo_syst
    use mo_config
    implicit none

    real(8) :: eta, eta_spring
    real(8), allocatable, dimension(:,:) :: radisorder

contains

    subroutine init_disorder( tradisorder, tcon, tseed )
        use mo_math, only: randuvec
        implicit none

        ! para list
        type(con_t) :: tcon
        real(8), allocatable, dimension(:,:) :: tradisorder
        integer :: tseed

        ! local
        real(8) :: temp
        integer :: i

        !temp = rand( tseed )
        call init_rand( tseed )
        tseed = 0

        associate( natom => tcon%natom )

            if( .not. allocated(radisorder) ) allocate( radisorder(free,natom) )

            do i=1, natom
                radisorder(:,i) = randuvec(free)
            end do

        end associate
    end subroutine

    subroutine make_disorder( tcon, tradisorder, teta )
        implicit none

        ! para list
        type(con_t) :: tcon
        real(8), allocatable, dimension(:,:) :: tradisorder
        real(8) :: teta

        associate( ra => tcon%ra )

            ra = ra + teta * tradisorder

        end associate
    end subroutine

end module
