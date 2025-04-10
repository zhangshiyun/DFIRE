module mo_corr
    implicit none

    type corr_t
        integer :: nbins
        real(8) :: wbin
        real(8) :: vmin, vmax
        real(8), allocatable, dimension(:) :: cum_n
        real(8), allocatable, dimension(:) :: cum_corr
    contains
        procedure :: x    => calc_xi
        procedure :: corr => calc_corr
    end type

    type corrab_t
        real(8), allocatable, dimension(:) :: a, b
        real(8) :: corrab
    contains
        procedure :: corr => calc_corrab
    end type

    type(corr_t)   :: kvcorr
    type(corrab_t) :: corrab

contains

    subroutine init_corr( tcorr, nbins, vmin, vmax )
        implicit none

        ! para list
        type(corr_t) :: tcorr
        integer :: nbins
        real(8) :: vmin, vmax

        tcorr%nbins = nbins
        tcorr%vmin  = vmin
        tcorr%vmax  = vmax
        tcorr%wbin  = (vmax-vmin) / nbins
        allocate( tcorr%cum_corr(nbins), tcorr%cum_n(nbins) );
    end subroutine

    pure function calc_xi( this, i ) result(re)
        implicit none

        ! para list
        class(corr_t), intent(in) :: this
        integer,       intent(in) :: i

        ! result
        real(8) :: re

        re = this%vmin + ( i - 0.5 ) * this%wbin
    end function

    pure function calc_corr( this, i ) result(re)
        implicit none

        ! para list
        class(corr_t), intent(in) :: this
        integer,       intent(in) :: i

        ! result
        real(8) :: re

        re = 0
        if ( this%cum_n(i) /= 0 ) then
            re = this%cum_corr(i) / this%cum_n(i)
        end if
    end function

    pure function calc_corrab( this ) result( re )
        use mo_math, only: corr
        implicit none

        ! para list
        class(corrab_t), intent(in) :: this

        ! result
        real(8) :: re

        re = corr(this%a, this%b)
    end function

end module
