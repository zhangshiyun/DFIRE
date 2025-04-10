module mo_integral
    implicit none

    integer, parameter, private :: default_bins = 10000

    abstract interface
            pure function abstract_func( tx ) result( tfunc ) ! the integrand function
            real(8), intent(in) :: tx
            real(8)             :: tfunc
        end function
    end interface

    procedure(abstract_func), pointer :: abstract_func_h => null()

contains

    function integral( abstract_func, tlow_x, tup_x, tnbins ) result( integration )
        implicit none

        ! para list
        real(8), external    :: abstract_func
        real(8), intent(in)  :: tlow_x, tup_x
        integer, intent(in), optional :: tnbins

        ! result
        real(8) :: integration

        ! local
        integer :: nbins
        real(8) :: dx, xmid, ymid
        integer :: i

        if( present(tnbins) ) then
            nbins = tnbins
        else
            nbins = default_bins
        end if

        integration = 0.d0
        dx   = ( tup_x - tlow_x ) / dble(nbins)
        xmid = tlow_x - dx*0.5d0

        do i=1, nbins
            xmid = xmid + dx
            ymid = abstract_func(xmid)
            integration = integration + ymid*dx
        end do
    end function

end module
