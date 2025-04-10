program main
    use mo_integral
    implicit none

    real(8) :: integration, comparison
    real(8) :: low_x = 4.d0
    real(8) :: up_x  = 9.d0
    integer :: nbins = 1000

    abstract_func_h => func_sqrtx

    ! appoint the integrand function as sqrt(), one could also give the name of
    ! the function directly: integration = integral( func_sqrtx, low_x, up_x ).
    ! Adjusting bin: integration = integral( func_sqrtx, low_x, up_x, nbins).
    integration = integral( abstract_func_h, low_x, up_x )

    comparison = 38.d0 / 3.d0

    print*, integration, comparison

contains

    pure function func_sqrtx( tx ) result( tfunc )
        implicit none

        ! para list
        real(8), intent(in) :: tx
        real(8)             :: tfunc

        tfunc = dsqrt( tx )
    end function

end program
