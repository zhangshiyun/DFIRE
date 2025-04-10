module mo_force
    use mo_syst
    use mo_config
    use mo_list
    use mo_network
    implicit none

    abstract interface
        subroutine abstract_force( tcon, tnb )
            import :: con_t, list_t
            type(con_t),  intent(inout)  :: tcon
            type(list_t), intent(in), optional :: tnb
        end subroutine
        logical function abstract_fun_force( tcon, tnb )
            import :: con_t, list_t
            type(con_t),  intent(inout)  :: tcon
            type(list_t), intent(in), optional :: tnb
        end function
    end interface

contains

    logical function calc_fun_force_lj( tcon, tnb )
        implicit none
        ! Warning:
        ! 1) pay attention to the range of list, make sure it is long enough to cover the interacted particles;
        ! 2) the texture is 50:50, which is different from the KA model!

        ! para list
        type(con_t),  intent(inout) :: tcon
        type(list_t), intent(in), optional :: tnb

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: lainv(free), ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        real(8) :: srinv, srinv3, srinv6, srinv12
        integer :: iround(free), cory
        integer :: i, j, k, jj

        ! lj
        !! V = econs * exx * [ (sxx/r)^12 - (sxx/r)^6 ]
        real(8), parameter :: econs = 4.d0 ! or 1.d0/72.d0
        real(8) :: exx, sxx                ! dummy coeffcient
        real(8), parameter :: eaa = 1.d0,  saa = 1.d0
        real(8), parameter :: ebb = 0.5d0, sbb = 0.88d0
        real(8), parameter :: eab = 1.5d0, sab = 0.8d0
        !! cut
        real(8), parameter :: rcut = 2.5d0
        real(8), parameter :: rcutinv = 1.d0/rcut
        !! smooth Vnew(r_cut) = 0
        !! smooth Fnew(r_cut) = 0
        !! Vnew(r) = V(r) - V(r_cut) - V'(r)*(r-r_cut)
        !! energy smooth V(r_cut) = Vcut * ess
        real(8), parameter :: Vcut = rcutinv**12 - rcutinv**6
        !! force  smooth // V'(r_cut) = - 1/sxx * [ 12/r_cut^13 - 6/r_cut^6 ] = Dvcut / sxx * ess
        real(8), parameter :: DVcut = - (12*rcutinv**13 - 6*rcutinv**7)


        associate(                    &
            natom    => tcon%natom,   &
            ra       => tcon%ra,      &
            fa       => tcon%fa,      &
            Ea       => tcon%Ea,      &
            la       => tcon%la,      &
            strain   => tcon%strain,  &
            stress   => tcon%stress,  &
            press    => tcon%press,   &
            pressxyz => tcon%pressxyz &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

        if ( present(tnb) ) then

            associate(list => tnb%list)

            do i=1, natom-1

                rai = ra(:,i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)

                    if (     i<=natom/2 .and. j<=natom/2 ) then
                        exx = eaa
                        sxx = saa
                    elseif ( i>natom/2  .and. j>natom/2  ) then
                        exx = ebb
                        sxx = sbb
                    else
                        exx = eab
                        sxx = sab
                    end if

                    raj = ra(:,j)
                    dra = raj - rai

                    iround = list(i)%iround(:,jj)
                    cory = list(i)%cory(jj)

                    dra = raj - rai
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )

                    if ( rij2 > (rcut*sxx)**2 ) cycle

                    rij = sqrt(rij2)

                    srinv   = sxx / rij
                    srinv3  = srinv  * srinv * srinv
                    srinv6  = srinv3 * srinv3
                    srinv12 = srinv6 * srinv6

                    Ea = Ea + exx * ( srinv12 - srinv6 - Vcut - (DVcut/sxx)*(rij-rcut*sxx) )

                    ! wij = fr * r = - V'(r) * r
                    wij  = exx * ( 12*srinv12 - 6*srinv6 + DVcut/sxx*rij )
                    wili = wili + wij

                    fr = wij / rij2

                    wilixyz = wilixyz + fr * dra(:)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            end associate

        else

            lainv = 1.d0 / la

            do i=1, natom-1

                rai = ra(:,i)

                do j=i+1, natom

                    if (     i<=natom/2 .and. j<=natom/2 ) then
                        exx = eaa
                        sxx = saa
                    elseif ( i>natom/2  .and. j>natom/2  ) then
                        exx = ebb
                        sxx = sbb
                    else
                        exx = eab
                        sxx = sab
                    end if

                    raj = ra(:,j)
                    dra = raj - rai

                    cory = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )

                    if ( rij2 > (rcut*sxx)**2 ) cycle

                    rij = sqrt(rij2)

                    srinv   = sxx / rij
                    srinv3  = srinv  * srinv * srinv
                    srinv6  = srinv3 * srinv3
                    srinv12 = srinv6 * srinv6

                    Ea = Ea + exx * ( srinv12 - srinv6 - Vcut - (DVcut/sxx)*(rij-rcut*sxx) )

                    ! wij = fr * r = - V'(r) * r
                    wij  = exx * ( 12*srinv12 - 6*srinv6 + DVcut/sxx*rij )
                    wili = wili + wij

                    fr = wij / rij2

                    wilixyz = wilixyz + fr * dra(:)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

        end if

        Ea       = econs * Ea
        fa       = econs * fa
        stress   = econs * stress  / product(la) / free
        press    = econs * wili    / product(la) / free
        pressxyz = econs * wilixyz / product(la)

        end associate

        calc_fun_force_lj = .true.

        return
    end function

    logical function calc_fun_force( tcon, tnb )
        implicit none

        ! para list
        type(con_t),  intent(inout) :: tcon
        type(list_t), intent(in), optional :: tnb

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        real(8) :: lainv(free)
        integer :: iround(free), cory
        integer :: i, j, k, jj

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz  &
            )

        Ea     = 0.d0
        fa     = 0.d0
        stress = 0.d0
        wili   = 0.d0; wilixyz = 0.d0

        if ( present(tnb) ) then

            associate(list => tnb%list)

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

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    wilixyz = wilixyz + fr * dra(:)**2

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            end associate

        else

            lainv = 1.d0 / la

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do j=i+1, natom

                    raj = ra(:,j)
                    rj  = radius(j)

                    dra = raj - rai

                    cory = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = ri + rj

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    wilixyz = wilixyz + fr * dra(:)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

        end if

        stress   = stress  / product(la) / free
        press    = wili    / product(la) / free
        pressxyz = wilixyz / product(la)

        end associate

        calc_fun_force = .true.

        return
    end function

    subroutine calc_force( tcon, tnb )
        implicit none

        ! para list
        type(con_t),  intent(inout) :: tcon
        type(list_t), intent(in), optional :: tnb

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        real(8) :: lainv(free)
        integer :: iround(free), cory
        integer :: i, j, k, jj

        associate(                    &
            natom    => tcon%natom,   &
            radius   => tcon%r,       &
            iiea     => tcon%iiea,       &
            ra       => tcon%ra,      &
            fa       => tcon%fa,      &
            Ea       => tcon%Ea,      &
            la       => tcon%la,      &
            strain   => tcon%strain,  &
            stress   => tcon%stress,  &
            press    => tcon%press,   &
            pressxyz => tcon%pressxyz &
            )

        Ea     = 0.d0
        fa     = 0.d0
        stress = 0.d0
        wili   = 0.d0; wilixyz = 0.d0
        iiea = 0.d0

        if ( present(tnb) ) then

            associate(list => tnb%list)

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

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    iiea(i) = 0.5d0 * ( 1.d0 - rij/dij )**alpha/alpha
                    iiea(j) = 0.5d0 * ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    wilixyz = wilixyz + fr * dra(:)**2

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            end associate

        else

            lainv = 1.d0 / la

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do j=i+1, natom

                    raj = ra(:,j)
                    rj  = radius(j)

                    dra = raj - rai

                    cory = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = ri + rj

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    wilixyz = wilixyz + fr * dra(:)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

        end if

        stress   = stress  / product(la) / free
        press    = wili    / product(la) / free
        pressxyz = wilixyz / product(la)

        end associate
    end subroutine

    subroutine calc_force_gel( tcon, tnb )
        implicit none

        ! para list
        type(con_t),  intent(inout) :: tcon
        type(list_t), intent(in), optional :: tnb

        ! local
        ! gel
        real(8), parameter :: rcut = 0.1d0
        real(8), parameter :: kout = 0.5d0
        real(8) :: rc

        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        real(8) :: dijrcut, dijrc
        integer :: iround(free), cory
        integer :: i, j, k, jj

        rc = rcut * kout / ( 1.d0 + kout )

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            list     => tnb%list       &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

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
                    dijrcut = dij * ( 1.d0 + rcut )
                    dijrc   = dij * ( 1.d0 + rc   )

                    if ( rij2 > dijrcut**2 ) cycle

                    rij = sqrt( rij2 )

                    if ( rij < dij ) then
                        Ea = Ea + ( 1.d0 - rij/dij )**2/2.d0
                        wij = (1.d0 - rij/dij) * rij / dij
                    elseif ( rij < dijrc ) then
                        Ea = Ea + kout*( 1.d0 - rij/dij )**2/2.d0
                        wij = kout*(1.d0 - rij/dij) * rij / dij
                    elseif ( rij < dijrcut ) then
                        wij = - kout * ( (dijrcut - rij)/dij ) * rij / dij
                    end if

                    wili = wili + wij

                    fr = wij / rij2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    wilixyz = wilixyz + fr * dra(:)**2

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            stress   = stress  / product(la) / free
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)
        end associate
    end subroutine

    subroutine calc_force_pin( tcon, tnb )
        implicit none

        ! para list
        type(con_t),  intent(inout) :: tcon
        type(list_t), intent(in), optional :: tnb

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        integer :: iround(free), cory
        integer :: i, j, k, jj

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            pinflag  => tcon%pinflag,  &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            list     => tnb%list       &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

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

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    if ( pinflag(i) == 0 .and. pinflag(j) == 0 ) then
                        wilixyz = wilixyz + fr * dra(:)**2

                        stress = stress - 2 * dra(1) * dra(2) * fr
                    end if

                end do

            end do

            stress   = stress  / product(la) / free
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)
        end associate
    end subroutine

    subroutine calc_force_lj( tcon, tnb )
        implicit none
        ! Warning:
        ! 1) pay attention to the range of list, make sure it is long enough to cover the interacted particles;
        ! 2) the texture is 50:50, which is different from the KA model!

        ! para list
        type(con_t),  intent(inout) :: tcon
        type(list_t), intent(in), optional :: tnb

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: lainv(free), ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        real(8) :: srinv, srinv3, srinv6, srinv12
        integer :: iround(free), cory
        integer :: i, j, k, jj

        ! lj
        !! V = econs * exx * [ (sxx/r)^12 - (sxx/r)^6 ]
        real(8), parameter :: econs = 4.d0 ! or 1.d0/72.d0
        real(8) :: exx, sxx                ! dummy coeffcient
        real(8), parameter :: eaa = 1.d0,  saa = 1.d0
        real(8), parameter :: ebb = 0.5d0, sbb = 0.88d0
        real(8), parameter :: eab = 1.5d0, sab = 0.8d0
        !! cut
        real(8), parameter :: rcut = 2.5d0
        real(8), parameter :: rcutinv = 1.d0/rcut
        !! smooth Vnew(r_cut) = 0
        !! smooth Fnew(r_cut) = 0
        !! Vnew(r) = V(r) - V(r_cut) - V'(r)*(r-r_cut)
        !! energy smooth V(r_cut) = Vcut * ess
        real(8), parameter :: Vcut = rcutinv**12 - rcutinv**6
        !! force  smooth // V'(r_cut) = - 1/sxx * [ 12/r_cut^13 - 6/r_cut^6 ] = Dvcut / sxx * ess
        real(8), parameter :: DVcut = - (12*rcutinv**13 - 6*rcutinv**7)


        associate(                    &
            natom    => tcon%natom,   &
            ra       => tcon%ra,      &
            fa       => tcon%fa,      &
            Ea       => tcon%Ea,      &
            la       => tcon%la,      &
            strain   => tcon%strain,  &
            stress   => tcon%stress,  &
            press    => tcon%press,   &
            pressxyz => tcon%pressxyz &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

       if ( present(tnb) ) then

            associate(list => tnb%list)

            do i=1, natom-1

                rai = ra(:,i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)

                    if (     tcon%r(i)==1.d0 .and. tcon%r(j)==1.d0 ) then
                        exx = eaa
                        sxx = saa
                    elseif (  tcon%r(i)==2.d0 .and. tcon%r(j)==2.d0 ) then
                        exx = ebb
                        sxx = sbb
                    else
                        exx = eab
                        sxx = sab
                    end if

                    raj = ra(:,j)
                    dra = raj - rai

                    iround = list(i)%iround(:,jj)
                    cory = list(i)%cory(jj)

                    dra = raj - rai
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )

                    if ( rij2 > (rcut*sxx)**2 ) cycle

                    rij = sqrt(rij2)

                    srinv   = sxx / rij
                    srinv3  = srinv  * srinv * srinv
                    srinv6  = srinv3 * srinv3
                    srinv12 = srinv6 * srinv6

                    Ea = Ea + exx * ( srinv12 - srinv6 - Vcut - (DVcut/sxx)*(rij-rcut*sxx) )

                    ! wij = fr * r = - V'(r) * r
                    wij  = exx * ( 12*srinv12 - 6*srinv6 + DVcut/sxx*rij )
                    wili = wili + wij

                    fr = wij / rij2

                    wilixyz = wilixyz + fr * dra(:)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            end associate

        else

            lainv = 1.d0 / la

            do i=1, natom-1

                rai = ra(:,i)

                do j=i+1, natom

                    if (     tcon%r(i)==1.d0 .and. tcon%r(j)==1.d0) then
                        exx = eaa
                        sxx = saa
                    elseif (  tcon%r(i)==2.d0 .and. tcon%r(j)==2.d0 ) then
                        exx = ebb
                        sxx = sbb
                    else
                        exx = eab
                        sxx = sab
                    end if

                    raj = ra(:,j)
                    dra = raj - rai

                    cory = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )

                    if ( rij2 > (rcut*sxx)**2 ) cycle

                    rij = sqrt(rij2)

                    srinv   = sxx / rij
                    srinv3  = srinv  * srinv * srinv
                    srinv6  = srinv3 * srinv3
                    srinv12 = srinv6 * srinv6

                    Ea = Ea + exx * ( srinv12 - srinv6 - Vcut - (DVcut/sxx)*(rij-rcut*sxx) )

                    ! wij = fr * r = - V'(r) * r
                    wij  = exx * ( 12*srinv12 - 6*srinv6 + DVcut/sxx*rij )
                    wili = wili + wij

                    fr = wij / rij2

                    wilixyz = wilixyz + fr * dra(:)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

        end if

        Ea       = econs * Ea
        fa       = econs * fa
        stress   = econs * stress  / product(la) / free
        press    = econs * wili    / product(la) / free
        pressxyz = econs * wilixyz / product(la)

        end associate
    end subroutine

    subroutine calc_force_spring( tcon, tnet )
        implicit none

        ! para list
        type(con_t),     intent(inout)          :: tcon
        type(network_t), intent(in), optional   :: tnet

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: rij2, rij, l0, fr, wij, wili, wilixyz(free), ks
        integer :: iround(free), cory
        integer :: i, j, k, ii

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            list     => tnet%sps,      &
            nlist    => tnet%nsps      &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

            do ii=1, nlist

                i      = list(ii)%i
                j      = list(ii)%j
                cory   = list(ii)%cory
                iround = list(ii)%iround
                l0     = list(ii)%l0
                ks     = list(ii)%ks

                rai = ra(:,i)
                raj = ra(:,j)

                dra = raj - rai
                dra(1) = dra(1) - cory * strain * la(free)

                do k=1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )
                rij  = sqrt( rij2 )

                Ea = Ea + 0.5d0 * ks * ( l0 - rij )**2

                wij  = ks * ( l0 - rij ) * rij
                wili = wili + wij

                fr = wij / rij2

                wilixyz = wilixyz + fr * dra(:)**2

                fa(:,j) = fa(:,j) + fr * dra
                fa(:,i) = fa(:,i) - fr * dra

                stress = stress - 2 * dra(1) * dra(2) * fr  ! 3d ?

            end do

            stress   = stress  / product(la) / free
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)
        end associate
    end subroutine

end module
