module  mo_extra_molecule
    use mo_syst
    use mo_config
    use mo_list
    use mo_force
    implicit none

    abstract interface
        logical function abstract_extra_force( tcon, tnb, ti )
            import :: con_t, list_t
            type(con_t),  intent(inout) :: tcon
            type(list_t), intent(in)    :: tnb
            integer,      intent(in)    :: ti
        end function
    end interface

    type tpextra_var
        real(8) :: forchi, chixi, vili
        real(8) :: dea
    end type
    type(tpextra_var) :: ex_var

contains



    logical function calc_extra_force( tcon, tnb, ti )
        implicit none

         ! para list
        type(con_t),  intent(inout) :: tcon
        type(list_t), intent(in)    :: tnb
        integer,      intent(in)    :: ti

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        integer :: iround(free), cory
        integer :: i, j, k, jj

        real(8), dimension(free) :: vai, vaj, vaij
        real(8) :: xxij, temp

        associate(                     &
            list     => tnb%list,      &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            la       => tcon%la,       &
            Ea       => tcon%Ea,       &
            !strain   => tcon%strain,   &
            !press    => tcon%press,    &
            !stress   => tcon%stress,   &
            !pressxyz => tcon%pressxyz, &

             va       => tcon%va,       &
             vili     => ex_var%vili,   &
             chixi    => ex_var%chixi,  &
             forchi   => ex_var%forchi, &
             dea      => ex_var%dea     &
        )

        if ( ti>0 .and. ti<=natom ) then

            dea = 0.d0
            i   = ti

            rai = ra(:,i)
            ri  = radius(i)

            do jj = 1, list(i)%nbsum

                j   = list(i)%nblist(jj)
                raj = ra(:,j)
                rj  = radius(j)

                dra = raj - rai

                do k = 1, free
                    dra(k) = dra(k) - idnint(dra(k) / la(k)) * la(k)
                end do

                rij2 = sum( dra**2 )
                dij  = ri + rj
                if ( rij2 > dij**2 ) cycle

                rij = sqrt( rij2 )
                dea = dea + ( 1.d0 - rij/dij )**alpha/alpha
            end do

            dea = dea * 0.5d0

        else

            !stress = 0.d0
            !wilixyz = 0.d0

            Ea     = 0.d0
            fa     = 0.d0
            wili   = 0.d0
            chixi  = 0.d0
            forchi = 0.d0

            do i = 1, natom
                rai = ra(:,i)
                ri  = radius(i)
                vai = va(:,i)
                do jj = 1, list(i)%nbsum
                    j = list(i)%nblist(jj)
                    if( j > i ) then
                        !iround = list(i)%iround(:,jj)
                        !cory = list(i)%cory(jj)

                        raj = ra(:,j)
                        rj  = radius(j)
                        vaj = va(:,j)

                        dra = raj - rai
                        !dra(1) = dra(1) - cory * strain * la(free)

                        do k = 1, free
                            !dra(k) = dra(k) - iround(k) * la(k)
                            dra(k) = dra(k) - idnint( dra(k) / la(k) ) * la(k)
                        end do

                        rij2 = sum( dra**2 )
                        dij  = ri + rj

                        if ( rij2 > dij**2 ) cycle

                        rij  = sqrt( rij2 )
                        temp = rij / dij

                        Ea   = Ea + ( 1.d0 - temp )**alpha/alpha
                        wij  = temp * (1.d0 - temp)**(alpha-1)
                        wili = wili + wij

                        fr      = wij / rij2
                        fa(:,j) = fa(:,j) + fr * dra
                        fa(:,i) = fa(:,i) - fr * dra

                        !wilixyz = wilixyz + fr * dra(:)**2
                        !stress = stress - 2 * dra(1) * dra(2) * fr
                        vaij   =  vaj - vai
                        xxij   = (alpha-1.d0) * temp**2 * (1.d0-temp)**(alpha-2.d0) - wij
                        chixi  = chixi  + xxij
                        forchi = forchi + xxij * sum( vaij*dra ) / rij2

                    end if
                end do
            end do

            !stress   = stress  / product(la) / free
            !pressxyz = wilixyz / product(la)
            !press    = wili    / product(la) / free
            vili = wili / dble(free)
            Ea   = Ea   / dble(natom)
        end if

        end associate

        calc_extra_force = .true.

        return

    end function

end module
