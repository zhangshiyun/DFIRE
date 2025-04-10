module mo_mode
    use mo_syst
    use mo_config
    use mo_network
    implicit none

    type matrix_t
        ! dynamic matrix; inverse matrix
        real(8), allocatable, dimension(:,:)   :: dymatrix, dymatrix0, invmatrix
        ! eigenvalues of matrix; participation ratio
        real(8), allocatable, dimension(:)     :: egdymatrix, pw
        ! vars related to works of Maloney
        real(8), allocatable, dimension(:)     :: varXi_x, varXi_y, varXi_s
        ! Tong Hua's \psi
        real(8), allocatable, dimension(:)     :: psi_th
        ! third order matrix
        real(8), allocatable, dimension(:,:,:) :: trimatrix
        ! ndim = free * m
        ! mdim = ndim + freedom of box
        integer :: natom, ndim, mdim
        ! flags of box's freedom
        logical :: boxflag   = .false.      ! can box length change
        logical :: xyflag    = .false.      ! can box length change seperately
        logical :: shearflag = .false.      ! can box be changed by shear
    contains
        ! calculate eigen-problem
        procedure :: solve   => solve_mode
        ! calculate participation ratio
        procedure :: calc_pw => calc_pw
        ! calculate inverse matrix
        procedure :: inv     => calc_inverse_matrix
    end type

    type(matrix_t) :: mode, mode0, mode1, mode2

contains
    subroutine init_mode( tmode, tcon, opboxflag, opxyflag, opshearflag )
        implicit none

        ! para list
        type(matrix_t), intent(inout)        :: tmode
        type(con_t),    intent(in)           :: tcon
        logical,        intent(in), optional :: opboxflag, opxyflag, opshearflag

        ! local

        associate(                       &
            natom     => tmode%natom,    &
            mdim      => tmode%mdim,     &
            ndim      => tmode%ndim,     &
            boxflag   => tmode%boxflag,  &
            xyflag    => tmode%xyflag,   &
            shearflag => tmode%shearflag &
            )

            boxflag   = .false.
            xyflag    = .false.
            shearflag = .false.
            if ( present(opboxflag) .and. present(opxyflag) ) then
                stop "boxflag and xyflag can not exist simultaneously"
            end if
            if ( present(opboxflag) )   boxflag   = opboxflag
            if ( present(opxyflag) )    xyflag    = opxyflag
            if ( present(opshearflag) ) shearflag = opshearflag

            natom = tcon%natom
            ndim  = free * tcon%natom
            mdim  = ndim

            if ( boxflag   ) mdim = mdim + 1
            if ( xyflag    ) mdim = mdim + free
            if ( shearflag ) mdim = mdim + 1

            if ( .not. allocated(tmode%dymatrix) ) then
                allocate( tmode%dymatrix(mdim,mdim), tmode%egdymatrix(mdim) )
               !allocate( tmode%dymatrix0(mdim,mdim) )
               !allocate( tmode%trimatrix(mdim,mdim,mdim) )
            end if

        end associate
    end subroutine

    subroutine kernel_matrix_fix( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: mij(free,free)

        ! \partial r_ij / \partial x_j
        rx = xij / rij
        ry = yij / rij

        ! \partial r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3

        ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
        ! = - \partial V_ij / [ \partial x_j \partial x_j ]
        ! =   - [ Vrr * rx**2   + vr * rxx ]
        ! yy: - [ Vrr * rx**2 + Vr * ryy ]
        ! xy: - [ Vrr * rx*ry + Vr * rxy ]
        mij(1,1) = - ( vrr*rx**2 + vr*rxx )
        mij(2,2) = - ( vrr*ry**2 + vr*ryy )
        mij(1,2) = - ( vrr*rx*ry + vr*rxy )
        mij(2,1) = mij(1,2)

          dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
        & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
          dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
        & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

          dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
          dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij
    end subroutine

    subroutine kernel_matrix_shear( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: res
        real(8) :: resx, resy
        real(8) :: reses
        real(8) :: xes, yes

        real(8) :: mij(free,free)
        real(8) :: mies(free)
        real(8) :: mee(1,1)

        ! \partial r_ij / \partial x_ij
        rx = xij / rij
        ry = yij / rij

        ! \partial r_ij / \partial \epsilon_x = rx * x_ij
        res = rx * xij

        ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]

        resx = rxx * yij  
        resy = rxy * yij  
        !
        reses = rxx * yij**2

        ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
        ! = - \partial V_ij / [ \partial x_j \partial x_j ]
        ! =   - [ Vrr * rx**2   + vr * rxx ]
        ! yy: - [ Vrr * rx**2 + Vr * ryy ]
        ! xy: - [ Vrr * rx*ry + Vr * rxy ]
        mij(1,1) = - ( vrr*rx**2 + vr*rxx )
        mij(2,2) = - ( vrr*ry**2 + vr*ryy )
        mij(1,2) = - ( vrr*rx*ry + vr*rxy )
        mij(2,1) = mij(1,2)

        mies(1) = - ( vrr*rx*res + vr*resx )
        mies(2) = - ( vrr*ry*res + vr*resy )

        mee(1,1) = vrr*res*res + vr*reses

        ! con ij
          dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
        & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
          dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
        & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

          dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
          dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

          ! ex ey
          dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) + mies
          dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+1, free*(i-1)+1:free*i ) + mies

          dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) - mies
          dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+1, free*(j-1)+1:free*j ) - mies

          dymatrix(free*natom+1:free*natom+1,free*natom+1:free*natom+1) = &
        & dymatrix(free*natom+1:free*natom+1,free*natom+1:free*natom+1) + mee
    end subroutine

    subroutine kernel_matrix_shear_angle( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij, theta )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij
        real(8), intent(in)    :: theta

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: res
        real(8) :: resx, resy
        real(8) :: reses
        real(8) :: xes, yes

        real(8) :: mij(free,free)
        real(8) :: mies(free)
        real(8) :: mee(1,1)

        ! \partial r_ij / \partial x_ij
        rx = xij / rij
        ry = yij / rij

        xes = (yij*cos(theta) - xij*sin(theta)) * cos(theta)
        yes = (yij*cos(theta) - xij*sin(theta)) * sin(theta)

        ! \partial r_ij / \partial \epsilon_x = rx * x_ij
        res = rx * xes + ry * yes

        ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]

        resx = rxx * xes - rx * sin(theta) * cos(theta) + rxy * yes - ry * sin(theta) * sin(theta)
        resy = rxy * xes + rx * cos(theta) * cos(theta) + ryy * yes + ry * cos(theta) * sin(theta)
        !
        reses = rxx * xes * xes + ryy * yes * yes

        ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
        ! = - \partial V_ij / [ \partial x_j \partial x_j ]
        ! =   - [ Vrr * rx**2   + vr * rxx ]
        ! yy: - [ Vrr * rx**2 + Vr * ryy ]
        ! xy: - [ Vrr * rx*ry + Vr * rxy ]
        mij(1,1) = - ( vrr*rx**2 + vr*rxx )
        mij(2,2) = - ( vrr*ry**2 + vr*ryy )
        mij(1,2) = - ( vrr*rx*ry + vr*rxy )
        mij(2,1) = mij(1,2)

        mies(1) = - ( vrr*rx*res + vr*resx )
        mies(2) = - ( vrr*ry*res + vr*resy )

        mee(1,1) = vrr*res*res + vr*reses

        ! con ij
          dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
        & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
          dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
        & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

          dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
          dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

          ! ex ey
          dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) + mies
          dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+1, free*(i-1)+1:free*i ) + mies

          dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) - mies
          dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+1, free*(j-1)+1:free*j ) - mies

          dymatrix(free*natom+1:free*natom+1,free*natom+1:free*natom+1) = &
        & dymatrix(free*natom+1:free*natom+1,free*natom+1:free*natom+1) + mee
    end subroutine

    subroutine kernel_matrix_compress_xy( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: rex, rey
        real(8) :: rexx, rexy, reyx, reyy
        real(8) :: rexex, rexey, reyey

        real(8) :: mij(free,free)
        real(8) :: miex(free), miey(free)
        real(8) :: mee(2,2)

        ! \partial r_ij / \partial x_ij
        rx = xij / rij
        ry = yij / rij
        ! \partial r_ij / \partial \epsilon_x = rx * x_ij
        rex = rx * xij
        rey = ry * yij

        ! \partial r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
        rexx = rxx * xij
        rexy = rxy * xij
        reyx = rxy * yij
        reyy = ryy * yij
        !
        rexex = rxx * xij**2
        reyey = ryy * yij**2
        rexey = rxy * xij*yij

        ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
        ! = - \partial V_ij / [ \partial x_j \partial x_j ]
        ! =   - [ Vrr * rx**2   + vr * rxx ]
        ! yy: - [ Vrr * rx**2 + Vr * ryy ]
        ! xy: - [ Vrr * rx*ry + Vr * rxy ]
        mij(1,1) = - ( vrr*rx**2 + vr*rxx )
        mij(2,2) = - ( vrr*ry**2 + vr*ryy )
        mij(1,2) = - ( vrr*rx*ry + vr*rxy )
        mij(2,1) = mij(1,2)


        miex(1) = - ( vrr*rx*rex + vr*rexx )
        miex(2) = - ( vrr*ry*rex + vr*rexy )
        miey(1) = - ( vrr*rx*rey + vr*reyx )
        miey(2) = - ( vrr*ry*rey + vr*reyy )

        mee(1,1) = vrr*rex*rex + vr*rexex
        mee(2,2) = vrr*rey*rey + vr*reyey
        mee(1,2) = vrr*rex*rey + vr*rexey
        mee(2,1) = mee(1,2)

        ! con ij
          dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
        & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
          dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
        & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

          dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
          dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

          ! ex ey
          dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) + miex
          dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+1, free*(i-1)+1:free*i ) + miex
          dymatrix( free*(i-1)+1:free*i, free*natom+2 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+2 ) + miey
          dymatrix( free*natom+2, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+2, free*(i-1)+1:free*i ) + miey

          dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) - miex
          dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+1, free*(j-1)+1:free*j ) - miex
          dymatrix( free*(j-1)+1:free*j, free*natom+2 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+2 ) - miey
          dymatrix( free*natom+2, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+2, free*(j-1)+1:free*j ) - miey

          dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) = &
        & dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) + mee
    end subroutine

    subroutine kernel_matrix_compress_box( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: rex, rey
        real(8) :: rexx, rexy, reyx, reyy
        real(8) :: rexex, rexey, reyey

        real(8) :: mij(free,free)
        real(8) :: miex(free), miey(free)
        !real(8) :: mee(2,2)
        real(8) :: mee

        ! \partial r_ij / \partial x_ij
        rx = xij / rij
        ry = yij / rij
        ! \partial r_ij / \partial \epsilon_x = rx * x_ij
        rex = rx * xij
        rey = ry * yij

        ! \partial r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
        rexx = rxx * xij
        rexy = rxy * xij
        reyx = rxy * yij
        reyy = ryy * yij
        !
        rexex = rxx * xij**2
        reyey = ryy * yij**2
        rexey = rxy * xij*yij

        ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
        ! = - \partial V_ij / [ \partial x_j \partial x_j ]
        ! =   - [ Vrr * rx**2   + vr * rxx ]
        ! yy: - [ Vrr * rx**2 + Vr * ryy ]
        ! xy: - [ Vrr * rx*ry + Vr * rxy ]
        mij(1,1) = - ( vrr*rx**2 + vr*rxx )
        mij(2,2) = - ( vrr*ry**2 + vr*ryy )
        mij(1,2) = - ( vrr*rx*ry + vr*rxy )
        mij(2,1) = mij(1,2)


        miex(1) = - ( vrr*rx*rex + vr*rexx )
        miex(2) = - ( vrr*ry*rex + vr*rexy )
        miey(1) = - ( vrr*rx*rey + vr*reyx )
        miey(2) = - ( vrr*ry*rey + vr*reyy )

        !mee(1,1) = vrr*rex*rex + vr*rexex
        !mee(2,2) = vrr*rey*rey + vr*reyey
        !mee(1,2) = vrr*rex*rey + vr*rexey
        !mee(2,1) = mee(1,2)
        mee = 0.d0
        mee = mee + vrr*rex*rex + vr*rexex
        mee = mee + vrr*rey*rey + vr*reyey
        mee = mee + 2 * ( vrr*rex*rey + vr*rexey )

        ! con ij
          dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
        & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
          dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
        & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

          dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
          dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

          ! ex ey
          dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) + miex + miey
          dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+1, free*(i-1)+1:free*i ) + miex + miey
        !   dymatrix( free*(i-1)+1:free*i, free*natom+2 ) = &
        ! & dymatrix( free*(i-1)+1:free*i, free*natom+2 ) + miey
        !   dymatrix( free*natom+2, free*(i-1)+1:free*i ) = &
        ! & dymatrix( free*natom+2, free*(i-1)+1:free*i ) + miey

          dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) - miex - miey
          dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+1, free*(j-1)+1:free*j ) - miex - miey
        !   dymatrix( free*(j-1)+1:free*j, free*natom+2 ) = &
        ! & dymatrix( free*(j-1)+1:free*j, free*natom+2 ) - miey
        !   dymatrix( free*natom+2, free*(j-1)+1:free*j ) = &
        ! & dymatrix( free*natom+2, free*(j-1)+1:free*j ) - miey

        !   dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) = &
        ! & dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) + mee
          dymatrix(free*natom+1,free*natom+1) = &
        & dymatrix(free*natom+1,free*natom+1) + mee
    end subroutine

    subroutine kernel_matrix_compress_shear( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: rex, rey, res
        real(8) :: rexx, rexy, reyx, reyy, resx, resy
        real(8) :: rexex, rexey, rexes, reyey, reyes, reses

        real(8) :: mij(free,free)
        real(8) :: miex(free), miey(free), mies(free)
        real(8) :: mee(3,3)

        ! \partial r_ij / \partial x_ij
        rx = xij / rij
        ry = yij / rij
        ! \partial r_ij / \partial \epsilon_x = rx * x_ij
        rex = rx * xij
        rey = ry * yij
        res = rx * yij

        ! \partial r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
        rexx = rxx * xij
        rexy = rxy * xij
        reyx = rxy * yij
        reyy = ryy * yij
        resx = rxx * yij
        resy = rxy * yij
        !
        rexex = rxx * xij**2
        rexey = rxy * xij*yij
        rexes = rxx * xij*yij  + rx * yij
        reyey = ryy * yij**2
        reyes = rxy * yij**2
        reses = rxx * yij**2

        ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
        ! = - \partial V_ij / [ \partial x_j \partial x_j ]
        ! =   - [ Vrr * rx**2   + vr * rxx ]
        ! yy: - [ Vrr * rx**2 + Vr * ryy ]
        ! xy: - [ Vrr * rx*ry + Vr * rxy ]
        mij(1,1) = - ( vrr*rx**2 + vr*rxx )
        mij(2,2) = - ( vrr*ry**2 + vr*ryy )
        mij(1,2) = - ( vrr*rx*ry + vr*rxy )
        mij(2,1) = mij(1,2)

        miex(1) = - ( vrr*rx*rex + vr*rexx )
        miex(2) = - ( vrr*ry*rex + vr*rexy )
        miey(1) = - ( vrr*rx*rey + vr*reyx )
        miey(2) = - ( vrr*ry*rey + vr*reyy )
        mies(1) = - ( vrr*rx*res + vr*resx )
        mies(2) = - ( vrr*ry*res + vr*resy )

        mee(1,1) = vrr*rex*rex + vr*rexex
        mee(2,2) = vrr*rey*rey + vr*reyey
        mee(3,3) = vrr*res*res + vr*reses
        mee(1,2) = vrr*rex*rey + vr*rexey
        mee(2,1) = mee(1,2)
        mee(1,3) = vrr*rex*res + vr*rexes
        mee(3,1) = mee(1,3)
        mee(2,3) = vrr*rey*res + vr*reyes
        mee(3,2) = mee(2,3)

        ! con ij
          dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
        & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
          dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
        & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

          dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
          dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

          ! ex ey
          dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) + miex
          dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+1, free*(i-1)+1:free*i ) + miex
          dymatrix( free*(i-1)+1:free*i, free*natom+2 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+2 ) + miey
          dymatrix( free*natom+2, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+2, free*(i-1)+1:free*i ) + miey
          dymatrix( free*(i-1)+1:free*i, free*natom+3 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+3 ) + mies
          dymatrix( free*natom+3, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+3, free*(i-1)+1:free*i ) + mies

          dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) - miex
          dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+1, free*(j-1)+1:free*j ) - miex
          dymatrix( free*(j-1)+1:free*j, free*natom+2 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+2 ) - miey
          dymatrix( free*natom+2, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+2, free*(j-1)+1:free*j ) - miey
          dymatrix( free*(j-1)+1:free*j, free*natom+3 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+3 ) - mies
          dymatrix( free*natom+3, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+3, free*(j-1)+1:free*j ) - mies

          dymatrix(free*natom+1:free*natom+3,free*natom+1:free*natom+3) = &
        & dymatrix(free*natom+1:free*natom+3,free*natom+1:free*natom+3) + mee
    end subroutine

    subroutine make_dymatrix( tmode, tcon, theta,  tnet, opflag )
        implicit none

        ! para list
        type(matrix_t),  intent(inout)        :: tmode
        type(con_t),     intent(in)           :: tcon
        real(8),     intent(in)                 :: theta
        type(network_t), intent(in), optional :: tnet
        integer,         intent(in), optional :: opflag

        ! local
        integer :: i, j, ii
        real(8) :: dra(free), rij, rij2, dij
        real(8) :: vr, vrr
        real(8) :: xij, yij
        real(8) :: ks, l0

        real(8) :: srinv, srinv3, srinv6, srinv12

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

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tmode%natom,    &
            mdim     => tmode%mdim,     &
            ndim     => tmode%ndim      &
            )

            dymatrix = 0.d0

            if ( .not. present( tnet ) ) then

                do i=1, natom-1
                    do j=i+1, natom

                        if ( tcon%r(i)==1.d0 .and. tcon%r(j)==1.d0 ) then
                            exx = eaa
                            sxx = saa
                        elseif ( tcon%r(i)==2.d0 .and. tcon%r(j)==2.d0  ) then
                            exx = ebb
                            sxx = sbb
                        else
                            exx = eab
                            sxx = sab
                        end if

                        dra  = calc_dra( tcon, i, j )
                        rij2 = sum(dra**2)

                        if ( rij2 > (rcut*sxx)**2 ) cycle

                        rij  = sqrt(rij2)

                        srinv   = sxx / rij
                        srinv3  = srinv  * srinv * srinv
                        srinv6  = srinv3 * srinv3
                        srinv12 = srinv6 * srinv6


                        ! \partial V_ij / \partial r_ij
                        vr = - econs * exx * (12*srinv12 - 6*srinv6 + DVcut/sxx*rij) / rij
                        ! \partial^2 V_ij / \partial r_ij^2
                        vrr = econs * exx * (12*13*srinv12 / rij2 - 42*srinv6 / rij2)

                        xij = dra(1)
                        yij = dra(2)

                        !if ( present( opflag ) ) then
                        !    write(*,*) "oh.."
                        !    stop
                        !end if
                        !call kernel_matrix_fix( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        if ( .not. present( opflag ) ) then
                            call kernel_matrix_fix( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        elseif ( opflag == 1 ) then
                            call kernel_matrix_compress_xy( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        elseif ( opflag == 11 ) then
                            call kernel_matrix_compress_box( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        elseif ( opflag == 2 ) then
                            call kernel_matrix_shear_angle( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij,theta )
                        elseif ( opflag == 3 ) then
                            call kernel_matrix_compress_shear( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        end if

                    end do
                end do

                if (present(opflag) .and. opflag ==3) then
                    dymatrix(free*natom+1,free*natom+2) = dymatrix(free*natom+1, free*natom+2) + tcon%pressxyz(1) * product(tcon%la)
                    dymatrix(free*natom+2,free*natom+1) = dymatrix(free*natom+2, free*natom+1) + tcon%pressxyz(2) * product(tcon%la)
                end if


            else
                associate(             &
                    list  => tnet%sps, &
                    nlist => tnet%nsps &
                    )

                    do ii=1, nlist

                        i  = list(ii)%i
                        j  = list(ii)%j
                        l0 = list(ii)%l0
                        ks = list(ii)%ks

                        dra = calc_dra( tcon, i, j )
                        rij2 = sum(dra**2)

                        rij  = sqrt(rij2)

                        ! \partial V_ij / \partial r_ij
                       !vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                        vr = 0.d0
                        ! \partial^2 V_ij / \partial r_ij^2
                       !vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2
                        vrr = ks

                        xij = dra(1)
                        yij = dra(2)

                        if ( .not. present( opflag ) ) then
                            call kernel_matrix_fix( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        elseif ( opflag == 1 ) then
                            call kernel_matrix_compress_xy( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        elseif ( opflag == 11 ) then
                            call kernel_matrix_compress_box( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        elseif ( opflag == 2 ) then
                            call kernel_matrix_shear( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        elseif ( opflag == 3 ) then
                            call kernel_matrix_compress_shear( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        end if

                    end do

                    if ( present( opflag ) ) then
                        if ( opflag == 1 ) then
                            tmode%varXi_x = tmode%dymatrix(1:ndim,ndim+1)
                            tmode%varXi_y = tmode%dymatrix(1:ndim,ndim+2)
                        elseif ( opflag == 2 ) then
                            tmode%varXi_s = tmode%dymatrix(1:ndim,ndim+1)
                        elseif ( opflag == 3 ) then
                            tmode%varXi_x = tmode%dymatrix(1:ndim,ndim+1)
                            tmode%varXi_y = tmode%dymatrix(1:ndim,ndim+2)
                            tmode%varXi_s = tmode%dymatrix(1:ndim,ndim+3)
                        end if
                    end if

                end associate

            end if

        end associate
    end subroutine

    subroutine calc_psi_th( tmode, istart )
        !
        !  \phi_th^i = \sum_\alpha \frac{1}{\omega_alpha^2} ||e_\omega_\alpha^i||^2
        !
        implicit none

        ! para list
        type(matrix_t), intent(inout) :: tmode
        integer, optional, intent(in) :: istart

        ! local
        integer :: lc_istart
        integer :: i, j
        real(8) :: temp

        if ( .not. allocated( tmode%psi_th ) ) then
            allocate( tmode%psi_th( tmode%natom ) )
        end if

        associate(                          &
            dymatrix   => tmode%dymatrix,   &
            egdymatrix => tmode%egdymatrix, &
            psi_th     => tmode%psi_th,     &
            ndim       => tmode%ndim,       &
            natom      => tmode%natom       &
            )

            psi_th = 0.d0

            if ( .not. present( istart ) ) then
                lc_istart = free + 1
            else
                lc_istart = istart
            end if

            do j=lc_istart, ndim
                if (egdymatrix(j) < 1.e-8) then
                    continue
                end if
                do i=1, natom
                    temp = 1.d0 / egdymatrix(j) * sum( dymatrix(free*(i-1)+1:free*i, j)**2 )
                    psi_th(i) = psi_th(i) + temp
                end do
            end do

        end associate
    end subroutine

    subroutine make_trimatrix( tmode, tcon )
        implicit none

        ! para list
        type(matrix_t), intent(inout) :: tmode
        type(con_t), intent(in) :: tcon

        ! local
        integer :: i, j
        real(8) :: dra(free), rij, rij2, rij3, rij5, dij, xi, xj, yi, yj
        real(8) :: vr, vrr!, vrr
        real(8) :: rx, ry
        real(8) :: rxx, ryy, rxy
        real(8) :: rxxx, rxxy, rxyy, ryyy
        real(8) :: mijk(free,free,free)   ! \p^2v / ( \pxi * \pxj ) ....

        associate(                        &
            trimatrix => tmode%trimatrix, &
            radius    => tcon%r,          &
            natom     => tcon%natom       &
            )

            trimatrix = 0.d0

            do i=1, natom-1
                do j=i+1, natom

                    dra = calc_dra( tcon, i, j )
                    rij2 = sum(dra**2)

                    dij = radius(i) + radius(j)

                    if ( rij2 > dij**2 ) cycle

                    rij  = sqrt(rij2)
                    rij3 = rij2*rij
                    rij5 = rij3*rij2

                    !vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                    !vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2
                    ! alpha = 2
                    vr   = - ( 1.d0 - rij/dij ) / dij
                    vrr  = 1.d0 / dij**2
                    !vrrr = 0.d0

                    xj = dra(1); yj = dra(2)
                    xi =-dra(1); yi =-dra(2)

                    rx = xj / rij
                    ry = yj / rij

                    rxx = 1.d0/rij - rx**2/rij3
                    ryy = 1.d0/rij - ry**2/rij3
                    rxy =          - rx*ry/rij3

                    ! \pr / ( \pxj \pxj \pxj )
                    rxxx = 3.d0*xj/rij3 + xj**3/rij5
                    ryyy = 3.d0*yj/rij3 + yj**3/rij5
                    rxxy = -yj/rij3 + 3.d0*xj**2*yj/rij5
                    rxyy = -xj/rij3 + 3.d0*xj*yj**2/rij5

                    ! \pv / ( \pxi \pxj \pxj )
                    mijk(1,1,1) = -( vrr*3.d0*rxx*rx - vr*rxxx )
                    mijk(1,1,2) = -( vrr*(2.d0*rxx*ry+rx*rxy) - vr*rxxy )
                    mijk(1,2,1) = mijk(1,1,2)
                    mijk(1,2,2) = -( vrr*(2.d0*ryy*rx+ry*rxy) - vr*rxyy )
                    mijk(2,1,1) = -mijk(1,1,2)
                    mijk(2,1,2) = -mijk(1,2,2)
                    mijk(2,2,1) = -mijk(1,2,2)
                    mijk(2,2,2) = -( vrr*3.d0*ryy*ry - vr*ryyy )

                    ! i
                    trimatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) =  mijk
                    trimatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = -mijk
                    trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = -mijk
                    ! j
                    trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(i-1)+1:free*i ) =  mijk
                    trimatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = -mijk
                    trimatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i, free*(j-1)+1:free*j ) =  mijk
                    ! iii, jjj
                    trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
                        trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) + mijk
                    trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
                        trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mijk

                end do
            end do

        end associate
    end subroutine

    subroutine solve_mode( this, oprange )
        implicit none

        ! para list
        class(matrix_t) :: this
        integer, optional :: oprange

        ! local
        integer :: rangevar

        rangevar = 0
        if ( present( oprange ) ) rangevar = oprange

        this%dymatrix0 = this%dymatrix

        associate(                        &
            mdim       => this%mdim,      &
            dymatrix   => this%dymatrix,  &
            egdymatrix => this%egdymatrix &
            )

            call solve_matrix( dymatrix, mdim, egdymatrix, rangevar )

        end associate
    end subroutine

    subroutine calc_inverse_matrix( this )
        implicit none

        ! para list
        class(matrix_t), intent(inout) :: this
        !integer, intent(in)            :: nu_ratter

        ! local
        integer :: i

        this%invmatrix = this%dymatrix

        ! m0 = Matrix
        ! m1 = eig(m0)
        ! m_inv = \sum_{i_real} eigenvalue(i) * e_{\omega_i} * e_{\omega_i}^T
        do i=1, this%mdim
            !if ( i<=free*(1) ) then
            if ( this%egdymatrix(i) < 1.d-10 ) then
                this%invmatrix(:,i) = 0.d0
            else
                this%invmatrix(:,i) = ( 1.d0 / this%egdymatrix(i) ) * this%invmatrix(:,i)
            end if
        end do

        this%invmatrix = matmul(this%invmatrix, transpose(this%dymatrix))
    end subroutine

    subroutine calc_pw( this )
        implicit none

        ! para list
        class(matrix_t) :: this

        ! local
        integer :: i, j
        real(8) :: sum4

        if ( .not. allocated(this%pw) ) allocate( this%pw(this%mdim) )

        associate(                    &
            mdim   => this%mdim,      &
            natom  => this%natom,     &
            pw     => this%pw,        &
            dymatrix => this%dymatrix &
            )

            do i=1, mdim
                sum4 = 0.d0
                do j=1, natom
                    sum4 = sum4 + sum( dymatrix(free*(j-1)+1:free*j,i)**2 )**2
                end do
                if ( mdim > natom * free ) then
                    sum4 = sum4 + sum( dymatrix(free*natom+1:mdim,i)**4 )
                end if
                pw(i) = 1.d0 / dble(natom) / sum4
            end do

        end associate
    end subroutine

    subroutine solve_matrix(a,order,b, rangevar)
        implicit none

        ! para list
        integer :: order
        real(8),dimension(1:order,1:order) :: a
        real(8),dimension(1:order) :: b
        integer :: rangevar

        ! local
        character :: jobz = 'V'
        character :: range = 'A'
        character :: uplo = 'U'

        integer :: n
        !real, dimension  :: a
        integer :: lda
        integer :: vl, vu, il, iu   ! will not referenced when range=A
        real(8) :: abstol           ! important variable
        integer :: m
        !real, dimension :: w        ! use b above ! output eigenvalues
        real(8), allocatable, dimension(:,:) :: z
        integer :: ldz
        integer, allocatable, dimension(:) :: isuppz
        real(8), allocatable, dimension(:) :: work
        integer :: lwork
        integer, allocatable, dimension(:) :: iwork
        integer :: liwork
        integer :: info

        n = order
        lda = order
        ldz= order

        if( rangevar == 0 ) then
            range = 'A'
        elseif( rangevar > 0 .and. rangevar <= order ) then
            range = 'I'
            il = 1
            iu = rangevar
        else
            write(*,*) "error rangevar set"
            stop
        end if

        !--
        abstol = -1
        !abstol = dlamch('S')        ! high precision

        !allocate(a(order,order)); allocate(b(order))
        allocate(z(order,order))
        allocate(isuppz(2*order))


        !- query the optimal workspace
        allocate(work(1000000)); allocate(iwork(1000000))
        lwork  = -1
        liwork = -1

        call dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,b,z,       &
                   & ldz, isuppz, work, lwork, iwork, liwork, info)

        lwork  = int(work(1))
        liwork = int(iwork(1))

        deallocate(work); deallocate(iwork)

        allocate(work(lwork)); allocate(iwork(liwork))

        !vvv- main

        call dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,b,z,       &
                   & ldz, isuppz, work, lwork, iwork, liwork, info)

        if(info .ne. 0) then
            write(*,*) '** Find error in subroutine diago'
            stop
        end if
        !^^^- done

        deallocate(work); deallocate(iwork)

        !- output results via a
        a = z

        !deallocate(a); deallocate(b)
        deallocate(z); deallocate(isuppz)
    end subroutine

end module
