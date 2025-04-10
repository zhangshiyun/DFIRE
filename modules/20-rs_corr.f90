module mo_rs_corr
    use mo_syst
    use mo_config
    use mo_math
    use mo_field
    implicit none
!This module is to compute two-point correlation between two (scalar, vector ...) fields $i$ and $j$ with the different/same particles/sites $id$ at the different/same time $t$ in the real space;    
!Assumption: 1) isotropy and translation symmetry. 2) point particles/sites;  
!Warning: 1) The output corr(r) is normalized by the averaged field. 2) The boundary conditions and length units of two fields shall be (approximately) same; 

!for visiting
!input
    type, public :: tpset_corr 
        logical  :: calc_flag    
        real(8)  :: unit_d                                !the unit length for the result field corr(r)            
        real(8)  :: rbin, rmin, rmax                      !the bin and cutoff normalized by unit_d 
        integer  :: idi_num, idj_num                      !the numbers of particles/sites of two fields respectively          
        integer, allocatable, dimension(:) :: idi         !the relevant particle id of the field i  
        integer, allocatable, dimension(:) :: idj         !## j  
        character(250) :: filename                        !the output file for corr(r)
        contains
        procedure :: init_setcorr  => init_setcorr 
        procedure :: clean_setcorr => clean_setcorr     
    end type

!output: the normalized correlation corr(r)
    type, public :: tpout_corr 
        real(8), allocatable, dimension(:) :: absr    
        real(8), allocatable, dimension(:) :: corr   
        contains
        procedure :: clean_outcorr => clean_outcorr    
    end type

!for computing
    type, abstract, public :: Base_corr    
! read from the settings
        real(8) :: unit_d             
        real(8) :: rbin, rmin, rmax   
        integer :: nbin               
        integer :: idi_num, idj_num         
        integer, allocatable, dimension(:) :: idi         !the relevant particle id of the field i  
        integer, allocatable, dimension(:) :: idj         !## j  
        
        integer :: ave_times                              !the effective times of calling body 

        real(8), allocatable, dimension(:) :: rvol  
        real(8), allocatable, dimension(:) :: accu_corr
    contains
        procedure :: initial               => initial
        procedure :: output                => output
        procedure :: clean                 => clean 
        procedure(abstr_body), deferred    :: body       
    end type


    type, extends(Base_corr), public :: tpcorr_vec
        contains
        procedure, pass  :: body  =>  body_corr_vec
    end type

    type, extends(Base_corr), public :: tpcorr_sca
        contains
        procedure, pass  :: body  =>  body_corr_sca
    end type

    abstract interface 

        subroutine abstr_body( this, tfield_i, tfield_j, tcalc_flag ) 
            import :: Base_corr, tpfield
            class(Base_corr),  intent(inout) :: this
            class(tpfield),    intent(in)    :: tfield_i
            class(tpfield),    intent(in)    :: tfield_j
            logical,         intent(in)      :: tcalc_flag
        end subroutine

    end interface


contains

    subroutine initial( this, tset_corr )
        implicit none

        ! para list
        class(Base_corr),  intent(inout) :: this
        type(tpset_corr),  intent(inout) :: tset_corr

        !local 
        integer :: i
        real(8) :: rad
        
        this%ave_times = 0

        associate(                         &
            unit_d     =>   this%unit_d,   &
            rbin       =>   this%rbin,     &   
            rmin       =>   this%rmin,     &
            rmax       =>   this%rmax,     &
            nbin       =>   this%nbin,     &     
            idi_num    =>   this%idi_num,  &
            idj_num    =>   this%idj_num   &
            )  
   
        unit_d    = tset_corr%unit_d
        rbin      = tset_corr%rbin   
        rmin      = tset_corr%rmin
        rmax      = tset_corr%rmax
        idi_num   = tset_corr%idi_num  
        idj_num   = tset_corr%idj_num  

        if ( .not. allocated(this%idi) ) allocate( this%idi(idi_num) )
        if ( .not. allocated(this%idj) ) allocate( this%idj(idj_num) )
        this%idi  =  tset_corr%idi
        this%idj  =  tset_corr%idj 
        call tset_corr%clean_setcorr()

        nbin = int( (rmax-rmin)/rbin )        
        if ( .not. allocated(this%rvol) )      allocate( this%rvol(nbin) )
        if ( .not. allocated(this%accu_corr) ) allocate( this%accu_corr(nbin) )

        associate(                         &
            rvol      =>   this%rvol,      &
            accu_corr =>   this%accu_corr  &        
            )

        rvol      = 0.d0
        accu_corr = 0.d0

        do i = 1, nbin
            rad     = rmin + dble(i-1) * rbin 
            rvol(i) = (rad+rbin)**free - rad**free       
        end do

        rvol = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * rvol 

        end associate   
        end associate  
  
    end subroutine

    subroutine body_corr_sca( this, tfield_i, tfield_j, tcalc_flag )
        implicit none

        ! para list
        class(tpcorr_sca), intent(inout) :: this
        class(tpfield),    intent(in)    :: tfield_i
        class(tpfield),    intent(in)    :: tfield_j
        logical,         intent(in)      :: tcalc_flag

        !local 
        integer :: i, j, k, step
        real(8) :: aryi, aryj, avei, avej
        real(8) :: dra(free), rai(free), raj(free), rij
        real(8) :: factor, vol             !factor = tfield%unit_d / unit_d, the volume in the unit of unit_d**free
        real(8) :: rbin, rmin                                   
        real(8), allocatable, dimension(:) :: corr, rvol 

        if( tcalc_flag .eqv. .false.) return 

        select type(tfield_i)        
        type is (tpfield_sca)
        select type(tfield_j)        
        type is (tpfield_sca)

        associate(                         &
            la        =>  tfield_i%la,     &
            funit_d   =>  tfield_i%unit_d, &
            ra_i      =>  tfield_i%ra,     &
            ra_j      =>  tfield_j%ra,     &
            ary_i     =>  tfield_i%ary,    &
            ary_j     =>  tfield_j%ary,    &
            idi_num   =>  this%idi_num,    &
            idj_num   =>  this%idj_num,    &   
            idi       =>  this%idi,        &
            idj       =>  this%idj,        &  
            unit_d    =>  this%unit_d,     &
            nbin      =>  this%nbin,       &  
            accu_corr =>  this%accu_corr   &        
            )

        allocate( corr(nbin) )
        allocate( rvol(nbin) )

        factor = unit_d / funit_d
        rvol   = this%rvol * factor**free 
        rbin   = this%rbin * factor 
        rmin   = this%rmin * factor 
        vol    = product( la )
        corr   = 0.d0
        avei   = 0.d0
        avej   = 0.d0

        do i = 1, idi_num
            rai  = ra_i(:,idi(i))
            aryi = ary_i(idi(i))
            avei = avei + aryi

            do j = 1, idj_num
                raj  = ra_j(:,idj(j))
                aryj = ary_j(idj(j))

                dra = rai - raj
                do k = 1, free
                    dra(k) = dra(k) - idnint(dra(k) / la(k)) * la(k)
                end do

                rij = sqrt( sum( dra**2 ) ) 
                rij = rij - rmin
                step = int( rij / rbin )  

                if( step.gt.0 .and. step.le.nbin ) then            
                    corr(step) = corr(step) + aryi * aryj
                endif                        
            end do
        end do

        do j = 1, idj_num
            aryj = ary_j(idj(j))
            avej = avej + aryj   
        end do

        corr = vol * corr / ( avei*avej )
        corr = corr / rvol

        accu_corr = accu_corr + corr

        end associate   
        end select 
        end select
    
        this%ave_times  = this%ave_times + 1

        deallocate( corr, rvol )
    
    end subroutine


    subroutine body_corr_vec( this, tfield_i, tfield_j, tcalc_flag )
        implicit none

        ! para list
        class(tpcorr_vec), intent(inout) :: this
        class(tpfield),    intent(in)    :: tfield_i
        class(tpfield),    intent(in)    :: tfield_j
        logical,         intent(in)      :: tcalc_flag

        !local 
        integer :: i, j, k, step
        real(8) :: dra(free), rai(free), raj(free), rij
        real(8), allocatable, dimension(:) :: aryi, aryj, avei, avej
        real(8) :: factor, vol             !factor = tfield%unit_d / unit_d, the volume in the unit of unit_d**free
        real(8) :: rbin, rmin                                   
        real(8), allocatable, dimension(:) :: corr, rvol 

        if( tcalc_flag .eqv. .false.) return 

        select type(tfield_i)        
        type is (tpfield_vec)
        select type(tfield_j)        
        type is (tpfield_vec)

        associate(                         &
            degree    =>  tfield_i%degree, &
            la        =>  tfield_i%la,     &
            funit_d   =>  tfield_i%unit_d, &
            ra_i      =>  tfield_i%ra,     &
            ra_j      =>  tfield_j%ra,     &
            ary_i     =>  tfield_i%ary,    &
            ary_j     =>  tfield_j%ary,    &
            idi_num   =>  this%idi_num,    &
            idj_num   =>  this%idj_num,    &   
            idi       =>  this%idi,        &
            idj       =>  this%idj,        &  
            unit_d    =>  this%unit_d,     &
            nbin      =>  this%nbin,       &     
            accu_corr =>  this%accu_corr   &        
            )

        allocate( corr(nbin) )
        allocate( rvol(nbin) )
        allocate( aryi(degree) )
        allocate( aryj(degree) )
        allocate( avei(degree) )
        allocate( avej(degree) )

        factor = unit_d / funit_d
        rvol   = this%rvol * factor**free 
        rbin   = this%rbin * factor 
        rmin   = this%rmin * factor 
        vol    = product( la )
        corr   = 0.d0
        avei   = 0.d0
        avej   = 0.d0

        do i = 1, idi_num
            rai  = ra_i(:,idi(i))
            aryi = ary_i(:,idi(i))
            avei = avei + aryi

            do j = 1, idj_num
                raj  = ra_j(:,idj(j))
                aryj = ary_j(:,idj(j))

                dra = rai - raj
                do k = 1, free
                    dra(k) = dra(k) - idnint(dra(k) / la(k)) * la(k)
                end do

                rij = sqrt( sum( dra**2 ) ) 
                rij = rij - rmin
                step = int( rij / rbin )  

                if( step.gt.0 .and. step.le.nbin ) then            
                    corr(step) = corr(step) + sum( aryi * aryj )
                endif                        
            end do
        end do

        do j = 1, idj_num
            aryj = ary_j(:,idj(j))
            avej = avej + aryj    
        end do

        corr = vol * corr / sum( avei * avej )
        corr = corr / rvol

        accu_corr = accu_corr + corr

        end associate   
        end select 
        end select
    
        this%ave_times  = this%ave_times + 1

        deallocate( corr, rvol )
        deallocate( aryi, aryj, avei, avej )
    
    end subroutine

    subroutine output( this, tout_corr, filename )
        implicit none

        ! para list
        class(Base_corr),  intent(in)        :: this
        type(tpout_corr), intent(inout)      :: tout_corr
        character(250), optional, intent(in) :: filename
     
        !local 
        integer :: i

         associate(                       &
            nbin      =>  this%nbin,      &
            rbin      =>  this%rbin,      &
            rmin      =>  this%rmin,      &
            ave_times =>  this%ave_times, & 
            accu_corr =>  this%accu_corr  &
            )  

        if(ave_times < 1) then
            print *, "fail to output, no calculation of correlation is done"
            return
        end if

        if ( .not. allocated(tout_corr%corr) )  allocate( tout_corr%corr(nbin) ) 
        if ( .not. allocated(tout_corr%absr) )  allocate( tout_corr%absr(nbin) )  

        associate(                    &
            corr  =>  tout_corr%corr, & 
            absr  =>  tout_corr%absr  &    
            )     
                
        do i = 1, nbin
            absr(i) = dble(i-1) * rbin + rmin     
        end do
        corr = accu_corr / dble(ave_times)

        if( present( filename ) ) then
        open( 1, file = filename )
        do i = 1, nbin
                write(1,*) absr(i), corr(i)
        end do
        close(1)
        end if         

        end associate
        end associate
    end subroutine 

    subroutine clean( this)
        implicit none

        ! para list
        class(Base_corr),  intent(inout) :: this
                    
        if ( allocated(this%idi) )       deallocate( this%idi )
        if ( allocated(this%idj) )       deallocate( this%idj )
        if ( allocated(this%rvol) )      deallocate( this%rvol )
        if ( allocated(this%accu_corr) ) deallocate( this%accu_corr )
    end subroutine


    subroutine init_setcorr( this )
        implicit none

        ! para list
        class(tpset_corr), intent(inout) :: this
                    
        if ( .not. allocated(this%idi) )  allocate( this%idi(this%idi_num) )
        if ( .not. allocated(this%idj) )  allocate( this%idj(this%idj_num) )
    end subroutine


    subroutine clean_setcorr( this)
        implicit none

        ! para list
        class(tpset_corr), intent(inout) :: this
                    
        if ( allocated(this%idi) ) deallocate( this%idi )
        if ( allocated(this%idj) ) deallocate( this%idj )
    end subroutine

    subroutine clean_outcorr( this)
        implicit none

        ! para list
        class(tpout_corr), intent(inout) :: this
                    
        if ( allocated(this%corr) ) deallocate( this%corr )
        if ( allocated(this%absr) ) deallocate( this%absr )
    end subroutine

end module
