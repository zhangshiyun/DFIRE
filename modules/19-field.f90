module mo_field
    use mo_syst

    type, public :: tpfield 
        integer  :: num                                         !the number of the points in the field 
        real(8)  :: unit_d                                      !the unit of length of the field, which will be 1, if the ra is normalized 
        logical  :: ra_flag = .false.
        logical  :: af_flag = .false.
        real(8), pointer, dimension(:)   :: la  => null()       !the boundry condition of the field
        real(8), pointer, dimension(:,:) :: ra  => null()
        contains            
        procedure :: init_field  => init_field
        procedure :: clean_field => clean_field
    end type

    type, extends(tpfield), public :: tpfield_vec
        integer  :: degree                                      !the dimensionality of the field
        real(8), allocatable, dimension(:) :: unit_f            !the unit of the field, which will be 1, if vector is normalized              
        real(8), pointer, dimension(:,:)   :: ary => null()  !the field 
        contains             
        procedure :: init_field  => init_field_vec
        procedure :: clean_field => clean_field_vec
    end type

    type, extends(tpfield), public :: tpfield_sca
        real(8) :: unit_f
        real(8), pointer, dimension(:) :: ary => null()  
        contains             
        procedure :: init_field  => init_field_sca
        procedure :: clean_field => clean_field_sca
    end type
    
contains

    subroutine init_field( this, all_flag )
        implicit none

        ! para list    
        class(tpfield), intent(inout) :: this
        logical, optional, intent(in) :: all_flag  

        !local        
        integer :: i

        this%ra_flag = .false.

        if( present( all_flag ) .and. all_flag .eqv. .true.) then   
            this%ra_flag = .true.      
            allocate( this%la(free) )
            allocate( this%ra(free, this%num) )
            this%la = 0.d0
            do i = 1, free 
                this%ra(:,i) = 0.d0
            end do
        endif     
        
    end subroutine


    subroutine init_field_sca( this, all_flag )
        implicit none

        ! para list    
        class(tpfield_sca), intent(inout) :: this
        logical, optional, intent(in)     :: all_flag  

        !local        
        integer :: i

        this%ra_flag = .false.
        this%af_flag = .true.

        if( present( all_flag ) .and. all_flag .eqv. .true.) then     
            this%ra_flag = .true.    
            allocate( this%la(free) )
            allocate( this%ra(free, this%num) )
            this%la = 0.d0
            do i = 1, free 
                this%ra(:,i) = 0.d0
            end do
        endif     

        allocate( this%ary(this%num) )
        this%ary  = 0.d0   
        this%unit_f = 0.d0

    end subroutine

    subroutine init_field_vec( this, all_flag )
        implicit none

        ! para list    
        class(tpfield_vec), intent(inout) :: this
        logical, optional, intent(in)     :: all_flag  

        !local        
        integer :: i
        
        this%ra_flag = .false.
        this%af_flag = .true.

        if( present( all_flag ) .and. all_flag .eqv. .true.) then     
            this%ra_flag = .true.       
            allocate( this%la(free) )
            allocate( this%ra(free, this%num) )
            this%la = 0.d0
            do i = 1, free 
                this%ra(:,i) = 0.d0
            end do
        endif     

        allocate( this%ary(this%degree, this%num) )    
        do i = 1, this%degree 
            this%ary(:,i) = 0.d0
        end do

        allocate( this%unit_f(this%degree) )
        this%unit_f = 0.d0    

    end subroutine
        
    subroutine clean_field( this )
        implicit none

        ! para list    
        class(tpfield), intent(inout) :: this

        if(this%ra_flag .eqv. .true.) then
            deallocate( this%la )
            deallocate( this%ra )
        else
            this%la => null()
            this%ra => null()
        end if

    end subroutine

    subroutine clean_field_sca( this )
        implicit none

        ! para list    
        class(tpfield_sca), intent(inout) :: this

        if(this%ra_flag .eqv. .true.) then
            deallocate( this%la )
            deallocate( this%ra )
        else
            this%la => null()
            this%ra => null()
        end if

        if(this%af_flag .eqv. .true.) then
            deallocate( this%ary )
        else
            this%ary => null()
        end if
    end subroutine

    subroutine clean_field_vec( this )
        implicit none

        ! para list    
        class(tpfield_vec), intent(inout) :: this

        if(this%ra_flag .eqv. .true.) then
            deallocate( this%la )
            deallocate( this%ra )
        else
            this%la => null()
            this%ra => null()
        end if

        if(this%af_flag .eqv. .true.) then
            deallocate( this%ary )
            deallocate( this%unit_f )
        else
            this%ary    => null()
        end if
    end subroutine

end module
