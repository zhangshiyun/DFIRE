module  mo_molecule
    use mo_syst
    use mo_config
    use mo_list
    use mo_force
    implicit none
!This molecule is for the template of molecule dynamical simulation
!Application scope: 1) Basically, it applies for NVT/NPT ensemble or simplier NVE/NPE ensembles with various dynamics 2) Via extensions, one could also simulate the grand canonical ensemble(open system), or add more functions, such as swap, etc.

! For visiting
! 1) input: for setting and resetting
    type, public :: tpset_molecule
        real(8)  :: dt             !the time step in the simulation
        real(8)  :: temper         !temperature
        real(8)  :: resc_prob      !probability of rescaling per step
    end type

    type, extends(tpset_molecule), public :: tpset_npt
        real(8)  :: pre            !pressure
    end type

    type, extends(tpset_npt), public :: tpset_npt_swap
        real(8)  :: swap_prob      !probability of swap per particle
    end type

! 2) output: instant qualities
    type, public :: tpout_molecule
        real(8)  :: Ea, Ek, Ev     !the kinetic energy, potential energy and total energy
        real(8)  :: T, press, phi  !temperature, pressure and volume fraction
    end type

    type, abstract, public :: Base_molecule
        contains
! for calling in the program main
        procedure(abstr_initial),      deferred :: initial
        procedure(abstr_body),         deferred :: body
        procedure(abstr_clean),        deferred :: clean
! black box
        !predic, correc, constri, rescale_p, rescale_t, ...
! external procedure
        !make_list, check_list, calc_force, ...
    end type

    !procedure(abstract_extra_force), pointer :: p_abstr_extra_force => null()
    !procedure(abstract_fun_force),   pointer :: p_abstr_fun_force   => null()
    !...

    abstract interface

        subroutine abstr_initial( this, tcon, tnb, tp_abstr_force, tset_molecule ) !Function:
!1) read in the parameters and the force type for the molecules_dynamics, via set_molecule and tp_abstr_force, respectively;
!2) initialize/restart the internal parameters via settings;
!3) initialize/restart the velocity, list and force;
            import :: con_t, list_t, Base_molecule, tpset_molecule
            class(Base_molecule), intent(inout) :: this
            type(con_t),          intent(inout) :: tcon
            type(list_t),         intent(inout) :: tnb
            class(tpset_molecule),intent(in)    :: tset_molecule
            logical,              external      :: tp_abstr_force
        end subroutine

        subroutine abstr_body( this, tcon, tnb, tp_abstr_force, tout_molecule ) !Function:
! performing the dynamics of N = natom particles and other functions.
! if tout_molecule is presented, it will output the instant qualities defined in out_molecule, which can be visitted via tout_molecule%
            import :: con_t, list_t, Base_molecule, tpout_molecule
            class(Base_molecule), intent(inout) :: this
            type(con_t),          intent(inout) :: tcon
            type(list_t),         intent(inout) :: tnb
            logical,              external      :: tp_abstr_force
            type(tpout_molecule), intent(inout), optional :: tout_molecule
        end subroutine

! release the memory
        subroutine abstr_clean( this )
            import :: Base_molecule
            class(Base_molecule), intent(inout) :: this
        end subroutine

    end interface

end module


