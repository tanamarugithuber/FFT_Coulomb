module FFT_coulomb
    use iso_fortran_env, only: real64
    use common_func_mod
    implicit none
    private
    integer, parameter :: dp = real64

    
    type :: FFT_coulomb_type
        ! The following arrays are allocated based on the size of the system (N) and the dimensions (3D). The rho_ref array is used to store the reference charge density, while rho_long and rho_short are used to store the long-range and short-range components of the charge density, respectively. The phi_rad array is used to store the potential calculated in real space, while phi_long and phi_short are used to store the potentials calculated in reciprocal space for the long-range and short-range components, respectively.
        !------------------------------------------------------------------
        ! charge density arrays
        !------------------------------------------------------------------
        real(dp), allocatable :: rho_l(:,:,:,:)
        real(dp), allocatable :: rho_long(:,:,:)
        real(dp), allocatable :: rho_short(:,:,:)
        !------------------------------------------------------------------
        ! potential arrays
        !------------------------------------------------------------------
        real(dp), allocatable :: phi_rad(:,:,:)
        real(dp), allocatable :: phi_long(:,:,:)
        real(dp), allocatable :: phi_short(:,:,:)
        !-----------------------------------------------------------------
        ! cofficients 
        !------------------------------------------------------------------
        real(dp), allocatable :: a_lm(:,:)
        real(dp), allocatable :: q_lm(:,:)
        !------------------------------------------------------------------
        ! parameters
        !------------------------------------------------------------------
        real(dp) :: a_skin ! skin thickness
        real(dp) :: R_nuc  ! radius of nuclear
        contains
        procedure :: compute_FFT_coulomb
        procedure :: compute_rho_l
        procedure :: compute_rho_0

    end type FFT_coulomb_type
    contains

    subroutine compute_FFT_coulomb(this)
        class(FFT_coulomb_type), intent(inout) :: this
         
        !-------------------------------------------------------------------
        ! long range term
        !-------------------------------------------------------------------
        


    
    end subroutine compute_FFT_coulomb

    function compute_rho_l(this,r,l) result(solution)
        class(FFT_coulomb_type), intent(in) :: this
        real(dp) :: solution
        real(dp), intent(in) :: r,l
        solution = 2.0_dp * r**l * this%compute_rho_0(r) / (this%R_nuc**l+r**l)

    end function compute_rho_l

    function compute_rho_0(this,r ) result(solution)
        class(FFT_coulomb_type), intent(in) :: this
        real(dp) :: solution
        real(dp), intent(in) :: r

        solution = 1.0_dp + exp((r-this%R_nuc)/this%a_skin)
        solution = 1.0_dp / solution

    end function compute_rho_0

end module FFT_coulomb