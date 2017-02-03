!*******************************************************************************
!>  @file gp_kernel.f
!>  @brief Contains module @ref gp_kernel.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref gpkernel_class. This module contains
!>  all the code necessary to define the kernel for nonparameterized profiles
!>  calculated using a gaussian process
!*******************************************************************************
      MODULE gp_kernel
      USE stel_kinds
      USE stel_constants
      USE v3_utilities
      USE mpi_inc
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  model module parameters
!*******************************************************************************
!>  Maximum size for kernel name lengths.
      INTEGER, PARAMETER, PUBLIC :: k_type_len=20

!>  Lower array bound for function profiles.
      INTEGER, PARAMETER, PUBLIC :: ilb_gp   = 0
!>  Upper array bound for function profiles.
      INTEGER, PARAMETER, PUBLIC :: iub_gp   = 21

!>  No profile type.
      INTEGER, PARAMETER, PRIVATE :: gpkernel_none_type         = -1
!>  One dimensional squared exponential type.
      INTEGER, PARAMETER, PRIVATE :: gpkernel_oned_sqexp_type   = 0
!>  Three dimensional squared exponential type. (Example)
!      INTEGER, PARAMETER, PRIVATE :: gpkernel_threed_sqexp_type = 1


!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) model class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a parameterized profile.
!-------------------------------------------------------------------------------
      TYPE gpkernel_class
!>  Gaussian process kernel type
         INTEGER :: k_type = gpkernel_none_type 
!>  Array of coefficients for the kernel functions.
         REAL(rprec), DIMENSION(ilb_gp:iub_gp) :: b = 0.0
      END TYPE

!-------------------------------------------------------------------------------
!>  Pointer to a gpkernel object. Used for creating arrays of gpkernel pointers.
!>  This is needed because fortran does not allow arrays of pointers directly.
!-------------------------------------------------------------------------------
      TYPE gpkernel_pointer
         TYPE (gpkernel_class), POINTER :: p => null()
      END TYPE

      PRIVATE oned_sqexp_kernel

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref pprofile_class.
!>
!>  Allocates memory and initializes a @ref pprofile_class.
!>
!>  @param[in] gp_type Kernel type discription.
!>  @param[in] b_arg  Array of kernel function coefficients
!>  @returns A pointer to a constructed @ref gpkernel_class object.
!-------------------------------------------------------------------------------
      FUNCTION gpkernel_construct(gp_type_arg, b_arg)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (gpkernel_class), POINTER :: gpkernel_construct
      CHARACTER (len=*), INTENT(in)         :: gp_type_arg
      REAL(rprec), DIMENSION(:), INTENT(in) :: b_arg

!  local variables
      CHARACTER (len = k_type_len)          :: gp_type_lc
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(gpkernel_construct)

! Store profile coeffients.
      gpkernel_construct%b = b_arg

!  Convert the type to lower case, and truncate
      gp_type_lc = gp_type_arg
      CALL tolower(gp_type_lc)
      SELECT CASE(TRIM(gp_type_lc))
      
         CASE ('sq_exp_1d')
            gpkernel_construct%k_type = gpkernel_oned_sqexp_type

         CASE ('none')
            gpkernel_construct%k_type = gpkernel_none_type

         CASE DEFAULT
            gpkernel_construct%k_type = gpkernel_none_type
            WRITE(*,*) 'Unrecognized k_type:', gp_type_lc
            WRITE(*,*) ' *** CHECK YOUR INPUT ***'
         
      END SELECT

      CALL profiler_set_stop_time('gpkernel_construct', start_time)

      END FUNCTION
!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref gpkernel_class object.
!>
!>  Deallocates memory and uninitializes a @ref gpkernel_class object.
!>
!>  @param[inout] this A @ref gpkernel_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE gpkernel_destruct(this)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (gpkernel_class), POINTER :: this

!  Start of executable code
      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Gets the value of a kernal at a radial position.
!>
!>  Evaluates the kernel using specifed radial positions. If the radial
!>  position is out side the plasma, it returns the baseline value. 
!>
!>  @param[in] this  A @ref gpkernel_class instance.
!>  @param[in] s1_arg First radial position used to evaluate the kernel.
!>  @param[in] s2_arg Second radial position used to evaluate the kernel.
!>  @returns The value of the kernel
!>  @todo This function makes no attempt to determine in the radial positions
!>  are in the domain. See pprofile_get_value
!-------------------------------------------------------------------------------
      FUNCTION gpkernel_get_value(this, s1_arg, s2_arg)

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec)                        :: gpkernel_get_value
      TYPE (gpkernel_class), INTENT (in) :: this
      REAL(rprec), INTENT(in)            :: s1_arg
      REAL(rprec), INTENT(in)            :: s2_arg

!  local variables
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      
      SELECT CASE(this%k_type)
      
         CASE (gpkernel_none_type)
            gpkernel_get_value = -1.E49_rprec

         CASE (gpkernel_oned_sqexp_type)
            gpkernel_get_value = zero
            gpkernel_get_value = oned_sqexp_kernel(s1_arg,s2_arg,this%b)
         
      END SELECT

! ECH TEMP
! Hard code 0 if s1, s2 >1 (FIX LATTER)
      IF ( (s1_arg > 1) .OR. (s2_arg > 1) ) THEN
         gpkernel_get_value = zero
      END IF
      CALL profiler_set_stop_time('gpkernel_get_value', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the name of the kernel type.
!>
!>  Convertes a kernel type to a string discription.
!>
!>  @param[in] this  A @ref gpkernel_class instance.
!>  @returns The name of the kernel type.
!-------------------------------------------------------------------------------
      FUNCTION gpkernel_get_k_type_name(this)

!  Declare Arguments
      CHARACTER(len=k_type_len)         :: gpkernel_get_k_type_name
      TYPE (gpkernel_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (this%k_type)

         CASE (gpkernel_oned_sqexp_type)
            gpkernel_get_k_type_name = 'sq_exp_1d'

         CASE DEFAULT
            gpkernel_get_k_type_name = 'none'

      END SELECT

      CALL profiler_set_stop_time('gpkernel_get_k_type_name',                  &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out the kernel to an output file.
!>
!>  This formats and writes out information about a kernel function.
!>
!>  @param[in] this A @ref gpkernel_class instance.
!>  @param[in] iou  Input/output unit of the output file.
!-------------------------------------------------------------------------------
      SUBROUTINE gpkernel_write(this, id, iou)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (gpkernel_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)     :: id
      INTEGER, INTENT(in)               :: iou

!  local variables
      INTEGER                           :: i
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Actual write. Could do a select-case on k_type.
      WRITE(iou,1100) id
      WRITE(iou,1200) TRIM(gpkernel_get_k_type_name(this))
      
      SELECT CASE(this%k_type)
         
         CASE (gpkernel_oned_sqexp_type)
            WRITE(iou,1210)
            WRITE(iou,1300) this%b(0:1)

      END SELECT

      CALL profiler_set_stop_time('gpkernel_write', start_time)

1100  FORMAT(/' Gaussian Process Kernel Write: id = ',a)
1200  FORMAT('   kernel_type = ',a)
1210  FORMAT(' b_0^2 * ( -( x1-x2)**2 / 2*b_1^2).    b(0:1) = ')
1300  FORMAT(2(2x,es15.8))

      END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Evaluate the one dimensional squared exponential kernel
!>
!>  b(0) * exp( - (x1- x2)**2 / 2*b(1) )
!>
!>  @param[in] x1 first position
!>  @param[in] x2 second position
!>  @param[in] b_array Array of coefficients
!>  @returns Value of the kernel
!-------------------------------------------------------------------------------

      FUNCTION oned_sqexp_kernel(x1,x2,b_array)

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec) :: oned_sqexp_kernel
      REAL(rprec), INTENT(in) :: x1
      REAL(rprec), INTENT(in) :: x2
      REAL(rprec), DIMENSION(ilb_gp:iub_gp), INTENT(in) :: b_array

!  local variables
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      oned_sqexp_kernel = b_array(0)**2*EXP(-1.0*(x1-x2)**2 /              &
     &                                   (2.0*b_array(1)**2))

      CALL profiler_set_stop_time('oned_sqexp_kernel', start_time)

      END FUNCTION

      END MODULE
