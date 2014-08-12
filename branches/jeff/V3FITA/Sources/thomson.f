!*******************************************************************************
!>  @file thomson.f
!>  @brief Contains module @ref thomson.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Implements thomson scattering diagnostic. Defines the base class of the type
!>  @ref thomson_class.
!>  @par Super Class:
!>  @ref diagnostic
!*******************************************************************************

      MODULE thomson

      USE stel_kinds, only: rprec

      IMPLICIT NONE
!*******************************************************************************
!  diagnostic module parameters
!*******************************************************************************
!>  Type descriptor for thomson no type.
      INTEGER, PARAMETER :: thomson_no_type = -1
!>  Type descriptor for thomson type measureing temperature.
      INTEGER, PARAMETER :: thomson_temperature_type = 0
!>  Type descriptor for intpol type measureing denisty.
      INTEGER, PARAMETER :: thomson_density_type = 1
!>  Type descriptor for intpol type measureing pressure.
      INTEGER, PARAMETER :: thomson_pressure_type = 2

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) thomson base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a thomson scattering signal.
!>  @par Super Class:
!>  @ref diagnostic
!-------------------------------------------------------------------------------
      TYPE thomson_class
!>  Type descirptor of the thomson type.
!>  @par Possible values are:
!>  * @ref thomson_no_type
!>  * @ref thomson_temperature_type
!>  * @ref thomson_density_type
!>  * @ref thomson_pressure_type
         INTEGER                    :: type
         REAL (rprec), DIMENSION(3) :: xcart
      END TYPE thomson_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref thomson_class object measureing temperature.
!>
!>  @param[in] xcart The measurement point.
!>  @returns A pointer to a constructed @ref thomson_class object.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct_te(xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct_te
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      thomson_construct_te =>                                                  &
     &   thomson_construct(thomson_temperature_type, xcart)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref thomson_class object measureing density.
!>
!>  @param[in] xcart The measurement point.
!>  @returns A pointer to a constructed @ref thomson_class object.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct_ne(xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct_ne
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      thomson_construct_ne => thomson_construct(thomson_density_type,          &
     &                                          xcart)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref thomson_class object measureing pressure.
!>
!>  @param[in] xcart The measurement point.
!>  @returns A pointer to a constructed @ref thomson_class object.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct_p(xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct_p
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      thomson_construct_p => thomson_construct(thomson_pressure_type,          &
     &                                         xcart)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref thomson_class object.
!>
!>  Allocates memory and initializes a @ref thomson_class object.
!>
!>  @param[in] type  The measurement type. Possible value are
!>                   @ref thomson_no_type, @ref thomson_temperature_type,
!>                   @ref thomson_density_type and @ref thomson_pressure_type.
!>  @param[in] xcart The measurement point.
!>  @returns A pointer to a constructed @ref thomson_class object.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct(type, xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct
      INTEGER, INTENT(in)                    :: type
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      ALLOCATE(thomson_construct)

      thomson_construct%type = type
      thomson_construct%xcart = xcart

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref intpol_class object.
!>
!>  Deallocates memory and uninitializes a @ref intpol_class object.
!>
!>  @param[inout] this A @ref intpol_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE thomson_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER :: this

!  Start of executable code
      this%type = thomson_no_type
      this%xcart = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  Calculates the modeled signal at a point. Temperature is provided by
!>  @ref model::model_get_te. Denisty is provided by @ref model::model_get_ne.
!>  Pressure is provided by @ref equilibrium::equilibrium_get_p.
!>
!>  @param[in]    this    A @ref thomson_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION thomson_get_modeled_signal(this, a_model, sigma)
      USE model
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: thomson_get_modeled_signal
      TYPE (thomson_class), INTENT(in)        :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Local variables

!  Start of executable code
      thomson_get_modeled_signal = 0.0

      SELECT CASE(this%type)

         CASE (thomson_temperature_type)
            thomson_get_modeled_signal(1) =                                    &
     &         model_get_te(a_model, this%xcart)

         CASE (thomson_density_type)
            thomson_get_modeled_signal(1) =                                    &
     &         model_get_ne(a_model, this%xcart)

         CASE (thomson_pressure_type)
            thomson_get_modeled_signal(1) =                                    &
     &         equilibrium_get_p(a_model%equilibrium, this%xcart)

      END SELECT

      sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the thomson type.
!>
!>  Returns a description of the thomson type for use when writting output files.
!>
!>  @param[in] this A @ref thomson_class instance.
!>  @returns A string describing the thomson type.
!-------------------------------------------------------------------------------
      FUNCTION thomson_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: thomson_get_signal_type
      TYPE (thomson_class), INTENT(in) :: this

!  Start of executable code
      thomson_get_signal_type = 'thscte'

      END FUNCTION

      END MODULE
