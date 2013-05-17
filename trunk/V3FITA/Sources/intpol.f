!*******************************************************************************
!  File magnetic.f
!  Contains module intpol
!  Defines the base class of the type intpol:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE intpol

      USE stel_kinds, only: rprec, dp
      USE integration_path
      USE coordinate_utilities

      IMPLICIT NONE

!*******************************************************************************
!  diagnostic module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: intpol_no_type = -1
      INTEGER, PARAMETER :: intpol_interferometry_type = 0
      INTEGER, PARAMETER :: intpol_polarimetry_type = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) intpol base class
!
!*******************************************************************************

      TYPE intpol_class
         INTEGER      :: type
         REAL (rprec) :: wavelength
         LOGICAL      :: in_degrees
         TYPE (path)  :: chord_path
      END TYPE intpol_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE intpol_construct
         MODULE PROCEDURE intpol_construct_int,                                &
     &                    intpol_construct_pol
      END INTERFACE

      PRIVATE :: int_function, pol_function

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a interferometer diagnostic
!-------------------------------------------------------------------------------
      FUNCTION intpol_construct_int(chord_path)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class)     :: intpol_construct_int
      TYPE (path), INTENT(in) :: chord_path

!  Start of executable code
      intpol_construct_int%type = intpol_interferometry_type
      intpol_construct_int%wavelength = 0.0
      intpol_construct_int%in_degrees = .false.
      intpol_construct_int%chord_path = chord_path

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a polarimetry diagnostic
!-------------------------------------------------------------------------------
      FUNCTION intpol_construct_pol(wavelength, in_degrees, chord_path)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class)      :: intpol_construct_pol
      REAL (rprec), INTENT(in) :: wavelength
      LOGICAL, INTENT(in)      :: in_degrees
      TYPE (path), INTENT(in)  :: chord_path

!  Start of executable code
      intpol_construct_pol%type = intpol_polarimetry_type
      intpol_construct_pol%wavelength = wavelength
      intpol_construct_pol%in_degrees = in_degrees
      intpol_construct_pol%chord_path = chord_path

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE intpol_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), INTENT(inout) :: this

!  Start of executable code
      this%type = intpol_no_type
      this%wavelength = 0.0
      this%in_degrees = .false.
      CALL path_destroy(this%chord_path)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2)         :: intpol_get_modeled_signal
      TYPE (intpol_class), INTENT(inout) :: this
      TYPE (model), INTENT(inout)        :: a_model

!  local variables
      CHARACTER(len=1), ALLOCATABLE :: context(:)
      INTEGER                       :: context_length

!  local parameters
      REAL(rprec), PARAMETER :: polar_constant = 2.62E-13_dp

!  Start of executable code
      IF (.not. a_model%eqstate%l_def_aux1) THEN
         CALL eq_aux_1_define(a_model%eqstate)
      END IF

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
      context_length = SIZE(TRANSFER(a_model, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(a_model, context)

      SELECT CASE (this%type)

         CASE (intpol_interferometry_type)
            intpol_get_modeled_signal(1) =                                     &
     &         path_integrate(this%chord_path, int_function, context)

         CASE (intpol_polarimetry_type)
            intpol_get_modeled_signal(1) = polar_constant *                    &
     &         (this%wavelength**2.0_dp) *                                     &
     &         path_integrate(this%chord_path, pol_function, context)
            IF (this%in_degrees) THEN
                intpol_get_modeled_signal(1) =                                 &
     &             intpol_get_modeled_signal(1)/degree
            ENDIF

      END SELECT

      DEALLOCATE(context)

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Interferometer Integration Function
!
!  lineAve = Int[ne*dl]
!
!-------------------------------------------------------------------------------
      FUNCTION int_function(context, xcart, dxcart, dx)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: int_function

!  local variables
      TYPE(model) :: a_model

!  Start of executable code
      a_model = TRANSFER(context, a_model)
      int_function = model_get_ne_xcyl(a_model, cart_to_cyl(xcart))*dx

      END FUNCTION

!-------------------------------------------------------------------------------
!  Polarimetry Integration Function
!
!  Phi = 2.62E-13*lamba^2*Int[ne*B.dl]
!
!-------------------------------------------------------------------------------
      FUNCTION pol_function(context, xcart, dxcart, dx)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: pol_function

!  local variables
      TYPE(model)               :: a_model
      REAL(rprec), DIMENSION(3) :: r_cyl, bcart

!  Start of executable code
      a_model = TRANSFER(context, a_model)
      r_cyl = cart_to_cyl(xcart)
      CALL eq_get_B_field(a_model%eqstate, r_cyl, bcart)
      pol_function = model_get_ne_xcyl(a_model, r_cyl)                         &
     &             * DOT_PRODUCT(bcart, dxcart)

      END FUNCTION

      END MODULE
