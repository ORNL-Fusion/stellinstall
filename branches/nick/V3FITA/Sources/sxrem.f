!*******************************************************************************
!  File magnetic.f
!  Contains module sxrem
!  Defines the base class of the type sxrem:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE sxrem

      USE stel_kinds, only: rprec
      USE integration_path

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) sxrem base class
!
!*******************************************************************************

      TYPE sxrem_class
         TYPE (path) :: chord_path
      END TYPE sxrem_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************

      PRIVATE :: sxr_function

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a sxrem.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_construct(chord_path)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_class), POINTER :: sxrem_construct
      TYPE (path), INTENT(in)     :: chord_path

!  Start of executable code
      ALLOCATE(sxrem_construct)

      sxrem_construct%chord_path = chord_path

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE sxrem_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      CALL path_destroy(this%chord_path)

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_modeled_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: sxrem_get_modeled_signal
      TYPE (sxrem_class), INTENT(inout)       :: this
      TYPE (model_class), INTENT(in)          :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

! local variables
      CHARACTER(len=1), ALLOCATABLE  :: context(:)
      INTEGER                        :: context_length

!  Start of executable code
      sxrem_get_modeled_signal = 0.0

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
      context_length = SIZE(TRANSFER(a_model, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(a_model, context)
      sxrem_get_modeled_signal(1) = path_integrate(this%chord_Path,            &
     &                                             sxr_function,               &
     &                                             context)

      sigma = 0.0
      DEALLOCATE(context)

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  SXRem Integration Function
!
!  lineAve = Int[sxrem*dl]
!
!-------------------------------------------------------------------------------
      FUNCTION sxr_function(context, xcart, dxcart, dx)
      USE model
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: sxr_function

! local variables
      TYPE(model_class)                     :: a_model

!  Start of executable code
      a_model = TRANSFER(context, a_model)
      sxr_function = model_get_sxrem(a_model, cart_to_cyl(xcart))*dx

      END FUNCTION

      END MODULE
