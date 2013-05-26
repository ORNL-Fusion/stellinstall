!*******************************************************************************
!  File equilibrium.f
!  Contains module equilibrium
!  Defines the base class of the type equilibrium:
!
!  Super Class: None
!  Sub Classes: vmec
!
!*******************************************************************************

      MODULE equilibrium
      USE vmec_equilibrium

      IMPLICIT NONE

!*******************************************************************************
!  equilibrium module parameters
!*******************************************************************************

      INTEGER, PARAMETER :: equilibrium_no_type = -1
      INTEGER, PARAMETER :: equilibrium_vmec_type = 0

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) equilibrium base class
!
!*******************************************************************************

      TYPE equilibrium_class
         INTEGER                    :: type = equilibrium_no_type
         TYPE (vmec_class), POINTER :: vmec => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE equilibrium_construct
         MODULE PROCEDURE equilibrium_construct_vmec
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct an signal equilibrium a vmec equilibrium.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_construct_vmec(vmec_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), POINTER :: equilibrium_construct_vmec
      TYPE (vmec_class), POINTER, INTENT(in) :: vmec_object

!  Start of executable code
      ALLOCATE(equilibrium_construct_vmec)

      equilibrium_construct_vmec%type = equilibrium_vmec_type
      equilibrium_construct_vmec%vmec => vmec_object

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE equilibrium_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%vmec)) THEN
         CALL vmec_destruct(this%vmec)
         this%vmec => null()
      END IF

      this%type = equilibrium_no_type

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the id number for a equilibrium parameter
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_param_id(this, param_name)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: equilibrium_get_param_id
      TYPE (equilibrium_class), INTENT(in)         :: this
      CHARACTER (len=data_name_length), INTENT(in) :: param_name

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_param_id =                                         &
     &         vmec_get_param_id(this%vmec, param_name)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!   Get Scaler Real value for an id number
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_scaler_value(this, id)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_scaler_value
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_scaler_value =                                     &
     &         vmec_get_scaler_value(this%vmec, id)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get equilibrium parameter name for an id number
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_param_name(this, id)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=data_name_length) :: equilibrium_get_param_name
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_param_name =                                       &
     &         vmec_get_param_name(this%vmec, id)

      END SELECT

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Query variable id for scaler type
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_scaler_value(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_scaler_value
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_scaler_value =                                      &
     &         vmec_is_scaler_value(this%vmec, id)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 1D array
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_1d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_1d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_1d_array = vmec_is_1d_array(this%vmec, id)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 2D array
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_2d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_2d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_2d_array = vmec_is_2d_array(this%vmec, id)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 3D array
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_3d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_3d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_3d_array = vmec_is_3d_array(this%vmec, id)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 4D array
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_4d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_4d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

      CASE (equilibrium_vmec_type)
         equilibrium_is_4d_array = vmec_is_4d_array(this%vmec, id)

      END SELECT

      END FUNCTION

      END MODULE
