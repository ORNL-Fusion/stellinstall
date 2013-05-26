!*******************************************************************************
!  File model.f
!  Contains module model
!  Defines the base class of the type model:
!
!  Super Class: None
!  Sub Classes: None
!
!*******************************************************************************

      MODULE model
      USE equilibrium

      IMPLICIT NONE

!*******************************************************************************
!  model module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: model_none_type     = -1

      INTEGER, PARAMETER :: model_ne_type       =  0 ! Straight ne profile.
      INTEGER, PARAMETER :: model_ne_te_p_type  =  1 ! derived from temp/press

      INTEGER, PARAMETER :: model_sxrem_type    =  0 ! Straight sxrem profile.

      INTEGER, PARAMETER :: model_te_type       =  0 ! Straight te profile.
      INTEGER, PARAMETER :: model_te_ne_p_type  =  1 ! derived from dens/press

!  Model parameter ids
      INTEGER, PARAMETER :: model_ne_unit_id           = 0
      INTEGER, PARAMETER :: model_ne_min_id            = 1
      INTEGER, PARAMETER :: model_pressure_fraction_id = 2

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) model class
!
!*******************************************************************************

      TYPE model_class
         INTEGER      :: ne_type = model_none_type
         INTEGER      :: sxrem_type = model_none_type
         INTEGER      :: te_type = model_none_type
         REAL (rprec) :: ne_unit = 0.0
         REAL (rprec) :: ne_min = 0.0
         REAL (rprec) :: pressure_fraction = 0.0
         TYPE (equilibrium_class), POINTER :: equilibrium => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION model_construct(ne_type, sxrem_type, te_type,                   &
     &                         ne_unit, ne_min, pressure_fraction,             &
     &                         equilibrium)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), POINTER                   :: model_construct
      CHARACTER (len=data_name_length), INTENT(in)  :: ne_type
      CHARACTER (len=data_name_length), INTENT(in)  :: sxrem_type
      CHARACTER (len=data_name_length), INTENT(in)  :: te_type
      REAL (rprec), INTENT(in)                      :: ne_unit
      REAL (rprec), INTENT(in)                      :: ne_min
      REAL (rprec), INTENT(in)                      :: pressure_fraction
      TYPE (equilibrium_class), POINTER, INTENT(in) :: equilibrium

!  Start of executable code
      ALLOCATE(model_construct)

      SELECT CASE (TRIM(ADJUSTL(ne_type)))

         CASE ('none')
            model_construct%ne_type = model_none_type

         CASE ('pp_ne')
            model_construct%ne_type = model_ne_type

         CASE ('pp_te_vmec_p')
            model_construct%ne_type = model_ne_te_p_type

         CASE DEFAULT
            WRITE (*,*) 'ne type: ' // TRIM(ADJUSTL(ne_type)) //               &
     &                  'setting type to none.'
            model_construct%ne_type = model_ne_type

      END SELECT

      SELECT CASE (TRIM(ADJUSTL(sxrem_type)))

         CASE ('none')
            model_construct%sxrem_type = model_none_type

         CASE ('pp_sxrem')
            model_construct%ne_type = model_sxrem_type

         CASE DEFAULT
            WRITE (*,*) 'sxrem type: ' // TRIM(ADJUSTL(ne_type)) //            &
     &                  'setting type to none.'
            model_construct%sxrem_type = model_none_type

      END SELECT

      SELECT CASE (TRIM(ADJUSTL(te_type)))

         CASE ('none')
            model_construct%ne_type = model_none_type

         CASE ('pp_te')
            model_construct%te_type = model_te_type

         CASE ('pp_ne_vmec_p')
            model_construct%te_type = model_te_ne_p_type

         CASE DEFAULT
            WRITE (*,*) 'te type: ' // TRIM(ADJUSTL(te_type)) //               &
     &                  'setting type to none.'
            model_construct%ne_type = model_none_type

      END SELECT

      model_construct%ne_unit = ne_unit
      model_construct%ne_min = ne_min
      model_construct%pressure_fraction = pressure_fraction
      model_construct%equilibrium => equilibrium

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE model_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%equilibrium)) THEN
         CALL equilibrium_destruct(this%equilibrium)
         this%equilibrium => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

      END MODULE