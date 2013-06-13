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
      INTEGER, PARAMETER, PRIVATE :: model_none_type     = -1

      INTEGER, PARAMETER, PRIVATE :: model_ne_type       =  0 ! Straight ne profile.
      INTEGER, PARAMETER, PRIVATE :: model_ne_te_p_type  =  1 ! derived from temp/press

      INTEGER, PARAMETER, PRIVATE :: model_sxrem_type    =  0 ! Straight sxrem profile.

      INTEGER, PARAMETER, PRIVATE :: model_te_type       =  0 ! Straight te profile.
      INTEGER, PARAMETER, PRIVATE :: model_te_ne_p_type  =  1 ! derived from dens/press

!  Model parameter ids
      INTEGER, PARAMETER :: model_ne_unit_id           = 0
      INTEGER, PARAMETER :: model_ne_min_id            = 1
      INTEGER, PARAMETER :: model_pressure_fraction_id = 2

      REAL(rprec), PARAMETER :: eV_per_Joule = one/1.602e-19_rprec

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

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
      SUBROUTINE model_set_param(this, id, i_index, j_index, value)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), INTENT(inout) :: this
      INTEGER, INTENT(in)               :: id
      INTEGER, INTENT(in)               :: i_index
      INTEGER, INTENT(in)               :: j_index
      REAL (rprec), INTENT(in)          :: value

!  Start of executable code
      SELECT CASE (id)

         CASE (model_ne_unit_id)
            this%ne_unit = value

         CASE (model_ne_min_id)
            this%ne_min = value

         CASE (model_pressure_fraction_id)
            this%pressure_fraction = value

         CASE DEFAULT
            CALL equilibrium_set_param(this%equilibrium, id,                   &
     &                                 i_index, j_index, value)

      END SELECT

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get model electron density at position
!-------------------------------------------------------------------------------
      FUNCTION model_get_ne(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: model_get_ne
      TYPE (model_class), INTENT(in)         :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%ne_type)

         CASE (model_ne_type)
            model_get_ne = equilibrium_get_ne(this%equilibrium, x_cart)

         CASE (model_ne_te_p_type)
            model_get_ne = eV_per_Joule*this%pressure_fraction                 &
     &                   * equilibrium_get_p(this%equilibrium, x_cart)         &
     &                   / equilibrium_get_te(this%equilibrium, x_cart)

         CASE DEFAULT
             model_get_ne = -2.0E49_rprec

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get model electron temperature at position
!-------------------------------------------------------------------------------
      FUNCTION model_get_te(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: model_get_te
      TYPE (model_class), INTENT(in)         :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%ne_type)

         CASE (model_te_type)
            model_get_te = equilibrium_get_te(this%equilibrium, x_cart)

         CASE (model_te_ne_p_type)
            model_get_te = eV_per_Joule*this%pressure_fraction                 &
     &                   * equilibrium_get_p(this%equilibrium, x_cart)         &
     &                   / equilibrium_get_ne(this%equilibrium, x_cart)

         CASE DEFAULT
            model_get_te = -2.0E49_rprec

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get model soft x-ray emissivity at position
!-------------------------------------------------------------------------------
      FUNCTION model_get_sxrem(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: model_get_sxrem
      TYPE (model_class), INTENT(in)         :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%sxrem_type)

         CASE (model_sxrem_type)
            model_get_sxrem =                                                  &
     &         equilibrium_get_sxrem(this%equilibrium, x_cart)

         CASE DEFAULT
            model_get_sxrem = -2.0E49_rprec

      END SELECT

      END FUNCTION

      END MODULE