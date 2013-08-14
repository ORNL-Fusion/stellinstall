!-------------------------------------------------------------------------------
!  The @header, @table_section, @table_subsection, @item and @end_table commands
!  are custom defined commands in Doxygen.in. They are defined under ALIASES.
!  For the page created here, the 80 column limit is exceeded. Arguments of
!  aliases are separated by ','. You intended ',' to be a string you must use an
!  escaped comma '\,'.
!
!>  @page model_sec Model Manual
!>
!>  @tableofcontents
!>  @section model_intro_sec Introduction
!>  This page documents the V3FIT interface to the model. This lists all the
!>  parameters associated with a model. Model parameters are common for all
!>  equilibrium types. All parameters are documented in a table of the following
!>  form.
!>  @header{Input variable, Description, Code Reference}
!>
!>  @section model_recon_param_sec Model Reconstruction Parameters
!>  @begin_table
!>     @item{ne_unit,           Electron density scaling.,                                @ref model::model_class::ne_unit}
!>     @item{ne_min,            Minimum electron density.,                                @ref model::model_class::ne_min}
!>     @item{pressure_fraction, Fractional portion of the pressure due to the electrons., @ref model::model_class::pressure_fraction}
!>  @end_table
!*******************************************************************************
!>  @file model.f
!>  @brief Contains module @ref model.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref model_class. The model contains
!>  information not specific to the equilibrium.
!*******************************************************************************

      MODULE model
      USE equilibrium

      IMPLICIT NONE

!*******************************************************************************
!  model module parameters
!*******************************************************************************
!>  No model.
      INTEGER, PARAMETER, PRIVATE :: model_none_type     = -1

!>  Electron denisty model.
      INTEGER, PARAMETER, PRIVATE :: model_ne_type       =  0
!>  Electron denisty model is derived from the temperature and pressure.
      INTEGER, PARAMETER, PRIVATE :: model_ne_te_p_type  =  1

!>  Soft X-ray Emissivity.
      INTEGER, PARAMETER, PRIVATE :: model_sxrem_type    =  0

!>  Electron temperature model.
      INTEGER, PARAMETER, PRIVATE :: model_te_type       =  0
!>  Electron temperature model is derived from the density and pressure.
      INTEGER, PARAMETER, PRIVATE :: model_te_ne_p_type  =  1

!  Model parameter ids
!>  Parameter id for the electron density units.
      INTEGER, PARAMETER :: model_ne_unit_id           = 0
!>  Parameter id for the minimum electron density.
      INTEGER, PARAMETER :: model_ne_min_id            = 1
!>  Parameter id for the electrion fraction of the pressure.
      INTEGER, PARAMETER :: model_pressure_fraction_id = 2

!>  Conversion from Joules to eV.
      REAL(rprec), PARAMETER :: eV_per_Joule = one/1.602e-19_rprec

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) model class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a model.
!-------------------------------------------------------------------------------
      TYPE model_class
!>  Type decritption of the electron denisty. Possible values are
!>  @ref model_none_type, @ref model_ne_type and @ref model_ne_te_p_type.
         INTEGER      :: ne_type = model_none_type
!>  Type decritption of the soft x-ray emissivity. Possible values are
!>  @ref model_none_type and @ref model_sxrem_type.
         INTEGER      :: sxrem_type = model_none_type
!>  Type decritption of the electron temperature. Possible values are
!>  @ref model_none_type, @ref model_te_type and @ref model_te_ne_p_type.
         INTEGER      :: te_type = model_none_type
!>  Unit scaling of the electron density.
         REAL (rprec) :: ne_unit = 0.0
!>  Minimum electron density.
         REAL (rprec) :: ne_min = 0.0
!>  Electron pressure fraction.
         REAL (rprec) :: pressure_fraction = 0.0
!>  The @ref equilibrium.
         TYPE (equilibrium_class), POINTER :: equilibrium => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref model_class containing an @ref equilibrium object.
!>
!>  Allocates memory and initializes a @ref model_class object.
!>  @param[in] ne_type           Model type for the electron density.
!>  @param[in] sxrem_type        Model type for the soft x-ray emissivity.
!>  @param[in] te_type           Model type for the electron temperature.
!>  @param[in] ne_unit           Scaleing of the electron denisty.
!>  @param[in] ne_min            Minimum electron density.
!>  @param[in] pressure_fraction Fraction of the pressure provided by electrons.
!>  @param[in] equilibrium       An instance of an @ref equilibrium.
!>  @returns A pointer to a constructed @ref model_class object.
!-------------------------------------------------------------------------------
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

         CASE ('pp_te_p')
            model_construct%ne_type = model_ne_te_p_type

         CASE DEFAULT
            WRITE (*,*) 'ne type: ' // TRIM(ADJUSTL(ne_type)) //               &
     &                  'setting type to none.'
            model_construct%ne_type = model_none_type

      END SELECT

      SELECT CASE (TRIM(ADJUSTL(sxrem_type)))

         CASE ('none')
            model_construct%sxrem_type = model_none_type

         CASE ('pp_sxrem')
            model_construct%sxrem_type = model_sxrem_type

         CASE DEFAULT
            WRITE (*,*) 'sxrem type: ' // TRIM(ADJUSTL(ne_type)) //            &
     &                  'setting type to none.'
            model_construct%sxrem_type = model_none_type

      END SELECT

      SELECT CASE (TRIM(ADJUSTL(te_type)))

         CASE ('none')
            model_construct%te_type = model_none_type

         CASE ('pp_te')
            model_construct%te_type = model_te_type

!>  @todo FIXME: Name here shouldn't have vmec in it
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
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref model_class object.
!>
!>  Deallocates memory and uninitializes a @ref model_class object.
!>
!>  @param[inout] this A @ref model_class instance.
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
!>  @brief Sets the value of a reconstruction model parameter.
!>
!>  Sets a model reconstruction parameter. If the id is not a model parameter,
!>  it is defered to the @ref equilibrium.
!>
!>  @param[inout] this    A @ref model_class instance.
!>  @param[in]    id      ID of the parameter.
!>  @param[in]    i_index The ith index of the parameter.
!>  @param[in]    j_index The jth index of the parameter.
!>  @param[in]    value   The value of the parameter.
!>  @note Model parameter do not use the indicies. However an equilibrium might.
!-------------------------------------------------------------------------------
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
!>  @brief Get the id for a reconstruction parameter.
!>
!>  Gets a model reconstruction parameter. If the id is not a model parameter,
!>  it is defered to the @ref equilibrium.
!>
!>  @param[in] this       A @ref model_class instance.
!>  @param[in] param_name Name of a reconstruction parameter.
!>  @returns The id for a reconstruction parameter.
!-------------------------------------------------------------------------------
      FUNCTION model_get_param_id(this, param_name)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                        :: model_get_param_id
      TYPE (model_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)  :: param_name

!  Start of executable code
      SELECT CASE (TRIM(param_name))

         CASE ('ne_unit')
            model_get_param_id = model_ne_unit_id

         CASE ('ne_min')
            model_get_param_id = model_ne_min_id

         CASE ('pressure_fraction')
            model_get_param_id = model_pressure_fraction_id

         CASE DEFAULT
            model_get_param_id =                                               &
     &         equilibrium_get_param_id(this%equilibrium, param_name)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the value of a reconstruction model parameter.
!>
!>  Gets a model reconstruction parameter value. If the id is not a model
!>  parameter, it is defered to the @ref equilibrium.
!>
!>  @param[in] this    A @ref model_class instance.
!>  @param[in] id      ID of the parameter.
!>  @param[in] i_index The ith index of the parameter.
!>  @param[in] j_index The jth index of the parameter.
!>  @returns The value of the parameter.
!>  @note Model parameter do not use the indicies. However an equilibrium might.
!-------------------------------------------------------------------------------
      FUNCTION model_get_param_value(this, id, i_index, j_index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: model_get_param_value
      TYPE (model_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: id
      INTEGER, INTENT(in)            :: i_index
      INTEGER, INTENT(in)            :: j_index

!  Start of executable code
      SELECT CASE (id)

         CASE (model_ne_unit_id)
            model_get_param_value = this%ne_unit

         CASE (model_ne_min_id)
            model_get_param_value = this%ne_min

         CASE (model_pressure_fraction_id)
            model_get_param_value = this%pressure_fraction

         CASE DEFAULT
            model_get_param_value =                                            &
     &         equilibrium_get_param_value(this%equilibrium, id,               &
     &                                     i_index, j_index)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the name of a reconstruction model parameter.
!>
!>  Gets a model reconstruction parameter name. If the id is not a model
!>  parameter, it is defered to the @ref equilibrium.
!>
!>  @param[in] this A @ref model_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns The name of the parameter.
!-------------------------------------------------------------------------------
      FUNCTION model_get_param_name(this, id)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=data_name_length) :: model_get_param_name
      TYPE (model_class), INTENT(in)  :: this
      INTEGER, INTENT(in)             :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (model_ne_unit_id)
            model_get_param_name = 'ne_unit'

         CASE (model_ne_min_id)
            model_get_param_name = 'ne_min'

         CASE (model_pressure_fraction_id)
            model_get_param_name = 'pressure_fraction'

         CASE DEFAULT
            model_get_param_name =                                             &
     &         equilibrium_get_param_name(this%equilibrium, id)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron density at a position.
!>
!>  Gets the electron density at a position. Electron density is computed based
!>  on the type of @ref model_class::ne_type.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @param[in] x_cart Cartesian position to get the electron density at.
!>  @returns The electron density at x_cart.
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
            model_get_ne = MAX(this%ne_min, this%ne_unit*model_get_ne)

         CASE (model_ne_te_p_type)
            model_get_ne = eV_per_Joule*this%pressure_fraction                 &
     &                   * equilibrium_get_p(this%equilibrium, x_cart)         &
     &                   / equilibrium_get_te(this%equilibrium, x_cart)

         CASE DEFAULT
             model_get_ne = -2.0E49_rprec

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron temperature at a position.
!>
!>  Gets the electron temperature at a position. Electron temperature is
!>  computed based on the type of @ref model_class::te_type.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @param[in] x_cart Cartesian position to get the electron temperature at.
!>  @returns The electron temperature at x_cart.
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
!>  @brief Gets the soft x-ray emissivity at a position.
!>
!>  Gets the soft x-ray emissivity at a position. Soft x-ray emissivity is
!>  computed based on the type of @ref model_class::sxrem_type.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @param[in] x_cart Cartesian position to get the soft x-ray emissivity at.
!>  @returns The soft x-ray emissivity at x_cart.
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

!-------------------------------------------------------------------------------
!>  @brief Gets the electron density type as a string.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @returns A string description of the electron density type.
!-------------------------------------------------------------------------------
      FUNCTION model_get_ne_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: model_get_ne_type
      TYPE (model_class), INTENT(in)   :: this

!  Start of executable code
      SELECT CASE (this%ne_type)

         CASE (model_ne_type)
            model_get_ne_type = 'pp_ne'

         CASE (model_ne_te_p_type)
            model_get_ne_type = 'pp_te_p'

         CASE DEFAULT
            model_get_ne_type = 'none'

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron temperature type as a string.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @returns A string description of the electron temperature type.
!-------------------------------------------------------------------------------
      FUNCTION model_get_te_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: model_get_te_type
      TYPE (model_class), INTENT(in)   :: this

!  Start of executable code
      SELECT CASE (this%te_type)

         CASE (model_te_type)
            model_get_te_type = 'pp_te'

         CASE (model_ne_te_p_type)
            model_get_te_type = 'pp_ne_vmec_p'

         CASE DEFAULT
            model_get_te_type = 'none'

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the soft x-ray emissivity type as a string.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @returns A string description of the soft x-ray emissivity type.
!-------------------------------------------------------------------------------
      FUNCTION model_get_sxrem_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: model_get_sxrem_type
      TYPE (model_class), INTENT(in)   :: this

!  Start of executable code
      SELECT CASE (this%sxrem_type)

      CASE (model_sxrem_type)
         model_get_sxrem_type = 'pp_sxrem'

      CASE DEFAULT
         model_get_sxrem_type = 'none'

      END SELECT

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Determines if a parameter id is a reconstruction parameter.
!>
!>  Determines if a parameter id is a reconstruction parameter. If the parameter
!>  is not a model parameter, it is defered to the @ref equilibrium.
!>
!>  @param[in] this A @ref model_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the id is a reconstruction parameter and false if not.
!-------------------------------------------------------------------------------
      FUNCTION model_is_recon_param(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                        :: model_is_recon_param
      TYPE (model_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (model_ne_unit_id, model_ne_min_id,                              &
     &         model_pressure_fraction_id)
            model_is_recon_param = .true.

         CASE DEFAULT
            model_is_recon_param =                                             &
     &         equilibrium_is_recon_param(this%equilibrium, id)

      END SELECT

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out the model to an output file.
!>
!>  Write out information about the model to an output file. This also then
!>  calls the @ref equilibrium_write subroutine to write out the equilibrium
!>  information aswell.
!>
!>  @param[in] this A @ref model_class instance.
!>  @param[in] iou  Input/output unit of the output file.
!-------------------------------------------------------------------------------
      SUBROUTINE model_write(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: iou

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,*) ' *** Model Parameters'
      WRITE (iou,*) 'model_ne_type is ',                                       &
     &              TRIM(model_get_ne_type(this))
      WRITE (iou,*) 'model_sxrem_type is ',                                    &
     &              TRIM(model_get_sxrem_type(this))
      WRITE (iou,*) 'model_te_type is ',                                       &
     &              TRIM(model_get_te_type(this))
      WRITE (iou, 1000) 'ne_pp_unit is ', this%ne_unit
      WRITE (iou, 1000) 'ne_min is ', this%ne_min
      WRITE (iou, 1000) 'e_pressure_fraction is ',                             &
     &                  this%pressure_fraction
1000  FORMAT(1x,a,2x,es12.5)

      CALL equilibrium_write(this%equilibrium, iou)

      END SUBROUTINE

      END MODULE