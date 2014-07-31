!-------------------------------------------------------------------------------
!  The @header, @table_section, @table_subsection, @item and @end_table commands
!  are custom defined commands in Doxygen.in. They are defined under ALIASES.
!  For the page created here, the 80 column limit is exceeded. Arguments of
!  aliases are separated by ','. If you intended ',' to be a string you must use
!  an escaped comma '\,'.
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
!>
!>  @section model_derived_param_sec Model Derived Parameters
!>  The following derived parameters are defined on a grid. The size of that
!>  grid will depend on equilibrium used.
!>  @begin_table
!>     @item{ne_grid,    Electron density grid.,      @ref model::model_class::ne_grid}
!>     @item{te_grid,    Electron temperature grid.,  @ref model::model_class::te_grid}
!>     @item{sxrem_grid, Soft X-ray emissivity grid., @ref model::model_class::sxrem_grid}
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

!  Derived Parameters
!>  Parameter id for the electrion density grid.
      INTEGER, PARAMETER :: model_ne_grid_id           = 3
!>  Parameter id for the electrion temperature grid.
      INTEGER, PARAMETER :: model_te_grid_id           = 4
!>  Parameter id for the soft x-ray emissivity grid.
      INTEGER, PARAMETER :: model_sxrem_grid_id        = 5

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
         INTEGER, DIMENSION(:), POINTER :: sxrem_type => null()
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

!>  Grid start
         REAL (rprec)                          :: grid_start
!>  Grid stop
         REAL (rprec)                          :: grid_step
!>  Girded electron density profile.
         REAL (rprec), DIMENSION(:), POINTER   :: ne_grid => null()
!>  Grided electron soft x-ray profile.
         REAL (rprec), DIMENSION(:,:), POINTER :: sxrem_grid => null()
!>  Grided electron temperature profile.
         REAL (rprec), DIMENSION(:), POINTER   :: te_grid => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the model density profile values.
!-------------------------------------------------------------------------------
      INTERFACE model_get_ne
         MODULE PROCEDURE model_get_ne_cart,                                   &
     &                    model_get_ne_radial
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for the model temperature profile values.
!-------------------------------------------------------------------------------
      INTERFACE model_get_te
         MODULE PROCEDURE model_get_te_cart,                                   &
     &                    model_get_te_radial
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for the model soft x-ray emissivity profile values.
!-------------------------------------------------------------------------------
      INTERFACE model_get_sxrem
         MODULE PROCEDURE model_get_sxrem_cart,                                &
     &                    model_get_sxrem_radial
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref model_class types using
!>  @ref model_construct_new or @ref signal_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE model_construct
         MODULE PROCEDURE model_construct_new,                                 &
     &                    model_construct_mpi
      END INTERFACE

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
      FUNCTION model_construct_new(ne_type, sxrem_type, te_type,               &
     &                             ne_unit, ne_min, pressure_fraction,         &
     &                             equilibrium)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), POINTER :: model_construct_new
      CHARACTER (len=data_name_length), INTENT(in)  :: ne_type
      CHARACTER (len=data_name_length), DIMENSION(:), INTENT(in) ::            &
     &   sxrem_type
      CHARACTER (len=data_name_length), INTENT(in)  :: te_type
      REAL (rprec), INTENT(in)                      :: ne_unit
      REAL (rprec), INTENT(in)                      :: ne_min
      REAL (rprec), INTENT(in)                      :: pressure_fraction
      TYPE (equilibrium_class), POINTER             :: equilibrium

!  local variables
      INTEGER                                       :: i, grid_size
      REAL (rprec)                                  :: grid_stop

!  Start of executable code
      ALLOCATE(model_construct_new)

      SELECT CASE (TRIM(ADJUSTL(ne_type)))

         CASE ('none')
            model_construct_new%ne_type = model_none_type

         CASE ('pp_ne')
            model_construct_new%ne_type = model_ne_type

         CASE ('pp_te_p')
            model_construct_new%ne_type = model_ne_te_p_type

         CASE DEFAULT
            WRITE (*,*) 'ne type: ' // TRIM(ADJUSTL(ne_type)) //               &
     &                  ' setting type to none.'
            model_construct_new%ne_type = model_none_type

      END SELECT

      ALLOCATE(model_construct_new%sxrem_type(SIZE(sxrem_type)))
      DO i = 1, SIZE(sxrem_type)

         SELECT CASE (TRIM(ADJUSTL(sxrem_type(i))))

            CASE ('none')
               model_construct_new%sxrem_type(i) = model_none_type

            CASE ('pp_sxrem')
               model_construct_new%sxrem_type(i) = model_sxrem_type

            CASE DEFAULT
               WRITE (*,*)                                                     &
     &            'sxrem type: ' // TRIM(ADJUSTL(sxrem_type(i))) //            &
     &            ' setting type to none.'
               model_construct_new%sxrem_type(i) = model_none_type

         END SELECT

      END DO

      SELECT CASE (TRIM(ADJUSTL(te_type)))

         CASE ('none')
            model_construct_new%te_type = model_none_type

         CASE ('pp_te')
            model_construct_new%te_type = model_te_type

!>  @todo FIXME: Name here shouldn't have vmec in it
         CASE ('pp_ne_vmec_p', 'pp_ne_p')
            model_construct_new%te_type = model_te_ne_p_type

         CASE DEFAULT
            WRITE (*,*) 'te type: ' // TRIM(ADJUSTL(te_type)) //               &
     &                  ' setting type to none.'
            model_construct_new%te_type = model_none_type

      END SELECT

!  Cannot have the temperature derived from the denisty at the same time the
!  density is derived from the temperature.
      CALL assert(.not.(                                                       &
     &        (model_construct_new%te_type .eq. model_te_ne_p_type)            &
     &        .and.                                                            &
     &        (model_construct_new%ne_type .eq. model_ne_te_p_type)),          &
     &        'Cannot derive both the denisty and temperature' //              &
     &        ' from pressure in the same model.')

      model_construct_new%ne_unit = ne_unit
      model_construct_new%ne_min = ne_min
      model_construct_new%pressure_fraction = pressure_fraction
      model_construct_new%equilibrium => equilibrium

      grid_size = equilibrium_get_grid_size(equilibrium)

      IF (grid_size .gt. 0) THEN
         model_construct_new%grid_start =                                      &
     &      equilibrium_get_grid_start(equilibrium)
         grid_stop = equilibrium_get_grid_stop(equilibrium)

         model_construct_new%grid_step =                                       &
     &      (grid_stop - model_construct_new%grid_start)/grid_size

         ALLOCATE(model_construct_new%ne_grid(grid_size))
         ALLOCATE(model_construct_new%sxrem_grid(SIZE(sxrem_type),             &
     &                                           grid_size))
         ALLOCATE(model_construct_new%te_grid(grid_size))
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref model_class containing an @ref equilibrium by
!>  syncing through mpi.
!>
!>  Allocates memory and initializes a @ref model_class object.
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref model_class object.
!-------------------------------------------------------------------------------
      FUNCTION model_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), POINTER :: model_construct_mpi
      INTEGER, INTENT(in)         :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                     :: error
      INTEGER                     :: sxrem_size
      INTEGER                     :: grid_size

!  Start of executable code
      ALLOCATE(model_construct_mpi)

      CALL MPI_BCAST(sxrem_size, 1, MPI_INTEGER, 0, intra_comm, error)
      ALLOCATE(model_construct_mpi%sxrem_type(sxrem_size))

      CALL MPI_BCAST(model_construct_mpi%ne_type, 1, MPI_INTEGER, 0,           &
     &               intra_comm, error)
      CALL MPI_BCAST(model_construct_mpi%sxrem_type, sxrem_size,               &
     &               MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(model_construct_mpi%te_type, 1, MPI_INTEGER, 0,           &
     &               intra_comm, error)

      CALL MPI_BCAST(model_construct_mpi%ne_unit, 1, MPI_REAL8, 0,             &
     &               intra_comm, error)
      CALL MPI_BCAST(model_construct_mpi%ne_min, 1, MPI_REAL8, 0,              &
     &               intra_comm, error)
      CALL MPI_BCAST(model_construct_mpi%pressure_fraction, 1,                 &
     &               MPI_REAL8, 0, intra_comm, error)

      model_construct_mpi%equilibrium =>                                       &
     &   equilibrium_construct(intra_comm)

      grid_size =                                                              &
     &   equilibrium_get_grid_size(model_construct_mpi%equilibrium)

      IF (grid_size .gt. 0) THEN
         model_construct_mpi%grid_start =                                      &
     &      equilibrium_get_grid_start(model_construct_mpi%equilibrium)
         model_construct_mpi%grid_step =                                       &
     &      equilibrium_get_grid_start(model_construct_mpi%equilibrium)

         ALLOCATE(model_construct_mpi%ne_grid(grid_size))
         ALLOCATE(model_construct_mpi%te_grid(grid_size))
         ALLOCATE(model_construct_mpi%sxrem_grid(sxrem_size,grid_size))

         CALL MPI_BCAST(model_construct_mpi%ne_grid, grid_size,                &
     &                  MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(model_construct_mpi%sxrem_grid,                        &
     &                  grid_size*sxrem_size, MPI_REAL8, 0, intra_comm,        &
     &                  error)
         CALL MPI_BCAST(model_construct_mpi%te_grid, grid_size,                &
     &                  MPI_REAL8, 0, intra_comm, error)
      END IF
#endif

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
      TYPE (model_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%sxrem_type)) THEN
         DEALLOCATE(this%sxrem_type)
         this%sxrem_type => null()
      END IF

      IF (ASSOCIATED(this%equilibrium)) THEN
         CALL equilibrium_destruct(this%equilibrium)
         this%equilibrium => null()
      END IF

      IF (ASSOCIATED(this%ne_grid)) THEN
         DEALLOCATE(this%ne_grid)
         this%ne_grid => null()
      END IF

      IF (ASSOCIATED(this%sxrem_grid)) THEN
         DEALLOCATE(this%sxrem_grid)
         this%sxrem_grid => null()
      END IF

      IF (ASSOCIATED(this%te_grid)) THEN
         DEALLOCATE(this%te_grid)
         this%te_grid => null()
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
!>  it is deferred to the @ref equilibrium.
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
!>  @brief Get the id for a parameter.
!>
!>  Gets a model reconstruction parameter. If the id is not a model parameter,
!>  it is deferred to the @ref equilibrium.
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

         CASE ('ne_grid')
            model_get_param_id = model_ne_grid_id

         CASE ('te_grid')
            model_get_param_id = model_te_grid_id

         CASE ('sxrem_grid')
            model_get_param_id = model_sxrem_grid_id

         CASE DEFAULT
            model_get_param_id =                                               &
     &         equilibrium_get_param_id(this%equilibrium, param_name)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the value of a model parameter.
!>
!>  Gets a model reconstruction parameter value. If the id is not a model
!>  parameter, it is deferred to the @ref equilibrium.
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

!  Not every equilibrium will provide grid for these profiles.
         CASE (model_ne_grid_id)
            IF (ASSOCIATED(this%ne_grid)) THEN
               model_get_param_value = this%ne_grid(i_index)
            ELSE
               model_get_param_value = 0.0
            END IF

         CASE (model_te_grid_id)
            IF (ASSOCIATED(this%te_grid)) THEN
               model_get_param_value = this%te_grid(i_index)
            ELSE
               model_get_param_value = 0.0
            END IF

         CASE (model_sxrem_grid_id)
            IF (ASSOCIATED(this%sxrem_grid)) THEN
               model_get_param_value = this%sxrem_grid(i_index,                &
     &                                                 j_index)
            ELSE
               model_get_param_value = 0.0
            END IF

         CASE DEFAULT
            model_get_param_value =                                            &
     &         equilibrium_get_param_value(this%equilibrium, id,               &
     &                                     i_index, j_index)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the name of a model parameter.
!>
!>  Gets a model reconstruction parameter name. If the id is not a model
!>  parameter, it is deferred to the @ref equilibrium.
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

         CASE (model_ne_grid_id)
            model_get_param_name = 'ne_grid'

         CASE (model_te_grid_id)
            model_get_param_name = 'te_grid'

         CASE (model_sxrem_grid_id)
            model_get_param_name = 'sxrem_grid'

         CASE DEFAULT
            model_get_param_name =                                             &
     &         equilibrium_get_param_name(this%equilibrium, id)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron density at a cartesian position.
!>
!>  Gets the electron density at a cartesian position. Electron density is
!>  computed based on the type of @ref model_class::ne_type.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @param[in] x_cart Cartesian position to get the electron density at.
!>  @returns The electron density at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION model_get_ne_cart(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: model_get_ne_cart
      TYPE (model_class), INTENT(in)         :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%ne_type)

         CASE (model_ne_type)
            model_get_ne_cart = equilibrium_get_ne(this%equilibrium,           &
     &                                             x_cart)
            model_get_ne_cart = MAX(this%ne_min,                               &
     &                              this%ne_unit*model_get_ne_cart)

         CASE (model_ne_te_p_type)
            model_get_ne_cart = eV_per_Joule*this%pressure_fraction            &
     &                        * equilibrium_get_p(this%equilibrium,            &
     &                                            x_cart)                      &
     &                        / model_get_te(this, x_cart)

         CASE DEFAULT
             model_get_ne_cart = -2.0E49_rprec

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron density at a radial position.
!>
!>  Gets the electron density at a radial position. Electron density is computed
!>  based on the type of @ref model_class::ne_type.
!>
!>  @param[in] this A @ref model_class instance.
!>  @param[in] r    Radial position to get the electron density at.
!>  @returns The electron density at r.
!-------------------------------------------------------------------------------
      FUNCTION model_get_ne_radial(this, r)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: model_get_ne_radial
      TYPE (model_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)       :: r

!  Start of executable code
      SELECT CASE (this%ne_type)

         CASE (model_ne_type)
            model_get_ne_radial = equilibrium_get_ne(this%equilibrium,         &
     &                                               r)
            model_get_ne_radial = MAX(this%ne_min,                             &
     &                                this%ne_unit*model_get_ne_radial)

         CASE (model_ne_te_p_type)
            model_get_ne_radial = eV_per_Joule*this%pressure_fraction          &
     &                          * equilibrium_get_p(this%equilibrium,          &
     &                                              r)                         &
     &                          / model_get_te(this, r)

         CASE DEFAULT
            model_get_ne_radial = -2.0E49_rprec

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron temperature at a cartesian position.
!>
!>  Gets the electron temperature at a cartesian position. Electron temperature
!>  is computed based on the type of @ref model_class::te_type.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @param[in] x_cart Cartesian position to get the electron temperature at.
!>  @returns The electron temperature at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION model_get_te_cart(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: model_get_te_cart
      TYPE (model_class), INTENT(in)         :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%te_type)

         CASE (model_te_type)
            model_get_te_cart = equilibrium_get_te(this%equilibrium,           &
     &                                             x_cart)

         CASE (model_te_ne_p_type)
            model_get_te_cart = eV_per_Joule*this%pressure_fraction            &
     &                        * equilibrium_get_p(this%equilibrium,            &
     &                                            x_cart)                      &
     &                        / model_get_ne(this, x_cart)

         CASE DEFAULT
            model_get_te_cart = -2.0E49_rprec

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron temperature at a radial position.
!>
!>  Gets the electron temperature at a radial position. Electron temperature is
!>  computed based on the type of @ref model_class::te_type.
!>
!>  @param[in] this A @ref model_class instance.
!>  @param[in] r    Radial position to get the electron temperature at.
!>  @returns The electron temperature at r.
!-------------------------------------------------------------------------------
      FUNCTION model_get_te_radial(this, r)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: model_get_te_radial
      TYPE (model_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)       :: r

!  Start of executable code
      SELECT CASE (this%te_type)

         CASE (model_te_type)
            model_get_te_radial = equilibrium_get_te(this%equilibrium,         &
     &                                               r)

         CASE (model_te_ne_p_type)
            model_get_te_radial = eV_per_Joule*this%pressure_fraction          &
     &                          * equilibrium_get_p(this%equilibrium,          &
     &                                              r)                         &
     &                          / model_get_ne(this, r)

         CASE DEFAULT
            model_get_te_radial = -2.0E49_rprec

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the soft x-ray emissivity at a cartesian position.
!>
!>  Gets the soft x-ray emissivity at a cartesian position. Soft x-ray
!>  emissivity is computed based on the type of @ref model_class::sxrem_type.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @param[in] x_cart Cartesian position to get the soft x-ray emissivity at.
!>  @param[in] index  Index of the soft x-ray emissivity model.
!>  @returns The soft x-ray emissivity at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION model_get_sxrem_cart(this, x_cart, index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: model_get_sxrem_cart
      TYPE (model_class), INTENT(in)         :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      INTEGER, INTENT(in)                    :: index

!  Start of executable code
      SELECT CASE (this%sxrem_type(index))

         CASE (model_sxrem_type)
            model_get_sxrem_cart =                                             &
     &         equilibrium_get_sxrem(this%equilibrium, x_cart, index)

         CASE DEFAULT
            model_get_sxrem_cart = -2.0E49_rprec

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the soft x-ray emissivity at a radial position.
!>
!>  Gets the soft x-ray emissivity at a radial position. Soft x-ray emissivity
!>  is computed based on the type of @ref model_class::sxrem_type.
!>
!>  @param[in] this  A @ref model_class instance.
!>  @param[in] r     Radial position to get the soft x-ray emissivity at.
!>  @param[in] index Index of the soft x-ray emissivity model.
!>  @returns The soft x-ray emissivity at r.
!-------------------------------------------------------------------------------
      FUNCTION model_get_sxrem_radial(this, r, index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: model_get_sxrem_radial
      TYPE (model_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)       :: r
      INTEGER, INTENT(in)            :: index

!  Start of executable code
      SELECT CASE (this%sxrem_type(index))

         CASE (model_sxrem_type)
            model_get_sxrem_radial =                                           &
     &         equilibrium_get_sxrem(this%equilibrium, r, index)

         CASE DEFAULT
            model_get_sxrem_radial = -2.0E49_rprec

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
            model_get_te_type = 'pp_ne_p'

         CASE DEFAULT
            model_get_te_type = 'none'

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the soft x-ray emissivity type as a string.
!>
!>  @param[in] this   A @ref model_class instance.
!>  @param[in] index  Index of the soft x-ray emissivity model.
!>  @returns A string description of the soft x-ray emissivity type.
!-------------------------------------------------------------------------------
      FUNCTION model_get_sxrem_type(this, index)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: model_get_sxrem_type
      TYPE (model_class), INTENT(in)   :: this
      INTEGER, INTENT(in)              :: index

!  Start of executable code
      SELECT CASE (this%sxrem_type(index))

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
!>  is not a model parameter, it is deferred to the @ref equilibrium.
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
!>  @brief Solves the model.
!>
!>  This solves the equilibirum model and computes modeled values. The modeled
!>  values computed are the grided model values.
!>
!>  @param[inout] this     A @ref model_class instance.
!>  @param[inout] num_iter Counter to track the number of iterations.
!>  @param[in]    iou      Input/output unit of the file to write logs to.
!>  @returns True if the convergece was sucessful and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION model_converge(this, num_iter, iou)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                           :: model_converge
      TYPE (model_class), INTENT(inout) :: this
      INTEGER, INTENT(inout)            :: num_iter
      INTEGER, INTENT(in)               :: iou

!  local variables
      REAL (rprec)                      :: r, dr
      INTEGER                           :: i, j

!  Start of executable code
      model_converge = equilibrium_converge(this%equilibrium, num_iter,        &
     &                                      iou)

      IF (model_converge) THEN

         IF (ASSOCIATED(this%ne_grid)) THEN
            DO i = 1, SIZE(this%ne_grid)
               r = (i - 1)*this%grid_step + this%grid_start
               this%ne_grid(i) = model_get_ne(this, r)
            END DO
         END IF

         IF (ASSOCIATED(this%sxrem_grid)) THEN
            DO i = 1, SIZE(this%sxrem_grid, 2)
               r = (i - 1)*this%grid_step + this%grid_start
               DO j = 1, SIZE(this%sxrem_grid, 1)
                   this%sxrem_grid(j,i) = model_get_sxrem(this, r, j)
               END DO
            END DO
         END IF

         IF (ASSOCIATED(this%te_grid)) THEN
            DO i = 1, SIZE(this%te_grid)
               r = (i - 1)*this%grid_step + this%grid_start
               this%te_grid(i) = model_get_te(this, r)
            END DO
         END IF

      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write out the model to an output file.
!>
!>  Write out information about the model to an output file. This also then
!>  calls the @ref equilibrium_write subroutine to write out the equilibrium
!>  information as well.
!>
!>  @param[in] this A @ref model_class instance.
!>  @param[in] iou  Input/output unit of the output file.
!-------------------------------------------------------------------------------
      SUBROUTINE model_write(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: iou

!  local variables
      INTEGER                        :: i

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,*) ' *** Model Parameters'
      WRITE (iou,*) 'model_ne_type is ',                                       &
     &              TRIM(model_get_ne_type(this))

      DO i = 1, SIZE(this%sxrem_type)
         WRITE (iou,*) 'model_sxrem_type is ',                                 &
     &                 TRIM(model_get_sxrem_type(this, i))
      END DO

      WRITE (iou,*) 'model_te_type is ',                                       &
     &              TRIM(model_get_te_type(this))
      WRITE (iou, 1000) 'ne_pp_unit is ', this%ne_unit
      WRITE (iou, 1000) 'ne_min is ', this%ne_min
      WRITE (iou, 1000) 'e_pressure_fraction is ',                             &
     &                  this%pressure_fraction
1000  FORMAT(1x,a,2x,es12.5)

      CALL equilibrium_write(this%equilibrium, iou)

      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!>  @page result_file_model Model Result File
!>
!>  @tableofcontents
!>  @section result_file_model_intro_sec Introduction
!>  This page documents the contents of a result NetCDF file contributed by the
!>  model. The remaining parts of the result file are documented in
!>  the @ref result_file_main page.
!>
!>  @section result_file_model_dim_sec Dimensions
!>  @header{Dimension, Description, Code Reference}
!>  @begin_table
!>     @item{model_grid_size, Size of the radial profile grids,          equilibrium::equilibrium_get_grid_size}
!>     @item{model_num_sxrem, Number of soft x-ray emissivity profiles., v3fit_input::v3fit_max_sxrem_profiles}
!>  @end_table
!>
!>  @section result_file_model_var_sec Variables
!>  @header{Variable(Dimensions), Description, Code Reference}
!>  @begin_table
!>      @item{model_ne_type (string_len),                      Discription of the electron density model.,      model::model_class::ne_type}
!>      @item{model_te_type (string_len),                      Discription of the electron temperature model.,  model::model_class::te_type}
!>      @item{model_sxrem_type (string_len\, model_num_sxrem), Discription of the soft x-ray emissivity model., model::model_class::sxrem_type}
!>      @item{ne_unit (maxnsteps),                             Electron density unit scaling,                   model::model_class::ne_unit}
!>      @item{ne_min (maxnsteps),                              Minimum electron density,                        model::model_class::ne_min}
!>      @item{pressure_fraction (maxnsteps),                   Fractional Pressure,                             model::model_class::pressure_fraction}
!>  @end_table
!>  @table_section{result_file_model_grid_sec, Model Profiles}
!>      @item{ne_grid (model_grid_size\, maxnsteps),                      Radial profile of the electron density,          model::model_class::ne_grid}
!>      @item{te_grid (model_grid_size\, maxnsteps),                      Radial profile of the electron temperature,      model::model_class::te_grid}
!>      @item{sxrem_grid (model_grid_size\, model_num_sxrem\, maxnsteps), Radial profile of the soft x-ray emissivity,     model::model_class::sxrem_grid}
!>  @end_table
!-------------------------------------------------------------------------------
!>  @brief Define NetCDF variables for the result file
!>
!>  Defines dimensions and variables for the model contribution of the result
!>  file. Multi dimensional arrays need to be transposed so arrays appear in the
!>  correct order in non fortran languages.
!>
!>  @param[in] this             A @ref model_class instance.
!>  @param[in] result_ncid      NetCDF file id of the result file.
!>  @param[in] maxnsetps_dim_id NetCDF dimension id of the number of steps
!>                              dimension.
!-------------------------------------------------------------------------------
      SUBROUTINE model_def_result(this, result_ncid, maxnsteps_dim_id,         &
     &                            string_len_dim_id)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: result_ncid
      INTEGER, INTENT(in)            :: maxnsteps_dim_id
      INTEGER, INTENT(in)            :: string_len_dim_id

!  local variables
      INTEGER                        :: status
      INTEGER                        :: model_grid_size_dim_id
      INTEGER                        :: model_num_sxrem_dim_id
      INTEGER                        :: ne_type_var_id
      INTEGER                        :: te_type_var_id
      INTEGER                        :: sxrem_type_var_id
      INTEGER                        :: ne_unit_var_id
      INTEGER                        :: ne_min_var_id
      INTEGER                        :: pressure_fraction_var_id
      INTEGER                        :: ne_grid_var_id
      INTEGER                        :: te_grid_var_id
      INTEGER                        :: sxrem_grid_var_id

!  Start of executable code
!  Define dimensions
      IF (ASSOCIATED(this%ne_grid) .or.                                        &
     &    ASSOCIATED(this%te_grid) .or.                                        &
     &    ASSOCIATED(this%sxrem_grid)) THEN
         status = nf_def_dim(result_ncid, 'model_grid_size',                   &
     &                       SIZE(this%ne_grid),                               &
     &                       model_grid_size_dim_id)
      ELSE
         status = nf_def_dim(result_ncid, 'model_grid_size', 1,                &
     &                       model_grid_size_dim_id)
      END IF
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%sxrem_type)) THEN
         status = nf_def_dim(result_ncid, 'model_num_sxrem',                   &
     &                       SIZE(this%sxrem_type, 1),                         &
     &                       model_num_sxrem_dim_id)
      END IF
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

!  Define variables
      status = nf_def_var(result_ncid, 'model_ne_type', nf_char, 1,            &
     &                    (/ string_len_dim_id /), ne_type_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'model_te_type', nf_char, 1,            &
     &                    (/ string_len_dim_id /), te_type_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%sxrem_type)) THEN
         status = nf_def_var(result_ncid, 'model_sxrem_type', nf_char,         &
     &                       2, (/ string_len_dim_id,                          &
     &                             model_num_sxrem_dim_id /),                  &
     &                       sxrem_type_var_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      status = nf_def_var(result_ncid, 'ne_unit', nf_double, 1,                &
     &                    (/ maxnsteps_dim_id /), ne_unit_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'ne_min', nf_double, 1,                 &
     &                    (/ maxnsteps_dim_id /), ne_min_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'pressure_fraction', nf_double,         &
     &                    1, (/ maxnsteps_dim_id /),                           &
     &                    pressure_fraction_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%ne_grid)) THEN
         status = nf_def_var(result_ncid, 'ne_grid', nf_double, 2,             &
     &                       (/ model_grid_size_dim_id,                        &
     &                          maxnsteps_dim_id /), ne_grid_var_id)
      END IF

      IF (ASSOCIATED(this%te_grid)) THEN
         status = nf_def_var(result_ncid, 'te_grid', nf_double, 2,             &
     &                       (/ model_grid_size_dim_id,                        &
     &                          maxnsteps_dim_id /), te_grid_var_id)
      END IF

      IF (ASSOCIATED(this%sxrem_grid)) THEN
         status = nf_def_var(result_ncid, 'sxrem_grid', nf_double, 3,          &
     &                       (/ model_grid_size_dim_id,                        &
     &                          model_num_sxrem_dim_id,                        &
     &                          maxnsteps_dim_id /), sxrem_grid_var_id)
      END IF

      CALL equilibrium_def_result(this%equilibrium, result_ncid,               &
     &                            maxnsteps_dim_id)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write inital data to NetCDF result file
!>
!>  This write the data that doesn't change.
!>
!>  @param[in] this        A @ref model_class instance.
!>  @param[in] result_ncid NetCDF file id of the result file.
!-------------------------------------------------------------------------------
      SUBROUTINE model_write_init_data(this, result_ncid)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: result_ncid

!  Local variables
      INTEGER                        :: i, status
      INTEGER                        :: ne_type_var_id
      INTEGER                        :: te_type_var_id
      INTEGER                        :: sxrem_type_var_id

!  Start of executable code
      status = nf_inq_varid(result_ncid, 'model_ne_type',                      &
     &                      ne_type_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var_text(result_ncid, ne_type_var_id,                    &
     &                         model_get_ne_type(this))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'model_te_type',                      &
     &                      te_type_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var_text(result_ncid, te_type_var_id,                    &
     &                         model_get_te_type(this))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%sxrem_type)) THEN
         status = nf_inq_varid(result_ncid, 'model_sxrem_type',                &
     &                         sxrem_type_var_id)
         DO i = 1, SIZE(this%sxrem_type)
            status = nf_put_vara_text(result_ncid, sxrem_type_var_id,          &
     &                                (/ 1, i /),                              &
     &                                (/ data_name_length, 1 /),               &
     &                                model_get_sxrem_type(this, i))
            CALL assert_eq(status, nf_noerr, nf_strerror(status))
         END DO
      END IF

      CALL equilibrium_write_init_data(this%equilibrium, result_ncid)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write step data to NetCDF result file
!>
!>  This write the data that can change each reconstruction step.
!>
!>  @param[in] this         A @ref model_class instance.
!>  @param[in] result_ncid  NetCDF file id of the result file.
!>  @param[in] current_step Step index to write variables to.
!-------------------------------------------------------------------------------
      SUBROUTINE model_write_step_data(this, result_ncid, current_step)
      USE ezcdf

!  Declare Arguments
      TYPE (model_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: result_ncid
      INTEGER, INTENT(in)            :: current_step

!  Local variables
      INTEGER                        :: i, status
      INTEGER                        :: ne_unit_var_id
      INTEGER                        :: ne_min_var_id
      INTEGER                        :: pressure_fraction_var_id
      INTEGER                        :: ne_grid_var_id
      INTEGER                        :: te_grid_var_id
      INTEGER                        :: sxrem_grid_var_id

!  Start of executable code
      status = nf_inq_varid(result_ncid, 'ne_unit', ne_unit_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var1_double(result_ncid, ne_unit_var_id,                 &
     &                            current_step, this%ne_unit)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'ne_min', ne_min_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var1_double(result_ncid, ne_min_var_id,                  &
     &                            current_step, this%ne_min)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'pressure_fraction',                  &
     &                      pressure_fraction_var_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var1_double(result_ncid, pressure_fraction_var_id,       &
     &                            current_step, this%pressure_fraction)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%ne_grid)) THEN
         status = nf_inq_varid(result_ncid, 'ne_grid', ne_grid_var_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
         status = nf_put_vara_double(result_ncid, ne_grid_var_id,              &
     &                               (/ 1, current_step /),                    &
     &                               (/ SIZE(this%ne_grid), 1 /),              &
     &                               this%ne_grid)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%te_grid)) THEN
         status = nf_inq_varid(result_ncid, 'te_grid', te_grid_var_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
         status = nf_put_vara_double(result_ncid, te_grid_var_id,              &
     &                               (/ 1, current_step /),                    &
     &                               (/ SIZE(this%te_grid), 1 /),              &
     &                               this%te_grid)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%sxrem_grid)) THEN
         status = nf_inq_varid(result_ncid, 'sxrem_grid',                      &
     &                         sxrem_grid_var_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
         status = nf_put_vara_double(result_ncid, sxrem_grid_var_id,           &
     &                               (/ 1, 1, current_step /),                 &
     &                               (/ SIZE(this%sxrem_grid, 2),              &
     &                                  SIZE(this%sxrem_grid, 1), 1 /),        &
     &                               this%sxrem_grid)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      CALL equilibrium_write_step_data(this%equilibrium, result_ncid,          &
     &                                 current_step)

      END SUBROUTINE

!*******************************************************************************
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize data to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref model_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE model_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), INTENT(inout) :: this
      INTEGER, INTENT(in)               :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                           :: error
      INTEGER                           :: sxrem_size
      INTEGER                           :: grid_size

!  Start of executable code
      sxrem_size = SIZE(this%sxrem_type)

      CALL MPI_BCAST(SIZE(this%sxrem_type), 1, MPI_INTEGER, 0,                 &
     &               intra_comm, error)

      CALL MPI_BCAST(this%ne_type, 1, MPI_INTEGER, 0, intra_comm,              &
     &               error)
      CALL MPI_BCAST(this%sxrem_type, sxrem_size, MPI_INTEGER, 0,              &
     &               intra_comm, error)
      CALL MPI_BCAST(this%te_type, 1, MPI_INTEGER, 0, intra_comm,              &
     &               error)

      CALL MPI_BCAST(this%ne_unit, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%ne_min, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%pressure_fraction, 1, MPI_REAL8, 0,                  &
     &               intra_comm, error)

      CALL equilibrium_sync_children(this%equilibrium, intra_comm)

      grid_size = equilibrium_get_grid_size(this%equilibrium)
      IF (grid_size .gt. 0) THEN
         CALL MPI_BCAST(this%ne_grid, grid_size, MPI_REAL8, 0,                 &
     &                  intra_comm, error)
         CALL MPI_BCAST(this%sxrem_grid, grid_size*sxrem_size,                 &
     &                  MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(this%te_grid, grid_size, MPI_REAL8, 0,                 &
     &                  intra_comm, error)
      END IF
#endif

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Syncronize the model state to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref model_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE model_sync_state(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model_class), INTENT(inout) :: this
      INTEGER, INTENT(in)               :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                           :: error
      INTEGER                           :: sxrem_size
      INTEGER                           :: grid_size

!  Start of executable code
      CALL equilibrium_sync_state(this%equilibrium, intra_comm)

      grid_size = equilibrium_get_grid_size(this%equilibrium)
      IF (grid_size .gt. 0) THEN
         CALL MPI_BCAST(this%ne_grid, grid_size, MPI_REAL8, 0,                 &
     &                  intra_comm, error)
         CALL MPI_BCAST(this%sxrem_grid, grid_size*sxrem_size,                 &
     &                  MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(this%te_grid, grid_size, MPI_REAL8, 0,                 &
     &                  intra_comm, error)
      END IF

#endif

      END SUBROUTINE

      END MODULE
