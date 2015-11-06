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
      USE mpi_inc
      USE profiler

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
!>  Position of the thomson scattering point.
         REAL (rprec), DIMENSION(3) :: xcart
      END TYPE thomson_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref thomson_class types using
!>  @ref thomson_construct_new or @ref thomson_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE thomson_construct
         MODULE PROCEDURE thomson_construct_new,                               &
     &                    thomson_construct_mpi
      END INTERFACE

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

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      thomson_construct_te =>                                                  &
     &   thomson_construct(thomson_temperature_type, xcart)

      CALL profiler_set_stop_time('thomson_construct_te', start_time)

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

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      thomson_construct_ne => thomson_construct(thomson_density_type,          &
     &                                          xcart)

      CALL profiler_set_stop_time('thomson_construct_ne', start_time)

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

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      thomson_construct_p => thomson_construct(thomson_pressure_type,          &
     &                                         xcart)

      CALL profiler_set_stop_time('thomson_construct_p', start_time)

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
      FUNCTION thomson_construct_new(type, xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct_new
      INTEGER, INTENT(in)                    :: type
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(thomson_construct_new)

      thomson_construct_new%type = type
      thomson_construct_new%xcart = xcart

      CALL profiler_set_stop_time('thomson_construct_new', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref thomson_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref thomson_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref thomson_class object.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER :: thomson_construct_mpi
      INTEGER, INTENT(in)           :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                       :: error

!  Start of executable code
      ALLOCATE(thomson_construct_mpi)

      CALL MPI_BCAST(thomson_construct_mpi%type, 1, MPI_INTEGER, 0,            &
     &               intra_comm, error)
      CALL MPI_BCAST(thomson_construct_mpi%xcart, 3, MPI_REAL8, 0,             &
     &               intra_comm, error)

#endif

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

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

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

      CALL profiler_set_stop_time('thomson_get_modeled_signal',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the thomson type.
!>
!>  Returns a description of the thomson type for use when writting output
!>  files.
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

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      thomson_get_signal_type = 'thscte'

      CALL profiler_set_stop_time('thomson_get_signal_type', start_time)

      END FUNCTION

!*******************************************************************************
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize data to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref thomson_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE thomson_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                 :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                            :: error
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL MPI_BCAST(this%type, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%xcart, 3, MPI_REAL8, 0, intra_comm, error)

      CALL profiler_set_stop_time('thomson_sync_children', start_time)

#endif

      END SUBROUTINE

      END MODULE
