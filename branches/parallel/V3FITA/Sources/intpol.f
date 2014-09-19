!*******************************************************************************
!>  @file intpol.f
!>  @brief Contains module @ref intpol.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Implements interferometry/polarimetry diagnostic. Defines the base class of
!>  the type @ref intpol_class.
!>  @par Super Class:
!>  @ref diagnostic
!*******************************************************************************

      MODULE intpol

      USE stel_kinds, only: rprec, dp
      USE integration_path
      USE coordinate_utilities
      USE profiler

      IMPLICIT NONE
!*******************************************************************************
!  diagnostic module parameters
!*******************************************************************************
!>  Type descriptor for intpol no type.
      INTEGER, PARAMETER :: intpol_no_type             = -1
!>  Type descriptor for intpol type represented a interferometer diagnostic.
      INTEGER, PARAMETER :: intpol_interferometry_type = 0
!>  Type descriptor for intpol type represented a polarmetry diagnostic.
      INTEGER, PARAMETER :: intpol_polarimetry_type    = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) intpol base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a intpol signal.
!>  @par Super Class:
!>  @ref diagnostic
!-------------------------------------------------------------------------------
      TYPE intpol_class
!>  Type descirptor of the intpol type.
!>  @par Possible values are:
!>  * @ref intpol_no_type
!>  * @ref intpol_interferometry_type
!>  * @ref intpol_polarimetry_type
         INTEGER                :: type
!>  Wavelength if the beam.
         REAL (rprec)           :: wavelength
!>  Controls if the results is computed in degrees or radians.
         LOGICAL                :: in_degrees
!>  The complete path of the chord.
         TYPE (vertex), POINTER :: chord_path => null()
      END TYPE intpol_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref intpol_class types using
!>  @ref intpol_construct_int, @ref intpol_construct_pol or
!>  @ref intpol_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE intpol_construct
         MODULE PROCEDURE intpol_construct_int,                                &
     &                    intpol_construct_pol,                                &
     &                    intpol_construct_mpi
      END INTERFACE

      PRIVATE :: int_function, pol_function

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref intpol_class object representing a interferometer
!>  diagnostic.
!>
!>  Allocates memory and initializes a @ref intpol_class object.
!>  @param[in] chord_paths The nodes of a multi point segment integration path.
!>  @returns A pointer to a constructed @ref intpol_class object.
!-------------------------------------------------------------------------------
      FUNCTION intpol_construct_int(chord_paths)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), POINTER :: intpol_construct_int
      REAL (rprec), DIMENSION(:,:) :: chord_paths

!  local variables
      INTEGER                      :: i
      REAL (rprec)                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(intpol_construct_int)

      DO i = 1, SIZE(chord_paths, 1)
         CALL path_append_vertex(intpol_construct_int%chord_path,              &
     &                           chord_paths(i,:))
      END DO

      intpol_construct_int%type = intpol_interferometry_type
      intpol_construct_int%wavelength = 0.0
      intpol_construct_int%in_degrees = .false.

      CALL profiler_set_stop_time('intpol_construct_int', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref intpol_class object representing a polarimetry
!>  diagnostic.
!>
!>  Allocates memory and initializes a @ref intpol_class object.
!>
!>  @param[in] wavelength  Wavelength of the polarimetry beam.
!>  @param[in] in_degrees  Specifies if the sigals is degrees or radians. True
!>                         specifies degree while false specifies radians.
!>  @param[in] chord_paths The nodes of a multi point segment integration path.
!>  @returns A pointer to a constructed @ref intpol_class object.
!-------------------------------------------------------------------------------
      FUNCTION intpol_construct_pol(wavelength, in_degrees, chord_paths)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), POINTER :: intpol_construct_pol
      REAL (rprec), INTENT(in)     :: wavelength
      LOGICAL, INTENT(in)          :: in_degrees
      REAL (rprec), DIMENSION(:,:) :: chord_paths

!  local variables
      INTEGER                      :: i
      REAL (rprec)                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(intpol_construct_pol)

      DO i = 1, SIZE(chord_paths, 1)
         CALL path_append_vertex(intpol_construct_pol%chord_path,              &
     &                           chord_paths(i,:))
      END DO

      intpol_construct_pol%type = intpol_polarimetry_type
      intpol_construct_pol%wavelength = wavelength
      intpol_construct_pol%in_degrees = in_degrees

      CALL profiler_set_stop_time('intpol_construct_pol', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref sxrem_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref intpol_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref intpol_class object.
!-------------------------------------------------------------------------------
      FUNCTION intpol_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), POINTER :: intpol_construct_mpi

      INTEGER, INTENT(in)          :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                      :: error

!  Start of executable code
      ALLOCATE(intpol_construct_mpi)

      intpol_construct_mpi%chord_path => path_construct(intra_comm)

      CALL MPI_BCAST(intpol_construct_mpi%type, 1, MPI_INTEGER, 0,             &
     &               intra_comm, error)
      CALL MPI_BCAST(intpol_construct_mpi%wavelength, 1, MPI_REAL8,            &
     &               0, intra_comm, error)
      CALL MPI_BCAST(intpol_construct_mpi%in_degrees, 1, MPI_LOGICAL,          &
     &               0, intra_comm, error)

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
      SUBROUTINE intpol_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), POINTER :: this

!  Start of executable code
      this%type = intpol_no_type
      this%wavelength = 0.0
      this%in_degrees = .false.

      IF (ASSOCIATED(this%chord_path)) THEN
         CALL path_destruct(this%chord_path)
         this%chord_path => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  Calculates the modeled signal by integrating along the chord path. The
!>  integration function is provided by @ref int_function for interferometer
!>  diagnostics and @ ref pol_function for polarimetry diagnostics.
!>
!>  @param[in]    this    A @ref intpol_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_modeled_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: intpol_get_modeled_signal
      TYPE (intpol_class), INTENT(inout)      :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  local variables
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      REAL (rprec)                            :: start_time

!  local parameters
      REAL(rprec), PARAMETER :: polar_constant = 2.62E-13_dp

!  Start of executable code
      start_time = profiler_get_start_time()

      intpol_get_modeled_signal = 0.0

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

      sigma = 0.0
      DEALLOCATE(context)

      CALL profiler_set_stop_time('intpol_get_modeled_signal',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the intpol type.
!>
!>  Returns a description of the intpol type for use when writting output files.
!>
!>  @param[in] this A @ref intpol_class instance.
!>  @returns A string describing the intpol type.
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: intpol_get_signal_type
      TYPE (intpol_class), INTENT(in)  :: this

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      intpol_get_signal_type = 'ipch'

      CALL profiler_set_stop_time('intpol_get_signal_type', start_time)

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
!>  @param[inout] this       A @ref intpol_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE intpol_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                            :: error
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL path_sync_children(this%chord_path, intra_comm)
      CALL MPI_BCAST(this%type, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%wavelength, 1, MPI_REAL8, 0, intra_comm,             &
     &               error)
      CALL MPI_BCAST(this%in_degrees, 1, MPI_LOGICAL, 0, intra_comm,           &
     &               error)

      CALL profiler_set_stop_time('intpol_sync_children', start_time)
#endif

       END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Interferometer callback function.
!>
!>  Returns the value of the denisty times the change in path length. This
!>  function is passed to @ref integration_path::path_integrate to act as a
!>  callback. The denisty is provided by @ref model::model_get_ne.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref the model.
!>  @param[in] xcart   A the integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] dx      A scalar change in path.
!>  @returns The ne(x)*dl at point x.
!>
!>  @note context in this case if not that same as @ref context. It refers to
!>  generic passed to the function. This is the fortran eqivalent to a void
!>  pointer in C.
!-------------------------------------------------------------------------------
      FUNCTION int_function(context, xcart, dxcart, dx)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: int_function

!  local variables
      TYPE(model_class)                     :: a_model
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      a_model = TRANSFER(context, a_model)
      int_function = model_get_ne(a_model, xcart)*dx

      CALL profiler_set_stop_time('int_function', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Polarmetry callback function.
!>
!>  Returns the value of the denisty times the dot product of the magnetic
!>  field vector and the path direction. This function is passed to
!>  @ref integration_path::path_integrate to act as a callback. The denisty is
!>  provided by @ref model::model_get_ne. The magnetic field vector is proved by
!>  @ref equilibrium::equilibrium_get_B_vec
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref the model.
!>  @param[in] xcart   A the integration point.
!>  @param[in] dxcart  A vector change in path.
!>  @param[in] dx      A scalar change in path. Not used in this function.
!>  @returns The ne(x)*dl at point x.
!>
!>  @note context in this case if not that same as @ref context. It refers to
!>  generic passed to the function. This is the fortran eqivalent to a void
!>  pointer in C.
!-------------------------------------------------------------------------------
      FUNCTION pol_function(context, xcart, dxcart, dx)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: pol_function

!  local variables
      TYPE(model_class)                     :: a_model
      REAL(rprec), DIMENSION(3)             :: bcart
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      a_model = TRANSFER(context, a_model)

      bcart = equilibrium_get_B_vec(a_model%equilibrium, xcart)
      pol_function =                                                           &
     &   model_get_ne(a_model, xcart)*DOT_PRODUCT(bcart, dxcart)

      CALL profiler_set_stop_time('pol_function', start_time)

      END FUNCTION

      END MODULE
