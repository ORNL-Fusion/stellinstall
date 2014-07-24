!*******************************************************************************
!>  @file geometric.f
!>  @brief Contains module @ref geometric
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref geometric_class.
!>  @par Super Class:
!>  @ref signal
!>  @par Sub Classes:
!>  @ref limiter
!*******************************************************************************

      MODULE geometric

      USE limiter

      IMPLICIT NONE

!*******************************************************************************
!  geometric module parameters
!*******************************************************************************
!>  Type descriptor for the @ref geometric_class no sub class.
      INTEGER, PARAMETER :: geometric_no_type      = -1
!>  Type descriptor for the @ref geometric_class subclass @ref limiter
      INTEGER, PARAMETER :: geometric_limiter_type = 0

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) geometric base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a geometric signal.
!>  @par Super class:
!>  @ref signal
!>  @par Sub classes:
!>  @ref limiter
!-------------------------------------------------------------------------------
      TYPE geometric_class
!>  Type descirptor of the geometric subclass type.
!>  @par Possible values are:
!>  * @ref geometric_no_type
!>  * @ref geometric_limiter_type,
         INTEGER                        :: type = geometric_no_type
!>  A reference to a limiter object.
         TYPE (limiter_class),  POINTER :: limiter => null()
      END TYPE geometric_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref geometric_class types either using
!>  @ref geometric_construct_limiter
!-------------------------------------------------------------------------------
      INTERFACE geometric_construct
         MODULE PROCEDURE geometric_construct_limiter,                         &
     &                    geometric_construct_mpi
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref geometric_class containing a @ref limiter object.
!>
!>  Allocates memory and initializes a @ref geometric_class object.
!>
!>  @param[in] limiter_object An instance of a @ref limiter subclass.
!>  @returns A pointer to a constructed @ref geometric_class object.
!-------------------------------------------------------------------------------
      FUNCTION geometric_construct_limiter(limiter_object)

      IMPLICIT NONE

!  Declare Arguments
	TYPE (geometric_class), POINTER :: geometric_construct_limiter
      TYPE (limiter_class), POINTER   :: limiter_object

!  Start of executable code
      ALLOCATE(geometric_construct_limiter)

      geometric_construct_limiter%type = geometric_limiter_type
      geometric_construct_limiter%limiter => limiter_object

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref geometric_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref geometric_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref geometric_class object.
!-------------------------------------------------------------------------------
      FUNCTION geometric_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (geometric_class), POINTER :: geometric_construct_mpi
      INTEGER, INTENT(in)             :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                         :: error

!  Start of executable code
      ALLOCATE(geometric_construct_mpi)

      CALL MPI_BCAST(geometric_construct_mpi%type, 1, MPI_INTEGER, 0,          &
     &               intra_comm, error)
      geometric_construct_mpi%limiter => limiter_construct(intra_comm)

#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref geometric_class object.
!>
!>  Deallocates memory and uninitializes a @ref geometric_class object.
!>
!>  @param[inout] this A @ref geometric_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE geometric_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (geometric_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%limiter)) THEN
         CALL limiter_destruct(this%limiter)
         this%limiter => null()
      END IF

      this%type = geometric_no_type

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_modeled_signal method of the subclass instance
!>  @ref geometric_class was constructed with.
!>
!>  @see limiter::limiter_get_modeled_signal
!>
!>  @param[in]    this    A @ref geometric_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION geometric_get_modeled_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)       :: geometric_get_modeled_signal
      TYPE (geometric_class), INTENT(in)      :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Start of executable code
      SELECT CASE(this%type)
         CASE (geometric_limiter_type)
            geometric_get_modeled_signal =                                     &
     &         limiter_get_modeled_signal(this%limiter, a_model,               &
     &                                    sigma)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the geometric type.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_signal_type method of the subclass instance
!>  @ref geometric_class was constructed with.
!>
!>  @see limiter::limiter_get_signal_type
!>
!>  @param[in] this A @ref geometric_class instance.
!>  @returns A string describing the geometric type.
!-------------------------------------------------------------------------------
      FUNCTION geometric_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: geometric_get_signal_type
      TYPE (geometric_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (geometric_limiter_type)
            geometric_get_signal_type = 'geometric ' //                        &
     &         TRIM(limiter_get_signal_type(this%limiter))

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_header method of the subclass instance @ref geometric_class
!>  was constructed with.
!>
!>  @see limiter::limiter_get_header
!>
!>  @param[in] this A @ref geometric_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION geometric_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   geometric_get_header
      TYPE (geometric_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (geometric_limiter_type)
            geometric_get_header = limiter_get_header(this%limiter)

         CASE DEFAULT
            geometric_get_header = 'N/A'

      END SELECT

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
!>  @param[inout] this       A @ref geometric_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE geometric_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (geometric_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                   :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                               :: error

!  Start of executable code
      CALL MPI_BCAST(this%type, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL limiter_sync_children(this%limiter, intra_comm)

#endif

      END SUBROUTINE

      END MODULE
