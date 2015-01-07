!*******************************************************************************
!>  @file prior_gaussian.f
!>  @brief Contains module prior_gaussian
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref prior_gaussian_class. This class
!>  implements priors of the type
!>
!>  P(B) = Exp(-((B - mu)/sigma)^2)
!>
!>  @par Super Class:
!>  @ref prior
!*******************************************************************************

      MODULE prior_gaussian

      USE data_parameters
      USE model
      USE stel_kinds, only : rprec, dp

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) prior_gaussian base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a prior_guassian signal.
!>  @par Super Class:
!>  @ref signal
!>  @par Sub Classes:
!>  @ref prior_guassian
!-------------------------------------------------------------------------------
      TYPE prior_gaussian_class
!>  Parameter id to check the guassian against.
         INTEGER               :: param_id = data_no_id
!>  Indices of the parameter.
         INTEGER, DIMENSION(2) :: indices = 0
      END TYPE prior_gaussian_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref prior_gaussian_class types either
!>  using
!>  @ref prior_gaussian_construct_new or @ref prior_gaussian_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE prior_gaussian_construct
         MODULE PROCEDURE prior_gaussian_construct_new,                        &
     &                    prior_gaussian_construct_mpi
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref prior_gaussian_class object.
!>
!>  Allocates memory and initializes a @ref prior_gaussian_class object.
!>
!>  @param[in] a_model    An instance of a @ref model object.
!>  @param[in] param_name The name the parameter.
!>  @param[in] indices    Indicies of the parameter.
!>  @returns A pointer to a constructed @ref prior_gaussian_class object.
!-------------------------------------------------------------------------------
      FUNCTION prior_gaussian_construct_new(a_model, param_name,               &
     &                                      indices)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), POINTER ::                                  &
     &   prior_gaussian_construct_new
      TYPE (model_class), INTENT(in)       :: a_model
      CHARACTER (len=*), INTENT(in)        :: param_name
      INTEGER, DIMENSION(2), INTENT(in)    :: indices

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(prior_gaussian_construct_new)

      prior_gaussian_construct_new%param_id =                                  &
     &   model_get_param_id(a_model, TRIM(param_name))
      prior_gaussian_construct_new%indices = indices

      CALL profiler_set_stop_time('prior_gaussian_construct_new',              &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref prior_gaussian_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref prior_gaussian_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref prior_gaussian_class object.
!-------------------------------------------------------------------------------
      FUNCTION prior_gaussian_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), POINTER ::                                  &
     &   prior_gaussian_construct_mpi
      INTEGER, INTENT(in)                  :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                              :: error

!  Start of executable code
      ALLOCATE(prior_gaussian_construct_mpi)

      CALL MPI_BCAST(prior_gaussian_construct_mpi%param_id, 1,                 &
     &               MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(prior_gaussian_construct_mpi%indices, 2,                  &
     &               MPI_INTEGER, 0, intra_comm, error)

#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref prior_gaussian_class object.
!>
!>  Deallocates memory and uninitializes a @ref prior_gaussian_class object.
!>
!>  @param[inout] this A @ref prior_gaussian_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE prior_gaussian_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), POINTER :: this

!  Start of executable code
      this%param_id = data_no_id
      this%indices = 0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTERS SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  The modeled signal is just the value of the parameter.
!>
!>  @param[in]    this    A @ref prior_gaussian_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION prior_gaussian_get_modeled_signal(this, a_model, sigma)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: prior_gaussian_get_modeled_signal
      TYPE (prior_gaussian_class), INTENT(in) :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(:), INTENT(out) :: sigma

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      prior_gaussian_get_modeled_signal =                                      &
     &   model_get_param_value(a_model, this%param_id,                         &
     &                         this%indices(1), this%indices(2))
      sigma = 0.0

      CALL profiler_set_stop_time('prior_gaussian_get_modeled_signal',         &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the prior type.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_signal_type method of the subclass instance
!>  @ref prior_class was constructed with.
!>
!>  @see prior_gaussian::prior_gaussia_get_signal_type
!>
!>  @param[in] this A @ref prior_gaussian_class instance.
!>  @returns A string describing the prior type.
!-------------------------------------------------------------------------------
      FUNCTION prior_gaussian_get_signal_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: prior_gaussian_get_signal_type
      TYPE (prior_gaussian_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      prior_gaussian_get_signal_type = 'guassian'

      CALL profiler_set_stop_time('prior_gaussian_get_signal_type',            &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out auxiliary signal information to an output file.
!>
!>  Writes out the s_name and coefficient of the combined signals.
!>
!>  @param[in] this    A @ref prior_gaussian_class instance.
!>  @param[in] iou     A input/output representing the file to write to.
!>  @param[in] a_model The equilibrium model.
!-------------------------------------------------------------------------------
      SUBROUTINE prior_gaussian_write_auxiliary(this, iou, a_model)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), INTENT(in) :: this
      INTEGER, INTENT(in)                     :: iou
      TYPE (model_class), INTENT(in)          :: a_model

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (iou,1000)
      WRITE (iou,1001) model_get_param_name(a_model, this%param_id),           &
     &                 this%indices

      CALL profiler_set_stop_time('prior_gaussian_write_auxiliary',            &
     &                            start_time)

1000  FORMAT('parameter name',8x,'inx1',2x,'inx2')
1001  FORMAT(a20,2(2x,i4))

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
!>  @param[inout] this       A @ref prior_gaussian_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE prior_gaussian_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                        :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                                    :: error
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL MPI_BCAST(this%param_id, 1, MPI_INTEGER, 0, intra_comm,             &
     &               error)
      CALL MPI_BCAST(this%indices, 2, MPI_INTEGER, 0, intra_comm, error)

      CALL profiler_set_stop_time('prior_gaussian_sync_children',              &
     &                            start_time)

#endif

      END SUBROUTINE

      END MODULE
