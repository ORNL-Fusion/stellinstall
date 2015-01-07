!*******************************************************************************
!>  @file sxrem_ratio.f
!>  @brief Contains module sxrem_ratio
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines a feedback signal based on the temperature based on the ration of
!>  the soft x-ray emissivity profile. This signal feeds back the temperature
!>  profile at a specified positon.
!>
!>  @par Super Class:
!>  @ref feedback
!*******************************************************************************

      MODULE sxrem_ratio

      USE stel_kinds, only: rprec, dp
      USE model

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) sxrem_ratio base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a sxrem_ratio signal.
!>  @par Super Class:
!>  @ref feedback
!-------------------------------------------------------------------------------
      TYPE sxrem_ratio_class
!>  Radial position to feedback the
         REAL (rprec), DIMENSION(3) :: x_cart
!>  Indicies of soft x-ray emissivity profiles to take the ratios of.
         INTEGER, DIMENSION(2)      :: indices
      END TYPE sxrem_ratio_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref sxrem_ratio_class types either
!>  using
!>  @ref sxrem_ratio_construct_new or @ref sxrem_ratio_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE sxrem_ratio_construct
         MODULE PROCEDURE sxrem_ratio_construct_new,                           &
     &                    sxrem_ratio_construct_mpi
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref sxrem_ratio_class object.
!>
!>  Allocates memory and initializes a @ref sxrem_ratio_class object.
!>
!>  @param[in] a_model    An instance of a @ref model object.
!>  @param[in] param_name The name the parameter.
!>  @param[in] indices    Indicies of the parameter.
!>  @returns A pointer to a constructed @ref sxrem_ratio_class object.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_ratio_construct_new(x_cart, indices)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_ratio_class), POINTER :: sxrem_ratio_construct_new
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      INTEGER, DIMENSION(2), INTENT(in)      :: indices

!  local variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(sxrem_ratio_construct_new)

      sxrem_ratio_construct_new%x_cart = x_cart
      sxrem_ratio_construct_new%indices = indices

      CALL profiler_set_stop_time('sxrem_ratio_construct_new',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref sxrem_ratio_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref sxrem_ratio_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref sxrem_ratio_class object.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_ratio_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_ratio_class), POINTER :: sxrem_ratio_construct_mpi
      INTEGER, INTENT(in)               :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                           :: error

!  Start of executable code
      ALLOCATE(sxrem_ratio_construct_mpi)

      CALL MPI_BCAST(sxrem_ratio_construct_mpi%x_cart, 3,                      &
     &               MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(sxrem_ratio_construct_mpi%indices, 2,                     &
     &               MPI_INTEGER, 0, intra_comm, error)

#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref sxrem_ratio_class object.
!>
!>  Deallocates memory and uninitializes a @ref sxrem_ratio_class object.
!>
!>  @param[inout] this A @ref sxrem_ratio_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE sxrem_ratio_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_ratio_class), POINTER :: this

!  Start of executable code
      this%x_cart = 0.0
      this%indices = 0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTERS SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  The modeled signal is the ratio of the soft x-ray emissivity profiles.
!>
!>  @param[in]    this    A @ref sxrem_ratio_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_ratio_get_modeled_signal(this, a_model, sigma)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: sxrem_ratio_get_modeled_signal
      TYPE (sxrem_ratio_class), INTENT(in)    :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(:), INTENT(out) :: sigma

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxrem_ratio_get_modeled_signal =                                         &
     &   model_get_sxrem_cart(a_model, this%x_cart, 1.0_dp,                    &
     &                        this%indices(1))/                                &
     &   model_get_sxrem_cart(a_model, this%x_cart, 1.0_dp,                    &
     &                        this%indices(2))
      sigma = 0.0

      CALL profiler_set_stop_time('sxrem_ratio_get_modeled_signal',            &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the observed signal.
!>
!>  The observed signal is the ratio computed as the ratio as a function of the
!>  temperature.
!>
!>  @param[in]    this    A @ref sxrem_ratio_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @returns The observed value.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_ratio_get_observed_signal(this, a_model)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: sxrem_ratio_get_observed_signal
      TYPE (sxrem_ratio_class), INTENT(in) :: this
      TYPE (model_class), INTENT(inout)    :: a_model

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxrem_ratio_get_observed_signal =                                        &
     &   model_get_sxrem_ratio(a_model, model_get_te_cart(a_model,             &
     &                                                    this%x_cart))

      CALL profiler_set_stop_time('sxrem_ratio_get_observed_signal',           &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the feedback type.
!>
!>  Returns a description of the feedback type for use when writting output
!>  files.
!>
!>  @param[in] this A @ref prior_gaussian_class instance.
!>  @returns A string describing the prior type.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_ratio_get_signal_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: sxrem_ratio_get_signal_type
      TYPE (sxrem_ratio_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxrem_ratio_get_signal_type = 'sxrem ratio'

      CALL profiler_set_stop_time('sxrem_ratio_get_signal_type',               &
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
!>  @param[in] this    A @ref sxrem_ratio_class instance.
!>  @param[in] iou     A input/output representing the file to write to.
!-------------------------------------------------------------------------------
      SUBROUTINE sxrem_ratio_write_auxiliary(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_ratio_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: iou

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (iou,1000) this%indices(1), this%indices(2)

      CALL profiler_set_stop_time('sxrem_ratio_write_auxiliary',               &
     &                            start_time)

1000  FORMAT('sxrem profile ratio',1x,i2,'/',1x,i2)

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
!>  @param[inout] this       A @ref sxrem_ratio_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE sxrem_ratio_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_ratio_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                                 :: error
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL MPI_BCAST(this%x_cart, 3, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%indices, 2, MPI_INTEGER, 0, intra_comm, error)

      CALL profiler_set_stop_time('prior_gaussian_sync_children',              &
     &                            start_time)

#endif

      END SUBROUTINE

      END MODULE
