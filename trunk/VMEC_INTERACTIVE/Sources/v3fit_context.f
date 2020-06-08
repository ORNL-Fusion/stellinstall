      MODULE v3fit_context
      
      use v3_utilities
      use stel_kinds
      use profiler
      
      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) v3fit context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a v3fit context. This contains all memory needed to
!>  operate v3fit.
!-------------------------------------------------------------------------------
      TYPE v3fit_context_class
!>  File input output unit for the runlog file. This file is a human readable
!>  output history of a v3fit run.
         INTEGER                                      :: runlog_iou

!--- MPI -----------------------------------------------------------------------
#if defined (MPI_OPT)
!>  MPI communicator reference for full domain.
         INTEGER :: global_comm = MPI_COMM_NULL
!>  MPI communicator reference for the equilibrium domain.
         INTEGER :: equilibrium_comm = MPI_COMM_NULL
!>  MPI communicator reference for the reconstruction domain.
         INTEGER :: reconstruction_comm = MPI_COMM_NULL
#endif
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref v3fit_context_class object.
!>
!>  Allocates memory and initializes a @ref v3fit_context_class object with an
!>  instance of a @ref commandline_parser object. This opens files for the
!>  recout, runlog, and result files. Over allocates space for the number of
!>  signals.
!>
!>  @returns A pointer to a constructed @ref v3fit_context_class object.
!-------------------------------------------------------------------------------
      FUNCTION v3fit_context_construct()

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), POINTER  :: v3fit_context_construct

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(v3fit_context_construct)

      CALL profiler_set_stop_time('v3fit_context_construct', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref v3fit_context_class object.
!>
!>  Deallocates memory and uninitializes a @ref v3fit_context_class object. This
!>  all the v3fit constructed objects are destroyed here.
!>
!>  @param[inout] this A @ref v3fit_context_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), POINTER :: this

!  Start of executable code
      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************

!-------------------------------------------------------------------------------
!>  @brief Create output files.
!>
!>  Creates output files that reconstructed results are written to.
!>
!>  @param[inout] this A @ref v3fit_context_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_create_files(this)
      USE safe_open_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: this

!  local variables
      INTEGER                                   :: status
      INTEGER                                   :: recon_rank
      REAL (rprec)                              :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

#if defined(MPI_OPT)
      CALL MPI_COMM_RANK(this%reconstruction_comm, recon_rank, status)
#else
      recon_rank = 0
#endif

!  Setup the runlog file.
!  Initalize a default value of the I\O unit. V3FIT increments from there.
      this%runlog_iou = 0
      CALL safe_open(this%runlog_iou, status,                                  &
     &               'runlog.vmec_interactive',                                &
     &               'replace', 'formatted', delim_in='none')
      CALL assert_eq(0, status, 'v3fit_context_construct: ' //                 &
     &               'Safe_open of runlog.vmec_interactive failed')

      WRITE (this%runlog_iou, *) 'VMEC_INTERACTIVE'

      CALL profiler_set_stop_time('v3fit_context_create_files',                &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Close output files.
!>
!>  Closes output files that reconstructed results are written to.
!>
!>  @param[inout] this A @ref v3fit_context_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_close_files(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: this

!  local variables
      INTEGER                                   :: status
      REAL (rprec)                              :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CLOSE(this%runlog_iou)

      CALL profiler_set_stop_time('v3fit_context_close_files',                 &
     &                            start_time)

      END SUBROUTINE

      END MODULE
