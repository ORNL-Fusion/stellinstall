      PROGRAM v3fit
      USE v3fit_context
      USE profiler
      use vmec_equilibrium

      IMPLICIT NONE

!  local variables
      TYPE (v3fit_context_class), POINTER :: context => null()
      TYPE (vmec_class), POINTER          :: vmec => null()
      REAL (rprec)                        :: time_start, time_end
      INTEGER                             :: eq_rank
      INTEGER                             :: recon_rank
      
      INTEGER                             :: eq_steps
      INTEGER                             :: restart_step
      INTEGER                             :: i, error
      LOGICAL                             :: eq_converged
      LOGICAL                             :: write_input
      REAL (rprec)                        :: start_time
      integer                             :: num_iter = 0
      integer                             :: am_id
      integer ::  state_flags = 0
      
      character(len=*), parameter :: vmec_nli_filename =                       &
     &                               'input.vmec.test'
      
!  Start of executable code

#if defined(MPI_OPT)
!  The Intel version of the MPI libraries does not provide a correct value for
!  time until after MPI_INIT is called. Make sure this is the first thing called
!  so that correct timing information can be used.
      CALL MPI_INIT(error)
#endif

      CALL second0(time_start)
      CALL profiler_construct

!  Create a context and read in the v3fit namelist input file.
      context => v3fit_context_construct()
      
!  Configure parallel computation.
      CALL config_parallelism(context)

      eq_rank = 0
#if defined(MPI_OPT)
      CALL MPI_COMM_RANK(context%equilibrium_comm, eq_rank, error)
#endif
      
      IF (eq_rank .eq. 0) THEN
!   open runlog_iou for logging
         CALL v3fit_context_create_files(context)
      END IF
      
      
!  use the reconstruction task as a testbed for re-converging VMEC
      write_input = .true.

      vmec => vmec_construct_eq(vmec_nli_filename,                             &
     &                          context%runlog_iou,                            &
#if defined(MPI_OPT)
     &                       context%equilibrium_comm,                      &
     &                       context%reconstruction_comm)
#else
     &                       0, 0)
#endif
      
      write(*,*) 'vmec wout filename: "',trim(vmec%wout_file_name),'"'
      
      
!  Solve equilibrium
      WRITE (*,*) ' *** Starting initial equilibrium ' //                   &
     &               'convergence (rec)'
      WRITE (context%runlog_iou,*) ' *** Starting initial ' //              &
     &                             'equilibrium convergence'
      
#if defined (MPI_OPT)
      CALL MPI_BCAST(mpi_equilibrium_task, 1, MPI_INTEGER, 0,               &
     &               context%equilibrium_comm, error)

      CALL MPI_BCAST(num_iter, 1, MPI_INTEGER, 0,                           &
     &               context%equilibrium_comm, error)
#endif

      eq_steps = 1
      eq_converged = vmec_converge(vmec, num_iter,                          &
     &                             context%runlog_iou,                      &
#if defined(MPI_OPT)
     &                             context%equilibrium_comm)
#else
     &                             0)
#endif
      
      
      CALL assert(eq_converged,                                             &
     &               'task_reconstruct no inital convergence')
     
     

#if defined(MPI_OPT)
!  Sync the number of equilibrium steps the inital convergence took so that the'
!  child processes start in the same place.
      CALL MPI_BCAST(eq_steps, 1, MPI_INTEGER, 0,                           &
     &               context%reconstruction_comm, error)
#endif

      WRITE (*,*)
      WRITE (*,1002) 0

      WRITE (context%runlog_iou,*)
      WRITE (context%runlog_iou,1002) 0
      
      IF (write_input) THEN
         CALL vmec_write_input(vmec, 0)
      END IF

      am_id = vmec_get_param_id(vmec, 'am')
      
      
      
      restart_step = 1
      DO i = restart_step, 10
         WRITE (*,*)
         WRITE (*,1002) i

         WRITE (context%runlog_iou,*)
         WRITE (context%runlog_iou,1002) i

!   change parameter: central pressure
         call vmec_set_param(vmec, am_id, 0, 0, real(i,8)*1000.0,              &
#if defined(MPI_OPT)
     &                           context%equilibrium_comm, state_flags)
#else
     &                           0, state_flags)
#endif
         
         call vmec_set_param(vmec, am_id, 1, 0, -real(i,8)*1000.0,             &
#if defined(MPI_OPT)
     &                           context%equilibrium_comm, state_flags)
#else
     &                           0, state_flags)
#endif

         
         write(*,*) 'central pressure: ',real(i,8),' kPa'
         
!   converge vmec
         IF (vmec_converge(vmec, num_iter, context%runlog_iou,                 &
#if defined(MPI_OPT)
     &                           context%equilibrium_comm)) THEN
#else
     &                           0)) THEN
#endif
            write(*,*) 'equilibrium converged :-)'
!  Step was successful. Post process the result.

            IF (write_input) THEN
               CALL vmec_write_input(vmec, i)
            END IF
         ELSE
            write(*,*) 'equilibrium not converged :-('
!  Step was not successful. Exit the loop.
            EXIT
         END IF
      END DO

#if defined(MPI_OPT)
      CALL MPI_BCAST(mpi_quit_task, 1, MPI_INTEGER, 0,                         &
     &               context%reconstruction_comm, error)
      CALL MPI_BCAST(mpi_quit_task, 1, MPI_INTEGER, 0,                         &
     &               context%equilibrium_comm, error)
#endif

      CALL profiler_set_stop_time('task_reconstruct', start_time)

      
      
      
      
      
      
!  Clean up memory
      IF (eq_rank .eq. 0) THEN
         recon_rank = 0
#if defined(MPI_OPT)
         CALL MPI_COMM_RANK(context%reconstruction_comm, recon_rank,           &
     &                      error)
#endif

         IF (recon_rank .eq. 0) THEN
            CALL second0(time_end)

            WRITE (context%runlog_iou,1000) time_end - time_start
            WRITE (*,1000) time_end - time_start

            CALL profiler_write(context%runlog_iou)

         END IF

         CALL v3fit_context_close_files(context)
      END IF
      CALL profiler_destruct

      CALL v3fit_context_destruct(context)

#if defined(MPI_OPT)
      CALL MPI_FINALIZE(error)
#endif

1000  FORMAT('  *** V3FIT run time (s): ',f10.2)
1001  FORMAT('nrstep 'i4,' completed. May be reconstructed.')
1002  FORMAT('  *** Reconstruction step ',i4)
1003  FORMAT('Initial g^2 = ',es12.5)

      END PROGRAM
      
!*******************************************************************************
!  INITIALIZE SUBROUTINES
!*******************************************************************************

!-------------------------------------------------------------------------------
!>  @brief Configure Parallelism
!>
!>  When using OpenMP, this sets the maximum number of threads to use. Threads
!>  operate within a process and share memory with the main thread. With in
!>  parallel constructs, it's important to ensure that all code is thread safe
!>  to avoid unexpected behavior and race conditions. The maximum number of
!>  threads is controled using the @fixed_width{-para} command line argument. An
!>  argument of a negative values specifies to use the default behavior.
!>
!>  When using MPI, this sets and configures the processes. A unique directory
!>  needs to be generated for each process to avoid overwriting output files.
!>
!>  The default way is to create and distribute the threads externally. In this
!>  case, the number of processes is fixed. Each process will need to be
!>  reassigned a working directory. One process will be used like a parent
!>  process. The remaining processed will be retasked as child processes.
!>
!>  The second way is to spawn sub processes from a single parent process. In
!>  this case, the parent process will spawn subprocesses to operate the child
!>  task in a sub directory. The number of sub processes can be configured
!>  automatically or manually. @see cl_parsing_sec
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
!  When using openmp the !$ allows conditional computation based on if openmp is
!  used or not.
      SUBROUTINE config_parallelism(context)
      USE v3fit_context
      USE file_opts
#if defined(MPI_OPT)
      USE data_parameters
#endif

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      INTEGER                                   :: num_threads
#if defined(MPI_OPT)
      INTEGER                                   :: error
      INTEGER                                   :: num_recon_processes
      INTEGER                                   :: num_eq_processes
      INTEGER                                   :: i
      INTEGER                                   :: mpi_rank
      INTEGER                                   :: mpi_size
      INTEGER                                   :: eq_rank
      INTEGER                                   :: eq_size
      INTEGER                                   :: recon_rank
      CHARACTER (len=dir_prefix_len)            :: directory_name
      CHARACTER (len=path_length)               :: temp_string
#endif
      REAL (rprec)                              :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (*,*) ' *** Configuring parallelism.'

!-- MPI ------------------------------------------------------------------------

#if defined (MPI_OPT)
!  Configure num processes to use for MPI
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, error)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, error)

!  Initialize the intra communicator to the world communicator. This is done so
!  global_comm has is value even when num_processes is 1 or the processes were
!  created externally.
      context%global_comm = MPI_COMM_WORLD

!  MPI parallelization is only implemented for specific tasks. If the root task
!  is not a reconstruction task stop the child processes.
      num_recon_processes = 1
      num_eq_processes = mpi_size

!  Wait for all processes to reach this point.
      CALL MPI_BARRIER(context%global_comm, error)

!  Partition the equilibrium communicators.
      recon_rank = MOD(mpi_rank, MAX(1, num_recon_processes))
      CALL MPI_COMM_SPLIT(context%global_comm, recon_rank, 0,                  &
     &                    context%equilibrium_comm, error)
      CALL MPI_COMM_RANK(context%equilibrium_comm, eq_rank, error)

!  Take the root process of each equilibrium communicator and partition them
!  from global to the reconstuction communicator.
      IF (eq_rank .gt. 0) THEN
         CALL MPI_COMM_SPLIT(context%global_comm, MPI_UNDEFINED, 0,            &
     &                       context%reconstruction_comm, error)
      ELSE
         CALL MPI_COMM_SPLIT(context%global_comm, 0, recon_rank,               &
     &                       context%reconstruction_comm, error)
      END IF

      CALL MPI_BARRIER(context%global_comm, error)
#endif

      CALL profiler_set_stop_time('config_parallelism', start_time)

1000  FORMAT('Using ',i4,' threads.')
1001  FORMAT('Configured global process ',i5,', recon process ',i5,            &
     &       ', eq process ',i5,' task ',a)

      END SUBROUTINE

