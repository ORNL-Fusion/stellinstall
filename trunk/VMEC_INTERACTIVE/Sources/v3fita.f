      module vmec_interactive
      
      use v3_utilities
      use stel_kinds
      use profiler
      USE vmec_equilibrium
      USE safe_open_mod
      USE file_opts

      implicit none
      integer                             :: i ! iteration / call number
      integer                             :: am_id
      TYPE (vmec_class), POINTER          :: vmec => null()
      REAL (rprec)                        :: time_start, time_end
      INTEGER                             :: eq_steps
      INTEGER                             :: restart_step
      INTEGER                             :: error
      LOGICAL                             :: eq_converged
      LOGICAL                             :: write_input
      REAL (rprec)                        :: start_time
      integer                             :: num_iter = 0
      logical                             :: is_changed = .false.
      character(len=*), parameter         :: vmec_nli_filename =        &
     &                               'input.vmec.test'
      
      INTEGER                             :: status
      INTEGER                             :: runlog_iou

!--- MPI -----------------------------------------------------------------------
      INTEGER                             :: mpi_rank = 0
#if defined (MPI_OPT)
      INTEGER                             :: mpi_size
      integer                             :: mpi_comm
#else
      integer, parameter                  :: mpi_comm = 0
#endif
      
      contains
      
! setup of Vmec Interactive
#if defined(MPI_OPT)
      subroutine vi_setup(my_mpi_comm) bind (C, name="vi_setup")
#else
      subroutine vi_setup() bind (C, name="vi_setup")
#endif

      use iso_c_binding
      
      implicit none
      
#if defined(MPI_OPT)
      integer, intent(in) :: my_mpi_comm
      
!  first executable statement for MPI case is here  !
      mpi_comm = my_mpi_comm
#endif
      
      
      CALL second0(time_start)
      CALL profiler_construct

      ! first iteration
      i = 0
      
#if defined (MPI_OPT)
      CALL MPI_COMM_SIZE(mpi_comm, mpi_size, error)
      CALL MPI_COMM_RANK(mpi_comm, mpi_rank, error)
#endif

      IF (mpi_rank .eq. 0) THEN
!   open runlog_iou for logging
!  Initalize a default value of the I\O unit. V3FIT increments from there.
         runlog_iou = 0
         CALL safe_open(runlog_iou, status,                             &
     &                  'runlog.vmec_interactive',                      &
     &                  'replace', 'formatted', delim_in='none')
         CALL assert_eq(0, status, 'v3fit_context_construct: ' //       &
     &                  'Safe_open of runlog.vmec_interactive failed')

         WRITE (runlog_iou, *) 'VMEC_INTERACTIVE'
      END IF
      write_input = .true.
      vmec => vmec_construct_eq(vmec_nli_filename,                      &
     &                          runlog_iou,                             &
     &                          mpi_comm)
     
      IF (mpi_rank .eq. 0) THEN
         write(*,*) 'vmec wout filename: "',                            &
     &              trim(vmec%wout_file_name),'"'

         WRITE (*,*) ' *** Starting initial equilibrium ' //            &
     &               'convergence (rec)'
         WRITE (runlog_iou,*) ' *** Starting initial ' //               &
     &                                'equilibrium convergence'
      end if
      
#if defined(MPI_OPT)
      CALL MPI_BCAST(num_iter, 1, MPI_INTEGER, 0,                       &
     &               mpi_comm, error)
#endif
      
      eq_converged = vmec_converge(vmec, num_iter,                      &
     &                             runlog_iou,                          &
     &                             mpi_comm)
     
      CALL assert(eq_converged,                                         &
     &               'task_reconstruct no inital convergence')
     
      if (mpi_rank .eq. 0) then
         WRITE (*,*)
         WRITE (*,1002) 0

         WRITE (runlog_iou,*)
         WRITE (runlog_iou,1002) 0
      
         IF (write_input) THEN
            CALL vmec_write_input(vmec, 0)
         END IF
      end if
      
1002  FORMAT('  *** Reconstruction step ',i4)
      
      end subroutine vi_setup

      
      
! set an input variable value
      subroutine vi_set(centralPressure) bind (C, name="vi_set")
      
      use iso_c_binding
      
      implicit none
      
      real(C_DOUBLE), intent(in) :: centralPressure

!  get ID of parameter to change      
      am_id = vmec_get_param_id(vmec, 'am')
      
!   change parameter: central pressure
      call vmec_set_param(vmec, am_id, 0, 0, centralPressure,           &
     &                    mpi_comm)

      call vmec_set_param(vmec, am_id, 1, 0, -centralPressure,          &
     &                    mpi_comm)

      if (mpi_rank .eq. 0) then
         write(*,*) 'central pressure: ',centralPressure,' Pa'
      end if
      
      is_changed = .true.
      
      end subroutine vi_set
      
      
      
      
      
! get an input or output variable value from VMEC
      subroutine vi_get(error) bind (C, name="vi_get")
      
      use iso_c_binding
      
      implicit none
      
      integer, intent(inout) :: error
      
      if (is_changed) then
         i = i+1
      
         if (mpi_rank .eq. 0) then
            WRITE (*,*)
            WRITE (*,1002) i

            WRITE (runlog_iou,*)
            WRITE (runlog_iou,1002) i
         end if


         
!   converge vmec
         IF (vmec_converge(vmec, num_iter, runlog_iou,                  &
     &                     mpi_comm)) THEN

            if (mpi_rank .eq. 0) then
               write(*,*) 'equilibrium converged :-)'
               IF (write_input) THEN
                  CALL vmec_write_input(vmec, i)
               END IF
            end if
            error = 0
         ELSE
            if (mpi_rank .eq. 0) then
               write(*,*) 'equilibrium not converged :-('
            end if
            error = 1
         END IF
      end if
      
1002  FORMAT('  *** Reconstruction step ',i4)

      end subroutine vi_get
      
!  controlled shutdown of the interactive VMEC service
      subroutine vi_exit() bind (C, name="vi_exit")
      
      use iso_c_binding
      
      implicit none
      
!  Clean up memory
      IF (mpi_rank .eq. 0) THEN
         CALL second0(time_end)

         WRITE (runlog_iou,1000) time_end - time_start
         WRITE (*,1000) time_end - time_start

         CALL profiler_write(runlog_iou)

         CLOSE(runlog_iou)
      END IF
      CALL profiler_destruct

1000  FORMAT('  *** V3FIT run time (s): ',f10.2)

      end subroutine vi_exit

      end module
