!*******************************************************************************
!>  @file bmw.f
!>  @brief Contains the main routines for Biot-Savart Magnetic Vmec
!>  Vector-Potential (BMW) code.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  BMW is a code for extending fields belond the VMEC domain in a manner that
!>  ensures divergence free fields. BMW does a Biot-Savart volume integration of
!>  of the equilibrium current density to obtain a continous vector potential
!>  every where on the mgrid grid.
!>
!>  @author Mark Cianciosa
!>
!>  Below is a brief discription of the major top level objects of the code. For
!>  discriptions of lower level objects consult the referenced top level
!>  objects.
!>
!>  @ref m_grid        Object containing the vaccum field information.
!>  @ref primed_grid   Object containing the plasma currents and primed grid
!>                     positions.
!>  @ref unprimed_grid Object containing the plasma vector potential response.
!>  @ref bmw_context   Defines the main pointers to all allocated memory and
!>                     objects.
!*******************************************************************************
!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief BMW main program.
!>
!>  Highest level BMW routine. This computes the vector potential on the mgrid
!>  grid.
!-------------------------------------------------------------------------------
      PROGRAM bmw
      USE bmw_context
      USE bmw_commandline_parser
      USE, INTRINSIC :: iso_fortran_env, Only : output_unit
      USE bmw_state_flags
      USE safe_open_mod

      IMPLICIT NONE

!  local variables
      TYPE (bmw_parallel_context_class), POINTER   :: parallel => null()
      TYPE (bmw_context_class), POINTER            :: context => null()
      TYPE (bmw_commandline_parser_class), POINTER ::                          &
     &   cl_parser => null()
      INTEGER                                      :: flags
      INTEGER                                      :: num_p
      INTEGER                                      :: io_unit
      INTEGER                                      :: status
      REAL (rprec)                                 :: start_time

!  Start of executable code
#if defined (MPI_OPT)
      CALL MPI_INIT(status)
#endif
      CALL profiler_construct

      start_time = profiler_get_start_time()

      parallel => bmw_parallel_context_construct(                              &
#if defined (MPI_OPT)
     &                                           MPI_COMM_WORLD                &
#endif
     &                                          )

      cl_parser => bmw_commandline_parser_construct(parallel)

!  Check if the required flags are set.
      IF (.not.bmw_commandline_parser_is_flag_set(cl_parser,                   &
     &    '-mgridf')) THEN
         WRITE (*,1001) '-mgridf'
         CALL bmw_commandline_parser_print_help
      END IF
      IF (.not.bmw_commandline_parser_is_flag_set(cl_parser,                   &
     &    '-woutf')) THEN
         WRITE (*,1001) '-woutf'
         CALL bmw_commandline_parser_print_help
      END IF
      IF (.not.bmw_commandline_parser_is_flag_set(cl_parser,                   &
     &    '-outf')) THEN
         WRITE (*,1001) '-outf'
         CALL bmw_commandline_parser_print_help
      END IF

      CALL bmw_parallel_context_set_threads(parallel,                          &
     &        bmw_commandline_parser_get_integer(cl_parser, '-para', 1))

      flags = bmw_state_flags_off
      IF (bmw_commandline_parser_is_flag_set(cl_parser, '-force')) THEN
         flags = IBSET(flags, bmw_state_flags_force)
      END IF
      IF (bmw_commandline_parser_is_flag_set(cl_parser, '-ju')) THEN
         flags = IBSET(flags, bmw_state_flags_ju)
      END IF
      IF (bmw_commandline_parser_is_flag_set(cl_parser, '-jv')) THEN
         flags = IBSET(flags, bmw_state_flags_jv)
         flags = IBCLR(flags, bmw_state_flags_ju)
      END IF
      IF (bmw_commandline_parser_is_flag_set(cl_parser,                        &
     &                                       '-siestaf')) THEN
         flags = IBSET(flags, bmw_state_flags_siesta)
         flags = IBCLR(flags, bmw_state_flags_ju)
         flags = IBCLR(flags, bmw_state_flags_jv)
      END IF

      flags = IBSET(flags, bmw_state_flags_mgrid)
      num_p = 0

      io_unit = output_unit
      IF (bmw_commandline_parser_is_flag_set(cl_parser, '-logf')) THEN
         CALL safe_open(io_unit, status,                                       &
     &      bmw_commandline_parser_get_string(cl_parser, '-logf'),             &
     &      'replace', 'formatted', delim_in='none')
      END IF

      IF (parallel%offset .eq. 0) THEN
         WRITE (io_unit,*)
         WRITE (io_unit,1000) series
         WRITE (io_unit,*)
      END IF
      CALL bmw_parallel_context_report(parallel, io_unit)

      context => bmw_context_construct(                                        &
     &   bmw_commandline_parser_get_string(cl_parser, '-mgridf'),              &
     &   bmw_commandline_parser_get_string(cl_parser, '-woutf'),               &
     &   bmw_commandline_parser_get_string(cl_parser, '-siestaf'),             &
     &   flags, num_p, parallel, io_unit)
      CALL bmw_context_set_up_grid(context,                                    &
     &   bmw_commandline_parser_get_integer(cl_parser, '-p_start', -1),        &
     &   bmw_commandline_parser_get_integer(cl_parser, '-p_end', -1),          &
     &   parallel, io_unit)
      CALL bmw_context_write(context,                                          &
     &        bmw_commandline_parser_get_string(cl_parser, '-outf'),           &
     &        parallel)

      CALL bmw_context_destruct(context)
      CALL bmw_commandline_parser_destruct(cl_parser)

      CALL profiler_set_stop_time('bmw_main', start_time)
      IF (parallel%offset .eq. 0) THEN
         CALL profiler_write(io_unit)
      END IF
      CALL profiler_destruct

      CLOSE(io_unit)

      CALL bmw_parallel_context_destruct(parallel                              &
#if defined (MPI_OPT)
     &                                   ,.true.                               &
#endif
     &                                  )

1000  FORMAT('BMW ',i4,' Series.')
1001  FORMAT('Required flag ',a,' not set.')

      END PROGRAM
