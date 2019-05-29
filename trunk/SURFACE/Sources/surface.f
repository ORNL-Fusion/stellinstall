!*******************************************************************************
!>  @file surface.f
!>  @brief Contains the main routines for SURFACE
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  SURFACE is a code for computing the virual surface on an arbitray wall
!>  structure. This uses the virtual casing to compute the interior magnetic
!>  field.
!>
!>  @author Mark Cianciosa
!*******************************************************************************
!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief SURFACE main program.
!>
!>  Highest level SURFACE routine.
!-------------------------------------------------------------------------------
      PROGRAM surface
      USE surface_context
      USE surface_commandline_parser
      USE bmw_parallel_context
      USE, INTRINSIC :: iso_fortran_env, Only : output_unit
      USE safe_open_mod

      IMPLICIT NONE

!  local variables
      TYPE (bmw_parallel_context_class), POINTER :: parallel => null()
      TYPE (surface_context_class), POINTER      :: context => null()
      TYPE (surface_commandline_parser_class), POINTER ::                      &
     &   cl_parser => null()
      INTEGER                                    :: io_unit
      INTEGER                                    :: status
      REAL (rprec)                               :: start_time

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

      cl_parser => surface_commandline_parser_construct(parallel)

!  Check if the required flags are set.
      IF (.not.surface_commandline_parser_is_flag_set(cl_parser,               &
     &                                                '-woutf')) THEN
         WRITE (*,1001) '-woutf'
         CALL surface_commandline_parser_print_help
      END IF
      IF (.not.surface_commandline_parser_is_flag_set(cl_parser,               &
     &                                                '-surff')) THEN
         WRITE (*,1001) '-surff'
         CALL surface_commandline_parser_print_help
      END IF
      IF (.not.surface_commandline_parser_is_flag_set(cl_parser,               &
     &    '-outf')) THEN
         WRITE (*,1001) '-outf'
         CALL surface_commandline_parser_print_help
      END IF

      CALL bmw_parallel_context_set_threads(parallel,                          &
     &        surface_commandline_parser_get_integer(cl_parser, '-para',       &
     &                                               1))

      io_unit = output_unit
      IF (surface_commandline_parser_is_flag_set(cl_parser,                    &
     &                                           '-logf')) THEN
         CALL safe_open(io_unit, status,                                       &
     &      surface_commandline_parser_get_string(cl_parser, '-logf'),         &
     &      'replace', 'formatted', delim_in='none')
      END IF


      IF (parallel%offset .eq. 0) THEN
         WRITE (io_unit,*)
         WRITE (io_unit,1000) 1
         WRITE (io_unit,*)
      END IF
      CALL bmw_parallel_context_report(parallel, io_unit)

      context => surface_context_construct(                                    &
     &   surface_commandline_parser_get_string(cl_parser, '-woutf'),           &
     &   surface_commandline_parser_get_string(cl_parser, '-surff'),           &
     &   parallel, io_unit)
      CALL surface_context_set_current(context, io_unit, parallel)
      CALL surface_context_write(context,                                      &
     &        surface_commandline_parser_get_string(cl_parser,                 &
     &                                              '-outf'),                  &
     &        io_unit, parallel)

      CALL surface_context_destruct(context)
      CALL surface_commandline_parser_destruct(cl_parser)

      CALL profiler_set_stop_time('surface_main', start_time)
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

1000  FORMAT('SURFACE verson ',i4,'.')
1001  FORMAT('Required flag ',a,' not set.')

      END PROGRAM
