!*******************************************************************************
!>  @file bmw.f
!>  @brief Contains the main routines for Biot-Savart Magnetic Vmec
!>  Vector-Potential code.
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
!>  Highest level BMW routine.
!-------------------------------------------------------------------------------
      PROGRAM bmw
!$    USE omp_lib
      USE bmw_context
      USE bmw_commandline_parser
      USE, INTRINSIC :: iso_fortran_env, Only : output_unit

      IMPLICIT NONE

!  local variables
      TYPE (bmw_context_class), POINTER            :: context => null()
      TYPE (bmw_commandline_parser_class), POINTER ::                          &
     &   cl_parser => null()
      INTEGER                                      :: num_threads
      REAL (rprec)                                 :: start_time
!  Start of executable code
      CALL profiler_construct

      start_time = profiler_get_start_time()

      WRITE (*,*)
      WRITE (*,1000) series
      WRITE (*,*)

      cl_parser => bmw_commandline_parser_construct()

!  Configure the number of threeds to use.
!$    num_threads = bmw_commandline_parser_get_integer(cl_parser,              &
!$   &                                                 '-para', 1)

!$    IF (num_threads .gt. 0) THEN
!$       CALL OMP_SET_NUM_THREADS(num_threads)
!$    END IF
!$OMP PARALLEL
!$    num_threads = OMP_GET_MAX_THREADS()
!$OMP END PARALLEL
!$    WRITE (*,1001) num_threads

      context => bmw_context_construct(                                        &
     &   bmw_commandline_parser_get_string(cl_parser, '-mgridf'),              &
     &   bmw_commandline_parser_get_string(cl_parser, '-woutf'))
      CALL bmw_context_write(context,                                          &
     &        bmw_commandline_parser_get_string(cl_parser, '-outf'))

      CALL bmw_context_destruct(context)
      CALL bmw_commandline_parser_destruct(cl_parser)

      CALL profiler_set_stop_time('bmw_main', start_time)
      CALL profiler_write(output_unit)
      CALL profiler_destruct

1000  FORMAT('BMW ',i4,' Series.')
1001  FORMAT('Using ',i4,' threads.')

      END PROGRAM
