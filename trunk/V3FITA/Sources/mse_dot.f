!-------------------------------------------------------------------------------
!  The @fixed_width, @begin_table, @item2 and @end_table commands are custom
!  defined commands in Doxygen.in. They are defined under ALIASES. For the page
!  created here, the 80 column limit is exceeded. Arguments of aliases are
!  separated by ','. If you intended ',' to be a string you must use an escaped
!  comma '\,'.
!
!>  @page mse_dot_sec Motional Stark Effect Diagnostic Dot File
!>
!>  @tableofcontents
!>  @section mse_dot_intro_sec Introduction
!>  Diagnostic dot files specify geometric and auxiliary information about each
!>  diagnostic. Parsing of a diagnostic dot file begins by searching for a valid
!>  keyword. All lines before a valid keyword is reached are ignored allowing
!>  for comments and other information to be inserted in the file. Once a
!>  keyword is reached, the parsing becomes strict. The diagnostic description
!>  format is dependant on the specific keyword. All diagnostic dot files must
!>  end with an @fixed_width{end_of_file} keyword.
!>
!>  This document contains information about the motional stark effect point
!>  chord specifications. Points are specified in a structured text file. This
!>  file is read in by the @ref v3fit_input::mse_dot_filename variable of the
!>  name list input file. MSE points are specified as a point, two vertices in
!>  defining the direction of the viewing cord and two vertices defining the
!>  direction of the neutral beam in either cartesian or cyclindrical
!>  coordinates.
!>
!>  @note The point of the MSE measurement could be infered from the
!>  intersection of the neutral beam vector and the viewing vector. However for
!>  numercial reasons the intersection point is expicity defined.
!>  @note The neutral beam vector could have been defined once however it's
!>  simpler to define it for each measurement point.
!>
!>  @section mse_dot_key_sec Keywords
!>  @begin_table
!>  @item2{@fixed_width{mse_point_XYZ_deg},      A MSE point specified in Cartesian Coordinates outputting the signal in degrees.}
!>  @item2{@fixed_width{mse_point_RPhiDegZ_deg}, A MSE point specified in Cylindrical Coordinates outputting the signal in degrees.}
!>  @item2{@fixed_width{mse_point_XYZ_rad},      A MSE point specified in Cartesian Coordinates outputting the signal in radians.}
!>  @item2{@fixed_width{mse_point_RPhiDegZ_rad}, A MSE point specified in Cylindrical Coordinates outputting the signal in radians.}
!>  @end_table
!>
!>  @section mse_dot_point_spec_sec Point Specification
!>  All motional stark effect points are specified in the following manner.
!>
!>  @fixed_width{key_word\n}
!>  @fixed_width{name(STRING len:@ref data_parameters::data_short_name_length)\n}
!>  @fixed_width{x(REAL) y(REAL) z(REAL)\n}
!>  @fixed_width{vx1(REAL) vy1(REAL) vz1(REAL)\n}
!>  @fixed_width{vx2(REAL) vy2(REAL) vz2(REAL)\n}
!>  @fixed_width{nx1(REAL) ny1(REAL) nz1(REAL)\n}
!>  @fixed_width{nx2(REAL) ny2(REAL) nz2(REAL)\n}
!>
!>  @section mse_dot_exam_sec Example File
!>  @code
!>  !  A motional stark effect chord specified in Cylindrical coordinates.
!>  mse_point_XYZ_deg
!>  mse001
!>  0.0 0.5 0.0
!>  0.0 0.0 0.0
!>  0.0 1.0 0.0
!>  -1.0 0.5 0.0
!>  1.0 0.5 0.0
!>
!>  !  All files must end with
!>  end_of_file
!>  @endcode
!>
!>  @section mse_dot_prog_ref_sec Programmers Reference
!>  Reference material for the coding to parse these files can be found in the
!>  @ref signal_dot and @ref mse_dot modules.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file mse_dot.f
!>  @brief Contains module @ref mse_dot
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module for opening and reading a 'mse.' file. The file format for these
!>  files are documented in @ref mse_dot_sec.
!*******************************************************************************

      MODULE mse_dot

      USE stel_kinds
      USE stel_constants
      USE signal_dot

      IMPLICIT NONE

      CONTAINS
!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Read an motional stark effect diagnostic dot file.
!>
!>  Parses the structure of a dignostic dot file and constructs a @ref mse
!>  object. Each point is parsed by @ref mse_dot_parse_point.
!>
!>  @param[in]    mse_file       Name of the mse. file.
!>  @param[inout] signals         Array of pointers to the a signal
!>  @param[inout] signals_created Number valid signals in the array, or index of
!>                                the last signal created.
!>  @param[in]    observed        Array of observed values.
!>  @param[in]    sigma           Array of observed sigma values
!>  @param[in]    weight          Array of signal weight values
!>  @param[inout] first_index     Index of the first instance of an intpol
!>                                signal.
!-------------------------------------------------------------------------------
      SUBROUTINE mse_dot_read(mse_file, signals, signals_created,              &
     &                        observed, sigma, weight, first_index)
      USE signal
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_filename_length), INTENT(in)   :: mse_file
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout) :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in)             :: observed
      REAL (rprec), DIMENSION(:), INTENT(in)             :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in)             :: weight
      INTEGER, INTENT(inout)                             :: first_index

!  local parameters
!  n_mse_keyword        integer - number of mse_keyword s.
      INTEGER, PARAMETER :: n_mse_keyword = 5

!  local variables
!  mse_keyword    character array - keywords for various ipch input types
      CHARACTER(len=data_name_length), DIMENSION(1:n_mse_keyword) ::           &
     &   mse_keyword
      TYPE (signal_dot_file) :: mse_dot_file_ref
      REAL (rprec)           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Initialize the keywords and keyword lengths
      mse_keyword(1) = 'mse_point_XYZ_deg'
      mse_keyword(2) = 'mse_point_RPhiDegZ_deg'
      mse_keyword(3) = 'mse_point_XYZ_rad'
      mse_keyword(4) = 'mse_point_RPhiDegZ_rad'
      mse_keyword(5) = 'end_of_file'

!  Open up the 'mse.' file
      mse_dot_file_ref = signal_dot_open(TRIM(mse_file), 'ipch')

!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(mse_dot_file_ref,                &
     &                                        mse_keyword))

            CASE DEFAULT ! This case should never fire.
               EXIT ! Exit out of infinte loop.

            CASE ('end_of_file')
               EXIT ! Exit out of infinte loop.

            CASE ('mse_point_XYZ_deg')

            CASE ('mse_point_RPhiDegZ_deg')

            CASE ('mse_point_XYZ_rad')

            CASE ('mse_point_RPhiDegZ_rad')

         END SELECT

      END DO

!  Close the 'mse.' file
      CALL signal_dot_close(mse_dot_file_ref)

      CALL profiler_set_stop_time('mse_dot_read', start_time)

      END SUBROUTINE

!*******************************************************************************
!  PARSING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Parse motional stark effect diagnostic point.
!>
!>  Parses a signle point. Points are laid out as,
!>
!>  * chord_name
!>  * Point of measurement as three reals.
!>  * Starting point of view as three reals.
!>  * Ending point of view as three reals.
!>  * Starting point of beam as three reals.
!>  * Ending point of beam as three reals.
!>
!>  @param[in]    mse_dot_file_ref Reference to the mse. file.
!>  @param[in]    coordinate_type  Specifies the cartesian or cylindrical
!>                                 coordinate type.
!>  @param[inout] signals          Array of pointers to the a signal
!>  @param[inout] signals_created  Number valid signals in the array, or index
!>                                 of the last signal created.
!>  @param[in]    observed         Array of observed values.
!>  @param[in]    sigma            Array of observed sigma values
!>  @param[in]    weight           Array of signal weight values
!>  @param[inout] first_index      Index of the first instance of an intpol
!>                                 signal.
!>  @param[in]    in_degrees       Signal measures degrees if true.
!-------------------------------------------------------------------------------
      SUBROUTINE mse_dot_parse_point(mse_dot_file_ref, coordinate_type,        &
     &                               signals, signals_created,                 &
     &                               observed, sigma, weight,                  &
     &                               current_profile, first_index,             &
     &                               in_degrees)
      USE signal
      USE v3fit_input, only: v3fit_max_diagnostics
      USE v3_utilities, only: err_fatal
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_dot_file), INTENT(inout)  :: mse_dot_file_ref
      CHARACTER (len=*), INTENT(in)          :: coordinate_type
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout)                 :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in) :: observed
      REAL (rprec), DIMENSION(:), INTENT(in) :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in) :: weight
      INTEGER, INTENT(in)                    :: current_profile
      INTEGER, INTENT(inout)                 :: first_index
      LOGICAL, INTENT(in)                    :: in_degrees

!  local variables
      REAL (rprec), DIMENSION(3)             :: xcart
      REAL (rprec), DIMENSION(3)             :: view_start
      REAL (rprec), DIMENSION(3)             :: view_end
      REAL (rprec), DIMENSION(3)             :: beam_start
      REAL (rprec), DIMENSION(3)             :: beam_end
      CHARACTER (len=data_short_name_length) :: point_name
      CHARACTER (len=signal_dot_line_len)    :: line
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Get the position and name of the point.
      CALL signal_dot_parse_chord(mse_dot_file_ref, coordinate_type,           &
     &                            point_name, xcart)

!  Get the start and end positons of the viewing chord.
      view_start = signal_dot_parse_3_real(mse_dot_file_ref,                   &
     &                                     'Failed to parse view ' //          &
     &                                     'start position.')
      view_end = signal_dot_parse_3_real(mse_dot_file_ref,                     &
     &                                   'Failed to parse view ' //            &
     &                                   'end position.')

!  Get the start and end positons of the beam chord.
      view_start = signal_dot_parse_3_real(mse_dot_file_ref,                   &
     &                                     'Failed to parse beam ' //          &
     &                                     'start position.')
      view_end = signal_dot_parse_3_real(mse_dot_file_ref,                     &
     &                                   'Failed to parse beam ' //            &
     &                                   'end position.')

!  Convert Coordinates if necessay.
      IF (TRIM(coordinate_type) .eq. 'RPHiDegZ') THEN
		 view_start(2) = view_start(2)*degree ! Convert from degrees to radians.
         view_start = cyl_to_cart(view_start)
         view_end(2) = view_end(2)*degree ! Convert from degrees to radians.
         view_end = cyl_to_cart(view_end)

         beam_start(2) = beam_start(2)*degree ! Convert from degrees to radians.
         beam_start = cyl_to_cart(beam_start)
         beam_end(2) = beam_end(2)*degree ! Convert from degrees to radians.
         view_end = cyl_to_cart(beam_end)
      END IF

      signals(signals_created + 1)%p => signal_construct(                      &
     &   diagnostic_construct(                                                 &
     &      mse_construct(xcart, view_start, view_end,                         &
     &                    beam_start, beam_end, in_degrees)),                  &
     &   point_name, point_name, 'arb', observed(signals_created + 1),         &
     &   sigma(signals_created + 1), weight(signals_created + 1))
      signals_created = signals_created + 1

!  At lease one mse signal was made. Set the first index. This should only be
!  run once.
      IF (first_index .eq. -1) THEN
         first_index = signals_created
      END IF

      CALL profiler_set_stop_time('mse_dot_parse_point', start_time)

      END SUBROUTINE

      END MODULE mse_dot
