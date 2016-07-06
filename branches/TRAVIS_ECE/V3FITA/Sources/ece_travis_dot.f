!-------------------------------------------------------------------------------
!  The @fixed_width, @begin_table, @item2 and @end_table commands are custom
!  defined commands in Doxygen.in. They are defined under ALIASES. For the page
!  created here, the 80 column limit is exceeded. Arguments of aliases are
!  separated by ','. If you intended ',' to be a string you must use an escaped
!  comma '\,'.
!
!>  @page ece_travis_dot_sec ECE TRAVIS Diagnostic Dot File
!>
!>  @tableofcontents
!>  @section ece_dot_intro_sec Introduction
!>  Diagnostic dot files specify geometric and auxiliary information about each
!>  diagnostic. Parsing of a diagnostic dot file begins by searching for a valid
!>  keyword. All lines before a valid keyword is reached are ignored allowing
!>  for comments and other information to be inserted in the file. Once a
!>  keyword is reached, the parsing becomes strict. The diagnostic description
!>  format is dependant on the specific keyword. All diagnostic dot files must
!>  end with an @fixed_width{end_of_file} keyword.
!>
!>  TRAVIS is an external Ray-tracing code for electron cyclotron calculations. 
!>  [N.B. Marushchenko et al., CPC 185(2015)]. This document contains information
!>  about the specification of ECE diagnostics that use TRAVIS to calculate the
!>  ECE signal.      
!>
!>  The ECE signals are specified in a structured text file. This file is read in by the
!>  @ref v3fit_input::ece_travis_dot_filename variable of the name list input file.
!>
!>  @section ece_travis_dot_key_sec Keywords
!>  @begin_table
!>  @item2{@fixed_width{ece_travis_channels}, A set of ECE channels specified by a TRAVIS ece_data file.}
!>  @end_table
!>
!>  @section ece_travis_dot_sig_spec_sec Channel Specification
!>  All TRAVIS ECE channels are specified in the following manner.
!>
!>  @fixed_width{key_word\n}
!>  @fixed_width{base_name(STRING len:@ref data_parameters::data_short_name_length)\n}
!>  @fixed_width{TRAVIS_data_file(STRING len:@ref file_opts::path_length)\n}
!>  @fixed_width{channels(INTEGER)\n}
!>
!> This will create a signal for each channel in the TRAVIS data file.
!> Channels must match the number of channels in the TRAVIS data file
!>
!>  @section ece_travis_dot_exam_sec Example File
!>  @code
!>  !  A set of ECE channels.
!>  ece_travis_channels
!>  eceray
!>  travisinput_ece.data
!>  32
!>      
!>  !  All files must end with
!>  end_of_file
!>  @endcode
!>
!>  @section ece_travis_dot_prog_ref_sec Programmers Reference
!>  Reference material for the coding to parse these files can be found in the
!>  @ref signal_dot and @ref ece_travis_dot modules.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file ece_travis_dot.f
!>  @brief Contains module @ref ece_dot
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module for opening and reading a 'ece_travis.' file. The file format for these
!>  files are documented in @ref ece_travis_dot_sec.
!*******************************************************************************

      MODULE ece_travis_dot

      USE stel_kinds
      USE stel_constants
      USE signal_dot

      
      IMPLICIT NONE

      CONTAINS
!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Read an ECE diagnostic dot file.
!>
!>  Parses the structure of a dignostic dot file and constructs a @ref ece
!>  object. Each chord is parsed by @ref ece_dot_parse_chord.
!>
!>  @param[in]    ece_file        Name of the ece. file.
!>  @param[inout] signals         Array of pointers to the a signal
!>  @param[inout] signals_created Number valid signals in the array, or index of
!>                                the last signal created.
!>  @param[in]    observed        Array of observed values.
!>  @param[in]    sigma           Array of observed sigma values
!>  @param[in]    weight          Array of signal weight values
!>  @param[inout] first_index     Index of the first instance of an intpol
!>                                signal.
!>  @todo FIXME: Remove the _new when removing the old code.
!-------------------------------------------------------------------------------
      SUBROUTINE ece_travis_dot_read(travis_file, signals,              &
     &                        signals_created,                          & 
     &                        observed, sigma, weight, first_index)
      USE signal
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=path_length), INTENT(in)            :: travis_file
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout) :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in)             :: observed
      REAL (rprec), DIMENSION(:), INTENT(in)             :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in)             :: weight
      INTEGER, INTENT(inout)                             :: first_index

!  local parameters
!  n_travis_keyword        integer - number of travis_keywords.
      INTEGER, PARAMETER :: n_travis_keyword = 2

!  local variables
!  travis_keyword    character array - keywords for various ece input types
      CHARACTER(len=data_name_length), DIMENSION(1:n_travis_keyword) ::           &
     &   travis_keyword
      TYPE(signal_dot_file)                 :: travis_dot_file_ref
      INTEGER                               :: current_profile
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Initialize the keywords
      travis_keyword(1) = 'ece_travis_channels'
      travis_keyword(2) = 'end_of_file'

!  Open up the 'ece_travis.' file
      travis_dot_file_ref = signal_dot_open(TRIM(travis_file),"travis")

!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(travis_dot_file_ref,         &
     &                                        travis_keyword))

            CASE DEFAULT ! This case should never fire.
               EXIT ! Exit out of infinte loop.

            CASE ('end_of_file')
               EXIT ! Exit out of infinte loop.

            CASE ('ece_travis_channels')
               CALL  ece_travis_dot_parse_chord(travis_dot_file_ref,     &
     &                                  signals, signals_created,        &
     &                                  observed, sigma, weight,         &
     &                                  first_index)

         END SELECT

      END DO

!  Close the 'ece.' file
      CALL signal_dot_close(travis_dot_file_ref)

      CALL profiler_set_stop_time('ece_travis_dot_read', start_time)

      END SUBROUTINE

!*******************************************************************************
!  PARSING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Parse ECE diagnostic chord.
!>
!>  Parses a single chord. Chords are laid out as,
!>
!>  * chord_name
!>  * Resonance as a single real.
!>  * Starting point as three reals.
!>  * Ending point as three reals.
!>
!>  @param[in]    ece_dot_file_ref Reference to the ece. file.
!>  @param[inout] signals          Array of pointers to the a signal
!>  @param[inout] signals_created  Number valid signals in the array, or index
!>                                 of the last signal created.
!>  @param[in]    observed         Array of observed values.
!>  @param[in]    sigma            Array of observed sigma values
!>  @param[in]    weight           Array of signal weight values
!>  @param[inout] first_index      Index of the first instance of an intpol
!>                                 signal.
!-------------------------------------------------------------------------------
      SUBROUTINE ece_travis_dot_parse_chord(travis_dot_file_ref,
     &                               signals, signals_created, observed,       &
     &                               sigma, weight, first_index)
      USE signal
      USE v3fit_input, only: v3fit_max_diagnostics
      USE v3_utilities, only: err_fatal

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_dot_file), INTENT(inout)  :: travis_dot_file_ref
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout)                 :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in) :: observed
      REAL (rprec), DIMENSION(:), INTENT(in) :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in) :: weight
      INTEGER, INTENT(inout)                 :: first_index

!  local variables
      CHARACTER (len=data_short_name_length) :: base_name
      CHARACTER (len=data_short_name_length) :: name
      CHARACTER (len=data_short_name_length) :: temp
      CHARACTER (len=path_length)            :: data_filename
      REAL (rprec)                           :: start_time
      INTEGER                                :: num_channels
      INTEGER                                :: index
      INTEGER, SAVE                          :: counter=0
      CHARACTER (len=8)                      :: fmt
!  Start of executable code
      start_time = profiler_get_start_time()

      base_name = signal_dot_parse_id(travis_dot_file_ref,                     &
     &               'Expected the id for TRAVIS ECE CHANNELS')

      data_filename = signal_dot_parse_filename(travis_dot_file_ref,           &
     &               'Expected the TRAVIS data filename for ECE')

      num_channels = signal_dot_parse_int(travis_dot_file_ref,                 &
     &               'Expected the number of channels for TRAVIS ECE')
      
!  Check if there are two many signals.
      IF (signals_created +num_channels .gt. v3fit_max_diagnostics) THEN
         CALL err_fatal('ece_travis_dot_parse_chord_new:' //                   &
     &                 'created signals exceeds v3fit_max_diagnostics')
      END IF

      fmt='(I3.3)'
      DO index = 1, num_channels
! Give each ece measurement a unique name
! Counter stores the last index of a ece diagnostic incase multiple 
! antennas are used with the same id
         counter = counter+1
         WRITE(temp,fmt) counter 
         name = trim(base_name)// trim(temp) 
         signals(signals_created + 1)%p => signal_construct(                    &
     &        diagnostic_construct(                                             &
     &        ece_travis_construct(data_filename,index)),                       &
     &        name, name, 'arb', observed(signals_created +1),                  &
     &        sigma(signals_created + 1), weight(signals_created + 1))
         signals_created = signals_created + 1

!  At lease one ece signal was made. Set the first index. This should only be
!  run once.
         IF (first_index .eq. -1) THEN
            first_index = signals_created
         END IF
      ENDDO   


      CALL profiler_set_stop_time('ece_dot_parse_chord', start_time)

      END SUBROUTINE

      END MODULE
