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
!>  about the specification of a ECE diagnostics that use TRAVIS to calculate the
!>  ECE signal.      
!>
!>  The ece_travis signals are specified in a structured text file. This file is read in by the
!>  @ref v3fit_input::ece_travis_dot_filename variable of the name list input file.
!>
!>  @section ece_travis_dot_key_sec Keywords
!>  @begin_table
!>  @item2{@fixed_width{ece_travis_channel}, A set of ECE channels specified by a TRAVIS ece_data file.}
!>  @end_table
!>
!>  @section ece_travis_dot_sig_spec_sec Channel Specification
!>  All TRAVIS ECE channels are specified in the following manner.
!>
!>  @fixed_width{key_word\n}
!>  @fixed_width{name(STRING len:@ref data_parameters::data_short_name_length)\n}
!>  @fixed_width{antenna index (INTEGER)\n}
!>  @fixed_width{mode type ( STRING len:1 (x/o))\n}
!>  @fixed_width{freq (GHz) (REAL)\n}
!>
!>
!>  @section ece_travis_dot_exam_sec Example File
!>  @code
!>  !  Sample TRAVIS ECE DOT FILE.
!>  ece_travis_channel
!>  ece001
!>  01
!>  x
!>  18.0
!>  ece_travis_channel
!>  ece002
!>  01
!>  o
!>  19.0
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
!>  @brief Contains module @ref ece_travis_dot
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
!>  @brief Read an ECE TRAVIS diagnostic dot file.
!>
!>  Parses the structure of a dignostic dot file and constructs a @ref ece
!>  object. Each chord is parsed by @ref ece_dot_parse_chord.
!>
!>  @param[in]    travis_file     Name of the ece_travis. file.
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
      travis_keyword(1) = 'ece_travis_channel'
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

            CASE ('ece_travis_channel')
               CALL  ece_travis_dot_parse_chord(travis_dot_file_ref,     &
     &                                  signals, signals_created,        &
     &                                  observed, sigma, weight,         &
     &                                  first_index)

         END SELECT

      END DO

!  Close the 'ece_travis.' file
      CALL signal_dot_close(travis_dot_file_ref)

      CALL profiler_set_stop_time('ece_travis_dot_read', start_time)

      END SUBROUTINE

!*******************************************************************************
!  PARSING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Parse ECE TRAVIS channel info.
!>
!>  Parses a single ECE Channel. Channels are laid out as,
!>
!>  * channel name
!>  * index of ece travis antenna as a integer
!>  * mode type as a character (x or o)
!>  * mode frequency in GHZ as a real
!>
!>  @param[in]    travis_dot_file_ref Reference to the ece_travis. file.
!>  @param[inout] signals             Array of pointers to the a signal
!>  @param[inout] signals_created     Number valid signals in the array, or index
!>                                    of the last signal created.
!>  @param[in]    observed            Array of observed values.
!>  @param[in]    sigma               Array of observed sigma values
!>  @param[in]    weight              Array of signal weight values
!>  @param[inout] first_index         Index of the first instance of an intpol
!>                                    signal.
!-------------------------------------------------------------------------------
      SUBROUTINE ece_travis_dot_parse_chord(travis_dot_file_ref,
     &                               signals, signals_created, observed,       &
     &                               sigma, weight, first_index)
      USE signal
      USE v3fit_input
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
      CHARACTER (len=data_short_name_length) :: name
      REAL (rprec)                           :: start_time
      INTEGER                                :: antenna_id
      INTEGER                                :: wmode
      CHARACTER (len=data_short_name_length) :: cmode
      REAL(rprec)                            :: freq

! Start of executable code
      start_time = profiler_get_start_time()


      
      name = signal_dot_parse_id(travis_dot_file_ref,                        &
     &               'Expected the id for TRAVIS ECE channel')

      antenna_id = signal_dot_parse_int(travis_dot_file_ref,                 &
     &               'Expected the TRAVIS ECE channel antenna')

      cmode =  signal_dot_parse_id(travis_dot_file_ref,                      &
     &     'Expected the ece mode type')


!> @TODO TRAVIS add check to reconize a space before x or o

      IF ((cmode(1:1)=='x') .OR. (cmode(1:1)=='X')) THEN
         wmode = 1
      ELSE IF ((cmode(1:1)=='o') .OR. (cmode(1:1)=='O')) THEN
         wmode = 0
      ELSE
         WRITE(*,*) "ERROR creating TRAVIS ECE CHANNEL"
         WRITE(*,*) "The mode type is not reconized "
         STOP
      END IF
      
      freq  =  signal_dot_parse_real(travis_dot_file_ref,                 &
     &           'Expected the TRAVIS channel frequeny')

! Check if there are two many signals.
      IF (signals_created +1 .gt. v3fit_max_diagnostics) THEN
         CALL err_fatal('ece_travis_dot_parse_chord_new:' //                   &
     &                 'created signals exceeds v3fit_max_diagnostics')
      END IF

         signals(signals_created + 1)%p => signal_construct(                    &
     &        diagnostic_construct(                                             &
     &     ece_travis_construct(antenna_id,wmode,freq) ),                       &
     &        name, name, 'arb', observed(signals_created +1),                  &
     &        sigma(signals_created + 1), weight(signals_created + 1))
         signals_created = signals_created + 1
!  Set the first index of an travis ece channel. This should only be run once.
         IF (first_index .eq. -1) THEN
            first_index = signals_created
#if defined(TRAV)
#else
         WRITE(*,*) "V3FIT must be linked to the TRAVIS library to use "       &
     &            //"the ece travis diagnostic! This build of V3FIT is "       &
     &           // " not linked to TRAVIS."
         STOP
#endif
         END IF


      CALL profiler_set_stop_time('ece_dot_parse_chord', start_time)

      END SUBROUTINE

      END MODULE
