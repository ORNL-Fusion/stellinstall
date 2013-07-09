!-------------------------------------------------------------------------------
!  The @fixed_width, @begin_table, @item2 and @end_table commands are custom
!  defined commands in Doxygen.in. They are defined under ALIASES. For the page
!  created here, the 8 column limit is exceeded. Arguments of aliases are
!  separated by ','. You intended ',' to be a string you must use an escaped
!  comma '\,'.
!
!>  @page sxrem_dot_sec Soft X-ray Diagnostic Dot File
!>
!>  @tableofcontents
!>  @section sxrem_dot_intro_sec Introduction
!>  Diagnostic dot files specify geometric and auxiliary information about each
!>  diagnostic. Parsing of a diagnostic dot file begins by searching for a valid
!>  keyword. All lines before a valid keyword is reached are ignored allowing
!>  for comments and other information to be inserted in the file. Once a
!>  keyword is reached, the parsing becomes strict. The diagnostic description
!>  format is dependant on the specific keyword. All diagnostic dot files must
!>  end with an @fixed_width{end_of_file} keyword.
!>
!>  This document contains information about the soft x-ray chord
!>  specifications. Chords are specified in a structured text file. This file is
!>  read in by the @ref v3fit_input::sxrch_dot_filename variable of the name
!>  list input file. Chords are specified as two verticies in either meters or
!>  radians.
!>
!>  @section sxrem_dot_key_sec Keywords
!>  @begin_table
!>  @item2{@fixed_width{sxr_chord_chord_XYZ},
!>  A soft x-ray chord specified in Cartesian coordinates.}
!>  @item2{@fixed_width{sxr_chord_chord_RPhiDegZ},
!>  A soft x-ray chord specified in Cylindrical coordinates.}
!>  @end_table
!>
!>  @section sxrem_dot_chord_spec_sec Chord Specification
!>  All soft x-ray chords are specified in the following manner.
!>
!>  @fixed_width{key_word\n}
!>  @fixed_width{name(STRING len:@ref data_parameters::data_short_name_length)\n}
!>  @fixed_width{x1(REAL) y1(REAL) z1(REAL)\n}
!>  @fixed_width{x2(REAL) y2(REAL) z2(REAL)\n}
!>
!>  @section sxrem_dot_exam_sec Example File
!>  @code
!>  !  A soft x-ray chord specified in Cartesian coordinates.
!>  sxr_chord_XYZ
!>  sxr001
!>  1.0 0.0 0.0
!>  0.5 0.0 0.0
!>
!>  !  A soft x-ray chord specified in Cylindrical coordinates.
!>  sxr_chord_RPhiDegZ
!>  sxr002
!>  1.0 0.0 0.0
!>  0.5 0.0 0.0
!>  
!>  !  All files must end with
!>  end_of_file
!>  @endcode
!>
!>  @section sxrem_dot_prog_ref_sec Programmers Reference
!>  Reference material for the coding to parse these files can be found in the
!>  @ref signal_dot and @ref sxrch_dot modules.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file sxrch_dot.f
!>  @brief Contains module @ref sxrch_dot
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module for opening and reading a 'sxrch.' file. The file format for these
!>  files are documented in @ref sxrem_dot_sec.
!*******************************************************************************

      MODULE sxrch_dot


      USE stel_kinds
      USE stel_constants
      USE sxrch_T
      USE signal_dot

      IMPLICIT NONE
          
!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for reading the soft x-ray diagnostic dot file using
!>  @ref sxrch_dot_read_old or @ref sxrch_dot_read_new.
!>  @todo. FIXME:  Remove inteface when old coding is removed.
!-------------------------------------------------------------------------------
      INTERFACE sxrch_dot_read
         MODULE PROCEDURE sxrch_dot_read_old,                                  &
     &                    sxrch_dot_read_new
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for parsing a soft x-ray chord using
!>  @ref sxrch_dot_parse_chord_old or @ref sxrch_dot_parse_chord_new.
!>  @todo. FIXME:  Remove inteface when old coding is removed.
!-------------------------------------------------------------------------------
      INTERFACE sxrch_dot_parse_chord
         MODULE PROCEDURE sxrch_dot_parse_chord_old,                           &
     &                    sxrch_dot_parse_chord_new
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!  Subroutine to parse the large structure of the sxrch_dot_file
!  Individual sxr chords are parsed by the sxrch_dot_parse_xxxx subroutines
!>  @todo Remove when removing old coding.
!-------------------------------------------------------------------------------
      SUBROUTINE sxrch_dot_read_old (sxrch_file,sxrch_desc_arr,n_chords)

      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in)            :: sxrch_file
      TYPE (sxrch_desc), DIMENSION(:), INTENT(inout) :: sxrch_desc_arr
      INTEGER, INTENT(inout)  :: n_chords


!  sxrch_file            Name of the sxrch. file.
!  sxrch_desc_arr        Array of sxrch_desc. Assumed to have been allocated
!                        before call to sxrch_dot_read
!  n_chords              Number of chords read from the file

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n_sxrchd_keyword = 3
!  n_sxrchd_keyword        integer - number of sxrchd_keyword s.

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=data_name_length), DIMENSION(1:n_sxrchd_keyword) ::        &
     &   sxrchd_keyword

      INTEGER               :: icount_chords
      INTEGER               :: n_sxrch_desc_arr

      TYPE(signal_dot_file) :: sxrch_dot_file_ref
!  Make sure that len for this_keyword is the same as len for sxrchd_keyword.

!  sxrchd_keyword    character array - keywords for various sxrch input types
!  n_char_key        integer array - lengths of the sxrchd_keyword strings

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Initialize the keywords and keyword lengths
      sxrchd_keyword(1) = 'sxr_chord_XYZ'
      sxrchd_keyword(2) = 'sxr_chord_RPhiDegZ'
      sxrchd_keyword(3) = 'end_of_file'

!  Find the length of the sxrch_desc_arr array
      n_sxrch_desc_arr = SIZE(sxrch_desc_arr)

!  Open up the 'sxrch.' file
      sxrch_dot_file_ref = signal_dot_open(TRIM(sxrch_file),"sxrch")

      icount_chords = 0
!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(sxrch_dot_file_ref,              &
     &                                        sxrchd_keyword))
         
            CASE DEFAULT ! This case should never fire.
               EXIT ! Exit out of infinte loop.

            CASE ('end_of_file')
               EXIT ! Exit out of infinte loop.
            
            CASE ('sxr_chord_XYZ')
               CALL sxrch_dot_parse_chord(sxrch_dot_file_ref,                  &
     &                                    'XYZ', icount_chords,                &
     &                                    n_sxrch_desc_arr,                    &
     &                                    sxrch_desc_arr)
               
            CASE ('sxr_chord_RPhiDegZ')
               CALL sxrch_dot_parse_chord(sxrch_dot_file_ref,                  &
     &                                    'RPHiDegZ', icount_chords,           &
     &                                    n_sxrch_desc_arr,                    &
     &                                    sxrch_desc_arr)

         END SELECT
         
!  Check for too many additions to the sxrch_desc_arr array
         IF  (icount_chords .gt. SIZE(sxrch_desc_arr)) THEN
            STOP 'sxrch_dot: Too many coils for sxrch_desc_arr '
         ENDIF
         
      END DO

!  Close the 'sxrch.' file
      CALL signal_dot_close(sxrch_dot_file_ref)
      
      n_chords = icount_chords
      
      RETURN
      END SUBROUTINE sxrch_dot_read_old

!-------------------------------------------------------------------------------
!>  @brief Read an soft x-ray diagnostic dot file.
!>
!>  Parses the structure of a dignostic dot file and constructs a @ref sxrem
!>  object. Each chord is parsed by @ref sxrch_dot_parse_chord.
!>
!>  @param[in]    sxrch_file      Name of the sxrch. file.
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
      SUBROUTINE sxrch_dot_read_new(sxrch_file, signals,                       &
     &                              signals_created, observed, sigma,          &
     &                              weight, first_index)
      USE signal
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_filename_length), INTENT(in)   :: sxrch_file
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout) :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in)             :: observed
      REAL (rprec), DIMENSION(:), INTENT(in)             :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in)             :: weight
      INTEGER, INTENT(inout)                             :: first_index

!  local parameters
!  n_sxrchd_keyword        integer - number of sxrchd_keyword s.
      INTEGER, PARAMETER :: n_sxrchd_keyword = 3

!  local variables
!  sxrchd_keyword    character array - keywords for various sxrch input types
      CHARACTER(len=data_name_length), DIMENSION(1:n_sxrchd_keyword) ::        &
     &   sxrchd_keyword
      TYPE(signal_dot_file)                 :: sxrch_dot_file_ref

!  Start of executable code
!  Initialize the keywords
      sxrchd_keyword(1) = 'sxr_chord_XYZ'
      sxrchd_keyword(2) = 'sxr_chord_RPhiDegZ'
      sxrchd_keyword(3) = 'end_of_file'

!  Open up the 'sxrch.' file
      sxrch_dot_file_ref = signal_dot_open(TRIM(sxrch_file),"sxrch")

!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(sxrch_dot_file_ref,              &
     &                                        sxrchd_keyword))

            CASE DEFAULT ! This case should never fire.
               EXIT ! Exit out of infinte loop.

            CASE ('end_of_file')
               EXIT ! Exit out of infinte loop.

            CASE ('sxr_chord_XYZ')
               CALL sxrch_dot_parse_chord(sxrch_dot_file_ref,                  &
     &                                    'XYZ', signals,                      &
     &                                    signals_created,                     &
     &                                    observed, sigma, weight,             &
     &                                    first_index)

            CASE ('sxr_chord_RPhiDegZ')
               CALL sxrch_dot_parse_chord(sxrch_dot_file_ref,                  &
     &                                    'RPHiDegZ', signals,                 &
     &                                    signals_created,                     &
     &                                    observed, sigma, weight,             &
     &                                    first_index)


         END SELECT

      END DO

!  Close the 'sxrch.' file
      CALL signal_dot_close(sxrch_dot_file_ref)

      END SUBROUTINE sxrch_dot_read_new

!*******************************************************************************
!  PARSING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!     Expected behaviour for the parsing routines:
!  On entry:
!     - character variable 'line' is defined
!     - logical variable sxrchd_l_done is .false.
!     - integer variable i_line is set
!  During Execution:
!     - integer variable i_line is kept current
!     - reads lines from file, looks for appropriate line structure
!     - presence of inappropriate line indicates end of chord
!     - constructs sxrch_desc_temp
!  On exit:
!     - 'line' contains the first inappropriate line
!     - End of File or other problem: sxrchd_l_done set to .true.
!     
!-------------------------------------------------------------------------------
!  Subroutine to parse an sxrch chord (redundant chord there)
!      sxrchd_keyword(1) = 'sxr_chord_XYZ' or 'sxr_chord_RPhiDegZ'
!>  @todo Remove when removing old coding.
!-------------------------------------------------------------------------------
      SUBROUTINE sxrch_dot_parse_chord_old(sxrch_dot_file_ref,                 &
     &                                     coordinate_type,                    &
     &                                     icount_chords,                      &
     &                                     n_sxrch_desc_arr,                   &
     &                                     sxrch_desc_arr)
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(inout)                        :: icount_chords
      INTEGER, INTENT(in)                           :: n_sxrch_desc_arr
      TYPE(sxrch_desc), DIMENSION(:), INTENT(inout) :: sxrch_desc_arr
      CHARACTER (len=*) :: coordinate_type
      TYPE(signal_dot_file), INTENT(inout) :: sxrch_dot_file_ref

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(3)      :: xcart_i, xcart_f
      CHARACTER (len=chord_name_len) :: chord_name
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Get the start, end and name of chord.
      CALL signal_dot_parse_chord(sxrch_dot_file_ref,                          &
     &                            coordinate_type,                             &
     &                            chord_name,                                  &
     &                            xcart_i, xcart_f)

      icount_chords = icount_chords + 1
      IF (icount_chords .gt. n_sxrch_desc_arr) THEN
         WRITE(*,*) 'TOO MANY SXRCH. Increase dimnesions'
      ELSE
         CALL sxrch_desc_construct(sxrch_desc_arr(icount_chords),              &
     &                             chord_name, xcart_i, xcart_f)
      ENDIF

      END SUBROUTINE sxrch_dot_parse_chord_old

!-------------------------------------------------------------------------------
!>  @brief Parse soft x-ray diagnostic chord.
!>
!>  Parses a signle chord. Chords are laid out as,
!>
!>  * chord_name
!>  * Starting point as three reals.
!>  * Ending point as three reals.
!>
!>  @param[in]    sxrch_file      Name of the sxrch. file.
!>  @param[in]    coordinate_type Specifies the cartesian or cylindrical
!>                                coordinate type.
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
      SUBROUTINE sxrch_dot_parse_chord_new(sxrch_dot_file_ref,                 &
     &                                     coordinate_type,                    &
     &                                     signals,                            &
     &                                     signals_created,                    &
     &                                     observed, sigma, weight,            &
     &                                     first_index)
      USE signal
      USE v3fit_input, only: v3fit_max_diagnostics
      USE v3_utilities, only: err_fatal

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_dot_file), INTENT(inout)  :: sxrch_dot_file_ref
      CHARACTER (len=*), INTENT(in)          :: coordinate_type
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout)                 :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in) :: observed
      REAL (rprec), DIMENSION(:), INTENT(in) :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in) :: weight
      INTEGER, INTENT(inout)                 :: first_index

!  local variables
      REAL(rprec), DIMENSION(3)      :: xcart_i, xcart_f
      CHARACTER (len=chord_name_len) :: chord_name

!  Start of executable code
!  Check if there are two many signals.
      IF (signals_created + 1 .gt. v3fit_max_diagnostics) THEN
         CALL err_fatal('sxrch_dot_parse_chord_new: created signals' //        &
     &                  ' exceeds v3fit_max_diagnostics')
      END IF

!  Get the start, end and name of chord.
      CALL signal_dot_parse_chord(sxrch_dot_file_ref,                          &
     &                            coordinate_type,                             &
     &                            chord_name,                                  &
     &                            xcart_i, xcart_f)

      signals(signals_created + 1)%p => signal_construct(                      &
     &   diagnostic_construct(                                                 &
     &      sxrem_construct(xcart_i, xcart_f)), chord_name, chord_name,        &
     &   'arb', observed(signals_created + 1),                                 &
     &   sigma(signals_created + 1), weight(signals_created + 1))
      signals_created = signals_created + 1

!  At lease one sxrem signal was made. Set the first index. This should only be
!  run once.
      IF (first_index .eq. -1) THEN
         first_index = signals_created
      END IF

      END SUBROUTINE sxrch_dot_parse_chord_new

      END MODULE sxrch_dot
