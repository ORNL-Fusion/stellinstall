!-------------------------------------------------------------------------------
!  The @fixed_width, @begin_table, @item2 and @end_table commands are custom
!  defined commands in Doxygen.in. They are defined under ALIASES. For the page
!  created here, the 80 column limit is exceeded. Arguments of aliases are
!  separated by ','. You intended ',' to be a string you must use an escaped
!  comma '\,'.
!
!>  @page intpol_dot_sec Inteferometry/Polarimetry Diagnostic Dot File
!>
!>  @tableofcontents
!>  @section intpol_dot_intro_sec Introduction
!>  Diagnostic dot files specify geometric and auxiliary information about each
!>  diagnostic. Parsing of a diagnostic dot file begins by searching for a valid
!>  keyword. All lines before a valid keyword is reached are ignored allowing
!>  for comments and other information to be inserted in the file. Once a
!>  keyword is reached, the parsing becomes strict. The diagnostic description
!>  format is dependant on the specific keyword. All diagnostic dot files must
!>  end with an @fixed_width{end_of_file} keyword.
!>
!>  This document contains information about the interferometry/polarimetry
!>  chord specifications. Chords are specified in a structured text file. This
!>  file is read in by the @ref v3fit_input::ipch_dot_filename variable of the
!>  name list input file. Chords are specified as two verticies in either meters
!>  or radians. For polarimetry, wavelengths are specified in meters.
!>
!>  @section intpol_dot_key_sec Keywords
!>  @begin_table
!>  @item2{@fixed_width{ip_chord_XYZ}\n @fixed_width{ip_chord_XYZ_int},           An interferometer chord specified in Cartesian coordinates.}
!>  @item2{@fixed_width{ip_chord_RPhiDegZ}\n @fixed_width{ip_chord_RPhiDegZ_int}, An interferometer chord specified in Cylindrical coordinates.}
!>  @item2{@fixed_width{ip_chord_XYZ_pol_rad},                                    An polarimetry chord specified in Cartesian coordinates outputting the signal in radians.}
!>  @item2{@fixed_width{ip_chord_RPhiDegZ_pol_rad},                               An polarimetry chord specified in Cylindrical coordinates outputting the signal in radians.}
!>  @item2{@fixed_width{ip_chord_XYZ_pol_deg},                                    An polarimetry chord specified in Cartesian coordinates outputting the signal in degress.}
!>  @item2{@fixed_width{ip_chord_RPhiDegZ_pol_deg},                               An polarimetry chord specified in Cylindrical coordinates outputting the signal in degress.}
!>  @end_table
!>
!>  @section intpol_dot_chord_spec_sec Chord Specification
!>  All interferometry chords are specified in the following manner.
!>
!>  @fixed_width{key_word\n}
!>  @fixed_width{name(STRING len:@ref data_parameters::data_short_name_length)\n}
!>  @fixed_width{x1(REAL) y1(REAL) z1(REAL)\n}
!>  @fixed_width{x2(REAL) y2(REAL) z2(REAL)\n}
!>
!>  All Polarimetry chords are specified in the following manner with the wave
!>  length .
!>
!>  @fixed_width{key_word\n}
!>  @fixed_width{name(STRING len:@ref data_parameters::data_short_name_length)\n}
!>  @fixed_width{x1(REAL) y1(REAL) z1(REAL)\n}
!>  @fixed_width{x2(REAL) y2(REAL) z2(REAL)\n}
!>  @fixed_width{wavelength(REAL)\n}
!>
!>  @section intpol_dot_exam_sec Example File
!>  @code
!>  !  An interferometer chord specified in Cartesian coordinates.
!>  ip_chord_XYZ_int
!>  int001
!>  1.0 0.0 0.0
!>  0.5 0.0 0.0
!>
!>  !  An interferometer chord specified in Cylindrical coordinates.
!>  ip_chord_RPhiDegZ_int
!>  int002
!>  1.0 0.0 0.0
!>  0.5 0.0 0.0
!>
!>  !  A polarimetry chord specified in Cartesian coordinates.
!>  ip_chord_XYZ_pol_deg
!>  pol001
!>  1.0 0.0 0.0
!>  0.5 0.0 0.0
!>  1.0E-5
!>
!>  !  A polarimetry chord specified in Cylindrical coordinates.
!>  ip_chord_RPhiDegZ_pol_deg
!>  int002
!>  1.0 0.0 0.0
!>  0.5 0.0 0.0
!>  1.0E-5
!>
!>  !  All files must end with
!>  end_of_file
!>  @endcode
!>
!>  @section intpol_dot_prog_ref_sec Programmers Reference
!>  Reference material for the coding to parse these files can be found in the
!>  @ref signal_dot and @ref ipch_dot modules.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file ipch_dot.f
!>  @brief Contains module @ref ipch_dot
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module for opening and reading a 'ipch.' file. The file format for these
!>  files are documented in @ref intpol_dot_sec.
!*******************************************************************************

      MODULE ipch_dot

      USE stel_kinds
      USE stel_constants
      USE ipch_T
      USE signal_dot

      IMPLICIT NONE
          
!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for reading the interferometry-polarimetry diagnostic dot file
!>  using @ref ipch_dot_read_old or @ref ipch_dot_read_new.
!>  @todo. FIXME:  Remove inteface when old coding is removed.
!-------------------------------------------------------------------------------
      INTERFACE ipch_dot_read
         MODULE PROCEDURE ipch_dot_read_old,                                   &
     &                    ipch_dot_read_new
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for parsing a interferometry-polarimetry chord using
!>  @ref ipch_dot_read_old or @ref ipch_dot_read_new.
!>  @todo. FIXME:  Remove inteface when old ocding is removed.
!-------------------------------------------------------------------------------
      INTERFACE ipch_dot_parse_chord
         MODULE PROCEDURE ipch_dot_parse_chord_old,                            &
     &                    ipch_dot_parse_chord_new
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!  Subroutine to parse the large structure of the ipch_dot_file
!  Individual ip chords are parsed by the ipch_dot_parse_xxxx subroutines
!>  @todo Remove when removing old coding.
!-------------------------------------------------------------------------------
      SUBROUTINE ipch_dot_read_old (ipch_file,ipch_desc_arr,n_chords)

      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in)            :: ipch_file
      TYPE (ipch_desc), DIMENSION(:), INTENT(inout) :: ipch_desc_arr
      INTEGER, INTENT(inout)  :: n_chords


!  ipch_file            Name of the ipch. file.
!  ipch_desc_arr        Array of ipch_desc. Assumed to have been allocated
!                        before call to ipch_dot_read
!  n_chords              Number of chords read from the file

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n_ipchd_keyword = 9
!  n_ipchd_keyword        integer - number of ipchd_keyword s.

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=data_name_length), DIMENSION(1:n_ipchd_keyword) ::         &
     &   ipchd_keyword

      INTEGER               :: icount_chords
      INTEGER               :: n_ipch_desc_arr

      TYPE(signal_dot_file) :: ipch_dot_file_ref
!  Make sure that len for this_keyword is the same as len for ipchd_keyword.

!  ipchd_keyword    character array - keywords for various ipch input types

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Initialize the keywords and keyword lengths
      ipchd_keyword(1) = 'ip_chord_XYZ'
      ipchd_keyword(2) = 'ip_chord_RPhiDegZ'
      ipchd_keyword(3) = 'ip_chord_XYZ_int'
      ipchd_keyword(4) = 'ip_chord_RPhiDegZ_int'
      ipchd_keyword(5) = 'ip_chord_XYZ_pol_rad'
      ipchd_keyword(6) = 'ip_chord_RPhiDegZ_pol_rad'
      ipchd_keyword(7) = 'ip_chord_XYZ_pol_deg'
      ipchd_keyword(8) = 'ip_chord_RPhiDegZ_pol_deg'
      ipchd_keyword(9) = 'end_of_file'

!  Find the length of the ipch_desc_arr array
      n_ipch_desc_arr = SIZE(ipch_desc_arr)

!  Open up the 'ipch.' file
      ipch_dot_file_ref = signal_dot_open(TRIM(ipch_file), "ipch")

      icount_chords = 0
!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(ipch_dot_file_ref,               &
     &                                        ipchd_keyword))
         
            CASE DEFAULT ! This case should never fire.
               EXIT ! Exit out of infinte loop.

            CASE ('end_of_file')
               EXIT ! Exit out of infinte loop.
            
            CASE ('ip_chord_XYZ','ip_chord_XYZ_int')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'XYZ', 'i', .false.,                  &
     &                                   icount_chords,                        &
     &                                   n_ipch_desc_arr,                      &
     &                                   ipch_desc_arr)
               
            CASE ('ip_chord_RPhiDegZ','ip_chord_RPhiDegZ_int')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'RPHiDegZ', 'i', .false.,             &
     &                                   icount_chords,                        &
     &                                   n_ipch_desc_arr,                      &
     &                                   ipch_desc_arr)

            CASE ('ip_chord_XYZ_pol_rad')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'XYZ', 'p', .false.,                  &
     &                                   icount_chords,                        &
     &                                   n_ipch_desc_arr,                      &
     &                                   ipch_desc_arr)

            CASE ('ip_chord_RPhiDegZ_pol_rad')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'RPHiDegZ', 'p', .false.,             &
     &                                   icount_chords,                        &
     &                                   n_ipch_desc_arr,                      &
     &                                   ipch_desc_arr)

            CASE ('ip_chord_XYZ_pol_deg')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'XYZ', 'p', .true.,                   &
     &                                   icount_chords,                        &
     &                                   n_ipch_desc_arr,                      &
     &                                   ipch_desc_arr)

            CASE ('ip_chord_RPhiDegZ_pol_deg')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'RPHiDegZ', 'p', .true.,              &
     &                                   icount_chords,                        &
     &                                   n_ipch_desc_arr,                      &
     &                                   ipch_desc_arr)

         END SELECT

!  Check for too many additions to the ipch_desc_arr array
         IF  (icount_chords .gt. SIZE(ipch_desc_arr)) THEN
            STOP 'ipch_dot: Too many coils for ipch_desc_arr '
         ENDIF
         
      END DO

!  Close the 'ipch.' file
      CALL signal_dot_close(ipch_dot_file_ref)
      
      n_chords = icount_chords
      
      RETURN
      END SUBROUTINE ipch_dot_read_old

!-------------------------------------------------------------------------------
!>  @brief Read an interferometry/polarimetry diagnostic dot file.
!>
!>  Parses the structure of a dignostic dot file and constructs a @ref intpol
!>  object. Each chord is parsed by @ref ipch_dot_parse_chord.
!>
!>  @param[in]    ipch_file       Name of the ipch. file.
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
      SUBROUTINE ipch_dot_read_new(ipch_file, signals,                         &
     &                             signals_created, observed, sigma,           &
     &                             weight, first_index)
      USE signal
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_filename_length), INTENT(in)   :: ipch_file
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout) :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in)             :: observed
      REAL (rprec), DIMENSION(:), INTENT(in)             :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in)             :: weight
      INTEGER, INTENT(inout)                             :: first_index

!  local parameters
!  n_ipchd_keyword        integer - number of ipchd_keyword s.
      INTEGER, PARAMETER :: n_ipchd_keyword = 9

!  local variables
!  ipchd_keyword    character array - keywords for various ipch input types
      CHARACTER(len=data_name_length), DIMENSION(1:n_ipchd_keyword) ::         &
     &   ipchd_keyword
      TYPE(signal_dot_file) :: ipch_dot_file_ref

!  Start of executable code
!  Initialize the keywords and keyword lengths
      ipchd_keyword(1) = 'ip_chord_XYZ'
      ipchd_keyword(2) = 'ip_chord_RPhiDegZ'
      ipchd_keyword(3) = 'ip_chord_XYZ_int'
      ipchd_keyword(4) = 'ip_chord_RPhiDegZ_int'
      ipchd_keyword(5) = 'ip_chord_XYZ_pol_rad'
      ipchd_keyword(6) = 'ip_chord_RPhiDegZ_pol_rad'
      ipchd_keyword(7) = 'ip_chord_XYZ_pol_deg'
      ipchd_keyword(8) = 'ip_chord_RPhiDegZ_pol_deg'
      ipchd_keyword(9) = 'end_of_file'

!  Open up the 'ipch.' file
      ipch_dot_file_ref = signal_dot_open(TRIM(ipch_file), "ipch")

!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(ipch_dot_file_ref,               &
     &                                        ipchd_keyword))

            CASE DEFAULT ! This case should never fire.
               EXIT ! Exit out of infinte loop.

            CASE ('end_of_file')
               EXIT ! Exit out of infinte loop.

            CASE ('ip_chord_XYZ','ip_chord_XYZ_int')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'XYZ', 'i', .false.,                  &
     &                                    signals, signals_created,            &
     &                                    observed, sigma, weight,             &
     &                                    first_index)


            CASE ('ip_chord_RPhiDegZ','ip_chord_RPhiDegZ_int')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'RPHiDegZ', 'i', .false.,             &
     &                                    signals, signals_created,            &
     &                                    observed, sigma, weight,             &
     &                                    first_index)


            CASE ('ip_chord_XYZ_pol_rad')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'XYZ', 'p', .false.,                  &
     &                                    signals, signals_created,            &
     &                                    observed, sigma, weight,             &
     &                                    first_index)


            CASE ('ip_chord_RPhiDegZ_pol_rad')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'RPHiDegZ', 'p', .false.,             &
     &                                    signals, signals_created,            &
     &                                    observed, sigma, weight,             &
     &                                    first_index)


            CASE ('ip_chord_XYZ_pol_deg')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'XYZ', 'p', .true.,                   &
     &                                    signals, signals_created,            &
     &                                    observed, sigma, weight,             &
     &                                    first_index)


            CASE ('ip_chord_RPhiDegZ_pol_deg')
               CALL ipch_dot_parse_chord(ipch_dot_file_ref,                    &
     &                                   'RPHiDegZ', 'p', .true.,              &
     &                                    signals, signals_created,            &
     &                                    observed, sigma, weight,             &
     &                                    first_index)


         END SELECT

      END DO

!  Close the 'ipch.' file
      CALL signal_dot_close(ipch_dot_file_ref)

      END SUBROUTINE ipch_dot_read_new

!*******************************************************************************
!  PARSING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!     Expected behaviour for the parsing routines:
!  On entry:
!     - character variable 'line' is defined
!     - logical variable ipchd_l_done is .false.
!     - integer variable i_line is set
!  During Execution:
!     - integer variable i_line is kept current
!     - reads lines from file, looks for appropriate line structure
!     - presence of inappropriate line indicates end of chord
!     - constructs ipch_desc_temp
!  On exit:
!     - 'line' contains the first inappropriate line
!     - End of File or other problem: ipchd_l_done set to .true.
!     
!-------------------------------------------------------------------------------
!  Subroutine to parse an ipch chord (redundant chord there)
!      ipchd_keyword(1) = 'ip_chord_XYZ' or 'ip_chord_RPhiDegZ'
!>  @todo Remove when removing old coding.
!-------------------------------------------------------------------------------
      SUBROUTINE ipch_dot_parse_chord_old(ipch_dot_file_ref,                   &
     &                                    coordinate_type,                     &
     &                                    chord_type,                          &
     &                                    inDegrees,                           &
     &                                    icount_chords,                       &
     &                                    n_ipch_desc_arr,                     &
     &                                    ipch_desc_arr)
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(inout)                       :: icount_chords
      INTEGER, INTENT(in)                          :: n_ipch_desc_arr
      TYPE(ipch_desc), DIMENSION(:), INTENT(inout) :: ipch_desc_arr
      CHARACTER (len=*) :: coordinate_type
      CHARACTER (len=1) :: chord_type
      LOGICAL, INTENT(in) :: inDegrees
      TYPE(signal_dot_file), INTENT(inout) :: ipch_dot_file_ref

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(3)        :: xcart_i, xcart_f
      CHARACTER (len=data_name_length) :: chord_name
      REAL(rprec)                      :: wavelength = 0.0
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Get the start, end and name of chord.
      CALL signal_dot_parse_chord(ipch_dot_file_ref,                           &
     &                            coordinate_type,                             &
     &                            chord_name,                                  &
     &                            xcart_i, xcart_f)

      IF (chord_type .eq. 'p') THEN
         wavelength = signal_dot_parse_real(ipch_dot_file_ref,                 &
     &                   "Expected wavelength for polarimetry chord")
      END IF

      icount_chords = icount_chords + 1
      IF (icount_chords .gt. n_ipch_desc_arr) THEN
         WRITE(*,*) 'TOO MANY ipCH. Increase dimensions'
      ELSE
         CALL ipch_desc_construct(ipch_desc_arr(icount_chords),                &
     &                            chord_name, xcart_i, xcart_f,                &
     &                            chord_type, wavelength, inDegrees)
      ENDIF

      END SUBROUTINE ipch_dot_parse_chord_old

!-------------------------------------------------------------------------------
!>  @brief Parse interferometry/polarimetry diagnostic chord.
!>
!>  Parses a signle chord. Chords are laid out as,
!>
!>  * chord_name
!>  * Starting point as three reals.
!>  * Ending point as three reals.
!>
!>  If the signal is a polarimetry chord, then also read a single real for the
!>  wavelength.
!>
!>  @param[in]    ipch_file       Name of the ipch. file.
!>  @param[in]    coordinate_type Specifies the cartesian or cylindrical
!>                                coordinate type.
!>  @param[in]    chord_type      Specifies if the cord is a 'i'nterferometer or
!>                                'p'olarimetry.
!>  @param[in]    inDegrees       Specifies if signal measures degrees or
!>                                radians.
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
      SUBROUTINE ipch_dot_parse_chord_new(ipch_dot_file_ref,                   &
     &                                    coordinate_type,                     &
     &                                    chord_type,                          &
     &                                    inDegrees,                           &
     &                                    signals,                             &
     &                                    signals_created,                     &
     &                                    observed, sigma, weight,             &
     &                                    first_index) 
      USE signal
      USE v3fit_input, only: v3fit_max_diagnostics
      USE v3_utilities, only: err_fatal
      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_dot_file), INTENT(inout)  :: ipch_dot_file_ref
      CHARACTER (len=*), INTENT(in)          :: coordinate_type
      CHARACTER (len=1), INTENT(in)          :: chord_type
      LOGICAL, INTENT(in)                    :: inDegrees
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout)                 :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in) :: observed
      REAL (rprec), DIMENSION(:), INTENT(in) :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in) :: weight
      INTEGER, INTENT(inout)                 :: first_index

!  local variables
      REAL (rprec), DIMENSION(2,3)     :: xcart
      CHARACTER (len=data_name_length) :: chord_name
      REAL (rprec)                     :: wavelength = 0.0

!  Start of executable code
!  Check if there are two many signals.
      IF (signals_created + 1 .gt. v3fit_max_diagnostics) THEN
         CALL err_fatal('ipch_dot_parse_chord_new: created signals' //         &
     &                  ' exceeds v3fit_max_diagnostics')
      END IF

!  Get the start, end and name of chord.
      CALL signal_dot_parse_chord(ipch_dot_file_ref,                           &
     &                            coordinate_type,                             &
     &                            chord_name,                                  &
     &                            xcart(1,:), xcart(2,:))

      IF (chord_type .eq. 'p') THEN
         wavelength = signal_dot_parse_real(ipch_dot_file_ref,                 &
     &                   "Expected wavelength for polarimetry chord")
         signals(signals_created + 1)%p => signal_construct(                   &
     &      diagnostic_construct(                                              &
     &         intpol_construct(wavelength, inDegrees, xcart)),                &
     &      chord_name, chord_name, 'm^-2',                                    &
     &      observed(signals_created + 1), sigma(signals_created + 1),         &
     &      weight(signals_created + 1))
      ELSE
         signals(signals_created + 1)%p => signal_construct(                   &
     &      diagnostic_construct(                                              &
     &         intpol_construct(xcart)), chord_name, chord_name,               &
     &      'm^-2', observed(signals_created + 1),                             &
     &      sigma(signals_created + 1), weight(signals_created + 1))
      END IF

      signals_created = signals_created + 1

!  At lease one intpol signal was made. Set the first index. This should only be
!  run once.
      IF (first_index .eq. -1) THEN
         first_index = signals_created
      END IF

      END SUBROUTINE ipch_dot_parse_chord_new

      END MODULE ipch_dot
