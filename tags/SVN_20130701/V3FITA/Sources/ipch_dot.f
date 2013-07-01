!*******************************************************************************
!  File ipch_dot.f
!  Contains module ipch_dot
!    Module for opening and reading a 'ipch.' file, and then placing the
!     interferometry-polarimetry chord information into an array of ipch_desc
!!  
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       stel_constants
!       ipch_T
!       safe_open_mod
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  Based on module sxrch_dot
!  First version JDH 2012-03-15
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!     Usage: Here is a sample of code that shows how to use this module
!
!      USE ipch_T
!      USE ipch_dot
!      TYPE(ipch_desc), DIMENSION(1000) :: ipch_descriptions_a
! ...
! !  Read in interferometry-polarimetry chords
!      CALL ipch_dot_read(file_name_ipch_dot,ipch_descriptions_a)
!
!-------------------------------------------------------------------------------
!   FILE STRUCTURE - ipch. file
!-------------------------------------------------------------------------------
!
!  ipch. File Structure:
! line 1: Character - identifier
!
!      (0 or more lines, no keyword present. All ignored by parser.)
!    keyword
!      lines containing data - number depends on keyword
!
!      (0 or more lines, no keyword present. All ignored by parser.)
!    keyword
!      lines containing data - number depends on keyword
!
!    ...
!
!    end
!      (0 or more lines, ignored by parser.)!    
!    
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
! 1)  Defined keywords are:
!       ipchd_keyword(1) = 'ip_chord_XYZ'
!       ipchd_keyword(2) = 'ip_chord_RPhiDegZ'
!       ipchd_keyword(3) = 'ip_chord_XYZ_int'
!       ipchd_keyword(4) = 'ip_chord_RPhiDegZ_int'
!       ipchd_keyword(5) = 'ip_chord_XYZ_pol_rad'
!       ipchd_keyword(6) = 'ip_chord_RPhiDegZ_pol_rad'
!       ipchd_keyword(7) = 'ip_chord_XYZ_pol_deg'
!       ipchd_keyword(8) = 'ip_chord_RPhiDegZ_pol_deg'
!       ipchd_keyword(9) = 'end_of_file' 
!
!*******************************************************************************

!*******************************************************************************
!  MODULE ipch_dot
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
! SECTION IV. PARSING SUBROUTINES
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************
      MODULE ipch_dot

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      USE stel_kinds
      USE stel_constants
      USE ipch_T
      USE signal_dot
      IMPLICIT NONE
          
!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE ipch_dot_read
         MODULE PROCEDURE ipch_dot_read_old,                                   &
     &                    ipch_dot_read_new
      END INTERFACE

      INTERFACE ipch_dot_parse_chord
         MODULE PROCEDURE ipch_dot_parse_chord_old,                            &
     &                    ipch_dot_parse_chord_new
      END INTERFACE

      CONTAINS
!*******************************************************************************
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
!*******************************************************************************
!  Subroutine to parse the large structure of the ipch_dot_file
!  Individual ip chords are parsed by the ipch_dot_parse_xxxx subroutines
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
      CHARACTER(len=len_key), DIMENSION(1:n_ipchd_keyword) ::                  &
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
!  Subroutine to parse the large structure of the ipch_dot_file
!  Individual ip chords are parsed by the ipch_dot_parse_xxxx subroutines
!-------------------------------------------------------------------------------
      SUBROUTINE ipch_dot_read_new(ipch_file, signals,                         &
     &                             signals_created, observed, sigma,           &
     &                             weight, first_index)
      USE signal
      USE data_parameters

      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!  ipch_file        Name of the ipch. file.
!  signals          Array of pointers to the a signal
!  signals_created  Number valid signals in the array, or index of the last
!                   created.
!  observed         Array of observed values.
!  sigma            Array of observed sigma values
!  weight           Array of signal weight values
      CHARACTER (len=data_filename_length), INTENT(in)   :: ipch_file
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout) :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in) :: observed
      REAL (rprec), DIMENSION(:), INTENT(in) :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in) :: weight
      INTEGER, INTENT(inout)                 :: first_index

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!  n_ipchd_keyword        integer - number of ipchd_keyword s.
      INTEGER, PARAMETER :: n_ipchd_keyword = 9

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!  ipchd_keyword    character array - keywords for various ipch input types
      CHARACTER(len=len_key), DIMENSION(1:n_ipchd_keyword) ::                  &
     &   ipchd_keyword
      TYPE(signal_dot_file) :: ipch_dot_file_ref

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
! SECTION IV. PARSING SUBROUTINES
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
!-------------------------------------------------------------------------------
      SUBROUTINE ipch_dot_parse_chord_old(ipch_dot_file_ref,                   &
     &                                    coordinate_type,                     &
     &                                    coord_type,                          &
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
      CHARACTER (len=1) :: coord_type
      LOGICAL, INTENT(in) :: inDegrees
      TYPE(signal_dot_file), INTENT(inout) :: ipch_dot_file_ref

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(3)  :: xcart_i, xcart_f
      CHARACTER (len=c_name_len) :: chord_name                               
      REAL(rprec)                :: wavelength = 0.0
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Get the start, end and name of chord.
      CALL signal_dot_parse_chord(ipch_dot_file_ref,                           &
     &                            coordinate_type,                             &
     &                            chord_name,                                  &
     &                            xcart_i, xcart_f)

      IF (coord_type .eq. 'p') THEN
         wavelength = signal_dot_parse_real(ipch_dot_file_ref,                 &
     &                   "Expected wavelength for polarimetry chord")
      END IF

      icount_chords = icount_chords + 1
      IF (icount_chords .gt. n_ipch_desc_arr) THEN
         WRITE(*,*) 'TOO MANY ipCH. Increase dimensions'
      ELSE
         CALL ipch_desc_construct(ipch_desc_arr(icount_chords),                &
     &                            chord_name, xcart_i, xcart_f,                &
     &                            coord_type, wavelength, inDegrees)
      ENDIF

      END SUBROUTINE ipch_dot_parse_chord_old

!-------------------------------------------------------------------------------
!  Subroutine to parse an ipch chord (redundant chord there)
!      ipchd_keyword(1) = 'ip_chord_XYZ' or 'ip_chord_RPhiDegZ'
!-------------------------------------------------------------------------------
      SUBROUTINE ipch_dot_parse_chord_new(ipch_dot_file_ref,                   &
     &                                    coordinate_type,                     &
     &                                    coord_type,                          &
     &                                    inDegrees,                           &
     &                                    signals,                             &
     &                                    signals_created,                     &
     &                                    observed, sigma, weight,             &
     &                                    first_index) 
      USE signal
      USE v3fit_input, only: v3fit_max_diagnostics
      USE v3_utilities, only: err_fatal
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (signal_dot_file), INTENT(inout)  :: ipch_dot_file_ref
      CHARACTER (len=*), INTENT(in)          :: coordinate_type
      CHARACTER (len=1), INTENT(in)          :: coord_type
      LOGICAL, INTENT(in)                    :: inDegrees
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout)                 :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in) :: observed
      REAL (rprec), DIMENSION(:), INTENT(in) :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in) :: weight
      INTEGER, INTENT(inout)                 :: first_index

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL (rprec), DIMENSION(2,3) :: xcart
      CHARACTER (len=c_name_len)   :: chord_name
      REAL (rprec)                 :: wavelength = 0.0

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
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

      IF (coord_type .eq. 'p') THEN
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

!*******************************************************************************
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************

      END MODULE ipch_dot
