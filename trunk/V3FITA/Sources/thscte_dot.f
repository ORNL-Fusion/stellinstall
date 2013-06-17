!*******************************************************************************
!  File thscte_dot.f
!  Contains module thscte_dot
!    Module for opening and reading a 'thscte.' file, and then placing the
!     Thomson Scattering information into an array of thscte_desc
!
!  Based on module sxrch_dot
!  First version JDH 2011-10-24
!  
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       stel_constants
!       thscte_T
!       safe_open_mod
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!     Usage: Here is a sample of code that shows how to use this module
!
!      USE thscte_T
!      USE thscte_dot
!      TYPE(thscte_desc), DIMENSION(1000) :: thscte_descriptions_a
! ...
! !  Read in Thomson scattering information
!      CALL thscte_dot_read(file_name_thscte_dot,thscte_descriptions_a)
!
!-------------------------------------------------------------------------------
!   FILE STRUCTURE - thscte. file
!-------------------------------------------------------------------------------
!
!  thscte. File Structure:
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
!       thscted_keyword(1) = 'thscte_XYZ'
!       thscted_keyword(2) = 'thscte_RPhiDegZ'
!       thscted_keyword(3) = 'thscte_XYZ_te'
!       thscted_keyword(4) = 'thscte_RPhiDegZ_te'
!       thscted_keyword(5) = 'thscte_XYZ_ne'
!       thscted_keyword(6) = 'thscte_RPhiDegZ_ne'
!       thscted_keyword(7) = 'thscte_XYZ_pres'
!       thscted_keyword(8) = 'thscte_RPhiDegZ_pres'
!       thscted_keyword(9) = 'end_of_file'
!
!*******************************************************************************

!*******************************************************************************
!  MODULE thscte_dot
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
! SECTION IV. PARSING SUBROUTINES
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************
      MODULE thscte_dot

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      USE stel_kinds
      USE stel_constants
      USE thscte_T
      USE signal_dot
      IMPLICIT NONE
          
!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE thscte_dot_read
         MODULE PROCEDURE thscte_dot_read_old,                                 &
     &                    thscte_dot_read_new
      END INTERFACE

      INTERFACE thscte_dot_parse_chord
         MODULE PROCEDURE thscte_dot_parse_chord_old,                          &
     &                    thscte_dot_parse_chord_new
      END INTERFACE

      CONTAINS
!*******************************************************************************
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
!*******************************************************************************
!  Subroutine to parse the large structure of the thscte_dot_file
!  Individual sxr chords are parsed by the thscte_dot_parse_xxxx subroutines
!-------------------------------------------------------------------------------
      SUBROUTINE thscte_dot_read_old(thscte_file,thscte_desc_arr,              &
     &                               n_points)

      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in)            :: thscte_file
      TYPE (thscte_desc), DIMENSION(:), INTENT(inout) :: thscte_desc_arr
      INTEGER, INTENT(inout)  :: n_points


!  thscte_file            Name of the thscte. file.
!  thscte_desc_arr        Array of thscte_desc. Assumed to have been allocated
!                        before call to thscte_dot_read
!  n_points              Number of points read from the file

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n_thscted_keyword = 9
!  n_thscted_keyword    integer - number of thscted_keyword s.

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=len_key), DIMENSION(1:n_thscted_keyword) ::                &
     &   thscted_keyword

      INTEGER               :: icount_points
      INTEGER               :: n_thscte_desc_arr
      CHARACTER (len=30)    :: s_name                                 
      CHARACTER (len=80)    :: l_name
      CHARACTER (len=len_key) :: this_keyword

      TYPE(signal_dot_file) :: thscte_dot_file_ref
!  Make sure that len for this_keyword is the same as len for thscted_keyword.

!  thscted_keyword   character array - keywords for various thscte input types

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Initialize the keywords and keyword lengths
      thscted_keyword(1) = 'thscte_XYZ'
      thscted_keyword(2) = 'thscte_RPhiDegZ'
      thscted_keyword(3) = 'thscte_XYZ_te'
      thscted_keyword(4) = 'thscte_RPhiDegZ_te'
      thscted_keyword(5) = 'thscte_XYZ_ne'
      thscted_keyword(6) = 'thscte_RPhiDegZ_ne'
      thscted_keyword(7) = 'thscte_XYZ_pres'
      thscted_keyword(8) = 'thscte_RPhiDegZ_pres'
      thscted_keyword(9) = 'end_of_file'

!  Find the length of the thscte_desc_arr array
      n_thscte_desc_arr = SIZE(thscte_desc_arr)

!  Open up the 'thscte.' file
      thscte_dot_file_ref = signal_dot_open(TRIM(thscte_file),                 &
     &                                      "thscte")

      icount_points = 0
!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(thscte_dot_file_ref,             &
     &                                        thscted_keyword))
         
            CASE DEFAULT ! This case should never fire.
               EXIT ! Exit out of infinte loop.

            CASE ('end_of_file')
               EXIT ! Exit out of infinte loop.
            
            CASE ('thscte_XYZ','thscte_XYZ_te')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'XYZ', 't', icount_points,          &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)
               
            CASE ('thscte_RPhiDegZ','thscte_RPhiDegZ_te')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'RPHiDegZ', 't',                    &
     &                                     icount_points,                      &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)

            CASE ('thscte_XYZ_ne')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'XYZ', 'd', icount_points,          &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)

            CASE ('thscte_RPhiDegZ_ne')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'RPHiDegZ', 'd',                    &
     &                                     icount_points,                      &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)

            CASE ('thscte_XYZ_pres')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'XYZ', 'p', icount_points,          &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)

            CASE ('thscte_RPhiDegZ_pres')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'RPHiDegZ', 'p',                    &
     &                                     icount_points,                      &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)

         END SELECT
         
!  Check for too many additions to the thscte_desc_arr array
         IF (icount_points .gt. SIZE(thscte_desc_arr)) THEN
            STOP 'thscte_dot: Too many coils for thscte_desc_arr '
         END IF
         
      END DO

!  Close the 'thscte.' file
      CALL signal_dot_close(thscte_dot_file_ref)
      
      n_points = icount_points
      
      RETURN
      END SUBROUTINE thscte_dot_read_old

!-------------------------------------------------------------------------------
!  Subroutine to parse the large structure of the thscte_dot_file
!  Individual sxr chords are parsed by the thscte_dot_parse_xxxx subroutines
!-------------------------------------------------------------------------------
      SUBROUTINE thscte_dot_read_new(thscte_file, signals,                     &
     &                               signals_created, observed, sigma,         &
     &                               weight, first_index)
      USE signal
      USE data_parameters

      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!  thscte_file      Name of the thscte. file.
!  signals          Array of pointers to the a signal
!  signals_created  Number valid signals in the array, or index of the last
!                   created.
!  observed         Array of observed values.
!  sigma            Array of observed sigma values
!  weight           Array of signal weight values
      CHARACTER (len=data_filename_length), INTENT(in)   :: thscte_file
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout) :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in) :: observed
      REAL (rprec), DIMENSION(:), INTENT(in) :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in) :: weight
      INTEGER, INTENT(inout)                 :: first_index

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!  n_thscted_keyword    integer - number of thscted_keyword s.
      INTEGER, PARAMETER :: n_thscted_keyword = 9

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!  thscted_keyword   character array - keywords for various thscte input types
      CHARACTER(len=len_key), DIMENSION(1:n_thscted_keyword) ::                &
     &   thscted_keyword
      TYPE(signal_dot_file) :: thscte_dot_file_ref

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Initialize the keywords and keyword lengths
      thscted_keyword(1) = 'thscte_XYZ'
      thscted_keyword(2) = 'thscte_RPhiDegZ'
      thscted_keyword(3) = 'thscte_XYZ_te'
      thscted_keyword(4) = 'thscte_RPhiDegZ_te'
      thscted_keyword(5) = 'thscte_XYZ_ne'
      thscted_keyword(6) = 'thscte_RPhiDegZ_ne'
      thscted_keyword(7) = 'thscte_XYZ_pres'
      thscted_keyword(8) = 'thscte_RPhiDegZ_pres'
      thscted_keyword(9) = 'end_of_file'

!  Open up the 'thscte.' file
      thscte_dot_file_ref = signal_dot_open(TRIM(thscte_file),                 &
     &                                      "thscte")

!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(thscte_dot_file_ref,             &
     &                                        thscted_keyword))

            CASE DEFAULT ! This case should never fire.
               EXIT ! Exit out of infinte loop.

            CASE ('end_of_file')
               EXIT ! Exit out of infinte loop.

            CASE ('thscte_XYZ','thscte_XYZ_te')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'XYZ', 't',                         &
     &                                     signals, signals_created,           &
     &                                     observed, sigma, weight,            &
     &                                     first_index)

            CASE ('thscte_RPhiDegZ','thscte_RPhiDegZ_te')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'RPHiDegZ', 't',                    &
     &                                     signals, signals_created,           &
     &                                     observed, sigma, weight,            &
     &                                     first_index)

            CASE ('thscte_XYZ_ne')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'XYZ', 'd',                         &
     &                                     signals, signals_created,           &
     &                                     observed, sigma, weight,            &
     &                                     first_index)

            CASE ('thscte_RPhiDegZ_ne')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'RPHiDegZ', 'd',                    &
     &                                     signals, signals_created,           &
     &                                     observed, sigma, weight,            &
     &                                     first_index)

            CASE ('thscte_XYZ_pres')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'XYZ', 'p',                         &
     &                                     signals, signals_created,           &
     &                                     observed, sigma, weight,            &
     &                                     first_index)

            CASE ('thscte_RPhiDegZ_pres')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'RPHiDegZ', 'p',                    &
     &                                     signals, signals_created,           &
     &                                     observed, sigma, weight,            &
     &                                     first_index)

         END SELECT

      END DO

!  Close the 'thscte.' file
      CALL signal_dot_close(thscte_dot_file_ref)

      END SUBROUTINE thscte_dot_read_new

!*******************************************************************************
! SECTION IV. PARSING SUBROUTINES
!*******************************************************************************

!-------------------------------------------------------------------------------
!     Expected behaviour for the parsing routines:
!  On entry:
!     - character variable 'line' is defined
!     - logical variable thscted_l_done is .false.
!     - integer variable i_line is set
!  During Execution:
!     - integer variable i_line is kept current
!     - reads lines from file, looks for appropriate line structure
!     - presence of inappropriate line indicates end of chord
!     - constructs thscte_desc_temp
!  On exit:
!     - 'line' contains the first inappropriate line
!     - End of File or other problem: thscted_l_done set to .true.
!     
!-------------------------------------------------------------------------------
!  Subroutine to parse an thscte chord (redundant chord there)
!      thscted_keyword(1) = 'thscte_XYZ' or 'thscte_RPhiDegZ'
!-------------------------------------------------------------------------------
      SUBROUTINE thscte_dot_parse_chord_old(thscte_dot_file_ref,               &
     &                                     coordinate_type,                    &
     &                                     point_type,                         &
     &                                     icount_points,                      &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(inout) :: icount_points
      INTEGER, INTENT(in)    :: n_thscte_desc_arr
      TYPE(thscte_desc), DIMENSION(:), INTENT(inout) :: thscte_desc_arr
      CHARACTER (len=*) :: coordinate_type
      CHARACTER (len=1) :: point_type
      TYPE(signal_dot_file), INTENT(inout) :: thscte_dot_file_ref

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(3)  :: xcart
      CHARACTER (len=c_name_len) :: chord_name
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Get the start, end and name of chord.
      CALL signal_dot_parse_chord(thscte_dot_file_ref,                         &
     &                            coordinate_type,                             &
     &                            chord_name,                                  &
     &                            xcart)

      icount_points = icount_points + 1
      IF (icount_points .gt. n_thscte_desc_arr) THEN
         WRITE(*,*) 'TOO MANY thscte. Increase dimensions'
      ELSE
         CALL thscte_desc_construct(thscte_desc_arr(icount_points),            &
     &                              chord_name, xcart, point_type)
      ENDIF

      END SUBROUTINE thscte_dot_parse_chord_old

!-------------------------------------------------------------------------------
!  Subroutine to parse an thscte chord (redundant chord there)
!      thscted_keyword(1) = 'thscte_XYZ' or 'thscte_RPhiDegZ'
!-------------------------------------------------------------------------------
      SUBROUTINE thscte_dot_parse_chord_new(thscte_dot_file_ref,               &
     &                                      coordinate_type,                   &
     &                                      point_type,                        &
     &                                      signals,                           &
     &                                      signals_created,                   &
     &                                      observed, sigma, weight,           &
     &                                      first_index)
      USE signal
      USE v3fit_input, only: v3fit_max_diagnostics
      USE v3_utilities, only: err_fatal
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE(signal_dot_file), INTENT(inout)   :: thscte_dot_file_ref
      CHARACTER (len=*)                      :: coordinate_type
      CHARACTER (len=1)                      :: point_type
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      INTEGER, INTENT(inout)                 :: signals_created
      REAL (rprec), DIMENSION(:), INTENT(in) :: observed
      REAL (rprec), DIMENSION(:), INTENT(in) :: sigma
      REAL (rprec), DIMENSION(:), INTENT(in) :: weight
      INTEGER, INTENT(inout)                 :: first_index

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL (rprec), DIMENSION(3)             :: xcart
      CHARACTER (len=c_name_len)             :: point_name
      TYPE (thomson_class), POINTER          :: thomson_obj => null()
      CHARACTER (len=signal_units_length)    :: units = ''

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Check if there are two many signals.
      IF (signals_created + 1 .gt. v3fit_max_diagnostics) THEN
         CALL err_fatal('ipch_dot_parse_chord_new: created signals' //         &
     &                  ' exceeds v3fit_max_diagnostics')
      END IF

!  Get the position and name of the point.
      CALL signal_dot_parse_chord(thscte_dot_file_ref,                         &
     &                            coordinate_type,                             &
     &                            point_name,                                  &
     &                            xcart)

      SELECT CASE (point_type)

         CASE ('t')
            thomson_obj => thomson_construct_te(xcart)
            units = 'ev'

         CASE ('d')
            thomson_obj => thomson_construct_ne(xcart)
            units = 'm^-3'

         CASE ('p')
            thomson_obj => thomson_construct_p(xcart)
            units = 'Pa'

      END SELECT

      signals(signals_created + 1)%p => signal_construct(                      &
     &   diagnostic_construct(thomson_obj), point_name, point_name,            &
     &   units, observed(signals_created + 1),                                 &
     &   sigma(signals_created + 1), weight(signals_created + 1))

      signals_created = signals_created + 1

!  At lease one intpol signal was made. Set the first index. This should only be
!  run once.
      IF (first_index .eq. -1) THEN
         first_index = signals_created
      END IF

      END SUBROUTINE thscte_dot_parse_chord_new

!*******************************************************************************
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************

      END MODULE thscte_dot
