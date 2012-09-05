!*******************************************************************************
!  File sxrch_dot.f
!  Contains module sxrch_dot
!    Module for opening and reading a 'sxrch.' file, and then placing the
!     soft x-ray chord information into an array of sxrch_desc
!
!  Based on module sxrch_dot
!  First version JDH 2011-10-14
!  
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       stel_constants
!       sxrch_T
!       safe_open_mod
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  Based on module sxrch_dot
!  First version JDH 2011-10-14
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!     Usage: Here is a sample of code that shows how to use this module
!
!      USE sxrch_T
!      USE sxrch_dot
!      TYPE(sxrch_desc), DIMENSION(1000) :: sxrch_descriptions_a
! ...
! !  Read in soft x-ray chords
!      CALL sxrch_dot_read(file_name_sxrch_dot,sxrch_descriptions_a)
!
!-------------------------------------------------------------------------------
!   FILE STRUCTURE - sxrch. file
!-------------------------------------------------------------------------------
!
!  sxrch. File Structure:
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
!       sxrchd_keyword(1) = 'sxr_chord_XYZ'
!       sxrchd_keyword(2) = 'sxr_chord_RPhiDegZ'
!       sxrchd_keyword(3) = 'end_of_file' 
!
!*******************************************************************************

!*******************************************************************************
!  MODULE sxrch_dot
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
! SECTION IV. PARSING SUBROUTINES
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************
      MODULE sxrch_dot

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      USE stel_kinds
      USE stel_constants
      USE sxrch_T
      USE signal_dot
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variables local to the module. Make them private.
!-------------------------------------------------------------------------------
      TYPE (sxrch_desc), SAVE :: sxrch_desc_temp
      PRIVATE sxrch_desc_temp

!  sxrch_desc_temp    A temporary sxrch_desc. SAVE attribute - see note in 
!                         diagostic_dot about Lahey Fujitsu compiler
          
!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
!*******************************************************************************

      SUBROUTINE sxrch_dot_read (sxrch_file,sxrch_desc_arr,n_chords)
      
!  Subroutine to parse the large structure of the sxrch_dot_file
!  Individual sxr chords are parsed by the sxrch_dot_parse_xxxx subroutines

      USE safe_open_mod
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
      CHARACTER(len=len_key), DIMENSION(1:n_sxrchd_keyword) ::                 &
     &   sxrchd_keyword

      INTEGER               :: icount_chords
      INTEGER               :: n_sxrch_desc_arr
      CHARACTER (len=30)    :: s_name                                 
      CHARACTER (len=80)    :: l_name

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

!  Clean up
      CALL sxrch_desc_destroy(sxrch_desc_temp)

!  Close the 'sxrch.' file
      CALL signal_dot_close(sxrch_dot_file_ref)
      
      n_chords = icount_chords
      
      RETURN
      END SUBROUTINE sxrch_dot_read

!*******************************************************************************
! SECTION IV. PARSING SUBROUTINES
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

      SUBROUTINE sxrch_dot_parse_chord(sxrch_dot_file_ref,                     &
     &                                 coordinate_type,                        &
     &                                 icount_chords,                          &
     &                                 n_sxrch_desc_arr,                       &
     &                                 sxrch_desc_arr)
!  Subroutine to parse an sxrch chord (redundant chord there)
!      sxrchd_keyword(1) = 'sxr_chord_XYZ' or 'sxr_chord_RPhiDegZ'
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
      CALL sxrch_desc_construct(sxrch_desc_temp,chord_name,                    &
     &                          xcart_i,xcart_f)

      icount_chords = icount_chords + 1
      IF (icount_chords .gt. n_sxrch_desc_arr) THEN
         WRITE(*,*) 'TOO MANY SXRCH. Increase dimnesions'
      ELSE
         sxrch_desc_arr(icount_chords) = sxrch_desc_temp
      ENDIF

      END SUBROUTINE sxrch_dot_parse_chord

!*******************************************************************************
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************

      END MODULE sxrch_dot
