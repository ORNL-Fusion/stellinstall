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
!       ipchd_keyword(3) = 'end_of_file' 
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

!-------------------------------------------------------------------------------
!  Variables local to the module. Make them private.
!-------------------------------------------------------------------------------
      INTEGER ::  ipchd_iou
      INTEGER :: i_line
      CHARACTER(len=200) :: line
      LOGICAL :: ipchd_l_done
      TYPE (ipch_desc), SAVE :: ipch_desc_temp

      PRIVATE ipchd_iou,i_line 
      PRIVATE line, ipchd_l_done, ipch_desc_temp

!  ipch_desc_temp    A temporary ipch_desc. SAVE attribute - see note in 
!                         diagostic_dot about Lahey Fujitsu compiler
!  ipchd_iou          i/o unit number for ipch_dot file
!  ipchd_l_done       logical - true if done reading file
!  i_line             integer - line counter for file
!  line               character - the last line read from the file
          
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

      SUBROUTINE ipch_dot_read (ipch_file,ipch_desc_arr,n_chords)
      
!  Subroutine to parse the large structure of the ipch_dot_file
!  Individual ip chords are parsed by the ipch_dot_parse_xxxx subroutines

      USE safe_open_mod
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
      INTEGER, PARAMETER :: n_ipchd_keyword = 3
!  n_ipchd_keyword        integer - number of ipchd_keyword s.

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=len_key), DIMENSION(1:n_ipchd_keyword) ::                  &
     &   ipchd_keyword
      INTEGER, DIMENSION(1:n_ipchd_keyword) :: n_char_key

      INTEGER               :: icount_chords
      INTEGER               :: n_ipch_desc_arr
      CHARACTER (len=30)    :: s_name                                 
      CHARACTER (len=80)    :: l_name
!  Make sure that len for this_keyword is the same as len for ipchd_keyword.

!  ipchd_keyword    character array - keywords for various ipch input types
!  n_char_key        integer array - lengths of the ipchd_keyword strings

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Initialize the keywords and keyword lengths
      ipchd_keyword(1) = 'ip_chord_XYZ'
      ipchd_keyword(2) = 'ip_chord_RPhiDegZ'
      ipchd_keyword(3) = 'end_of_file'
      n_char_key(1:n_ipchd_keyword) =                                          &
     &   len_trim(ipchd_keyword(1:n_ipchd_keyword))

!  Find the length of the ipch_desc_arr array
      n_ipch_desc_arr = SIZE(ipch_desc_arr)

!  Open up the 'ipch.' file
      ipchd_iou = signal_dot_open(TRIM(ipch_file), "ipch", line)
      i_line = 2
      ipchd_l_done = .false.

      icount_chords = 0

!  Infinite Loop
!  Character variable line should be defined on entry
      DO
!  Branch on the keyword
         SELECT CASE (signal_dot_read_keyword(n_ipchd_keyword,                 &
     &                                        n_char_key, line,                &
     &                                        ipchd_keyword))
         
            CASE DEFAULT ! Did  not find a keyword. Read another line.
               ipchd_l_done = signal_dot_check_next_line(ipchd_iou,            &
     &                                                   i_line,               &
     &                                                   line)
            
            CASE ('ip_chord_XYZ')
               CALL ipch_dot_parse_chord('XYZ', icount_chords,                 &
     &                                   n_ipch_desc_arr,                      &
     &                                   ipch_desc_arr)
               
            CASE ('ip_chord_RPhiDegZ')
               CALL ipch_dot_parse_chord('RPHiDegZ', icount_chords,            &
     &                                   n_ipch_desc_arr,                      &
     &                                   ipch_desc_arr)
                           
            CASE ('end_of_file')
               ipchd_l_done = .true.

         END SELECT
         
!  Check for too many additions to the ipch_desc_arr array
         IF  (icount_chords .gt. SIZE(ipch_desc_arr)) THEN
            STOP 'ipch_dot: Too many coils for ipch_desc_arr '
         ENDIF
        
         IF (ipchd_l_done) EXIT ! the infinite loop
         
      END DO

!  Clean up
      CALL ipch_desc_destroy(ipch_desc_temp)

!  Close the 'ipch.' file
      CLOSE(ipchd_iou)
      
      n_chords = icount_chords
      
      RETURN
      END SUBROUTINE ipch_dot_read

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

      SUBROUTINE ipch_dot_parse_chord(coordinate_type,                         &
     &                                icount_chords,                           &
     &                                n_ipch_desc_arr,                         &
     &                                ipch_desc_arr)
!  Subroutine to parse an ipch chord (redundant chord there)
!      ipchd_keyword(1) = 'ip_chord_XYZ' or 'ip_chord_RPhiDegZ'
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(inout)                       :: icount_chords
      INTEGER, INTENT(in)                          :: n_ipch_desc_arr
      TYPE(ipch_desc), DIMENSION(:), INTENT(inout) :: ipch_desc_arr
      CHARACTER (len=*) :: coordinate_type
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(3)   :: xcart_i, xcart_f
      CHARACTER (len=30)          :: chord_name                               
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Get the start, end and name of chord.
      CALL signal_dot_parse_chord(ipchd_iou,                                   &
     &                            coordinate_type,                             &
     &                            "ip",                                        &
     &                            line, i_line,                                &
     &                            chord_name,                                  &
     &                            xcart_i, xcart_f)
      CALL ipch_desc_construct(ipch_desc_temp,chord_name,                      &
     &                         xcart_i,xcart_f)

!  Read the next line, in preparation for returning to main parsing loop
!-----------------------------------------------
      ipchd_l_done = signal_dot_check_next_line(ipchd_iou, i_line, line)

      icount_chords = icount_chords + 1
      IF (icount_chords .gt. n_ipch_desc_arr) THEN
         WRITE(*,*) 'TOO MANY ipCH. Increase dimensions'
      ELSE
         ipch_desc_arr(icount_chords) = ipch_desc_temp
      ENDIF

      END SUBROUTINE ipch_dot_parse_chord


!*******************************************************************************
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************

      END MODULE ipch_dot
