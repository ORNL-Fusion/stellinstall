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
!       thscted_keyword(3) = 'end_of_file' 
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

!-------------------------------------------------------------------------------
!  Variables local to the module. Make them private.
!-------------------------------------------------------------------------------
      TYPE (thscte_desc), SAVE :: thscte_desc_temp
      PRIVATE thscte_desc_temp

!  thscte_desc_temp     A temporary thscte_desc. SAVE attribute - see note in 
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

      SUBROUTINE thscte_dot_read (thscte_file,thscte_desc_arr,n_points)
      
!  Subroutine to parse the large structure of the thscte_dot_file
!  Individual sxr chords are parsed by the thscte_dot_parse_xxxx subroutines

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
      INTEGER, PARAMETER :: n_thscted_keyword = 3
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
      thscted_keyword(3) = 'end_of_file'

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
            
            CASE ('thscte_XYZ')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'XYZ', icount_points,               &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)
               
            CASE ('thscte_RPhiDegZ')
               CALL thscte_dot_parse_chord(thscte_dot_file_ref,                &
     &                                     'RPHiDegZ', icount_points,          &
     &                                     n_thscte_desc_arr,                  &
     &                                     thscte_desc_arr)

         END SELECT
         
!  Check for too many additions to the thscte_desc_arr array
         IF  (icount_points .gt. SIZE(thscte_desc_arr)) THEN
            STOP 'thscte_dot: Too many coils for thscte_desc_arr '
         ENDIF
         
      END DO

!  Clean up
      CALL thscte_desc_destroy(thscte_desc_temp)

!  Close the 'thscte.' file
      CALL signal_dot_close(thscte_dot_file_ref)
      
      n_points = icount_points
      
      RETURN
      END SUBROUTINE thscte_dot_read

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

      SUBROUTINE thscte_dot_parse_chord(thscte_dot_file_ref,                   &
     &                                  coordinate_type,                       &
     &                                  icount_points,                         &
     &                                  n_thscte_desc_arr,                     &
     &                                  thscte_desc_arr)
!  Subroutine to parse an thscte chord (redundant chord there)
!      thscted_keyword(1) = 'thscte_XYZ' or 'thscte_RPhiDegZ'
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(inout) :: icount_points
      INTEGER, INTENT(in)    :: n_thscte_desc_arr
      TYPE(thscte_desc), DIMENSION(:), INTENT(inout) :: thscte_desc_arr
      CHARACTER (len=*) :: coordinate_type
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
      CALL thscte_desc_construct(thscte_desc_temp,chord_name,                  &
     &                           xcart)

      icount_points = icount_points + 1
      IF (icount_points .gt. n_thscte_desc_arr) THEN
         WRITE(*,*) 'TOO MANY thscte. Increase dimensions'
      ELSE
         thscte_desc_arr(icount_points) = thscte_desc_temp
      ENDIF

      END SUBROUTINE thscte_dot_parse_chord

!*******************************************************************************
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************

      END MODULE thscte_dot
