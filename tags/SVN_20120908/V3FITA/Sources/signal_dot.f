!*******************************************************************************
!  File signal_dot.f
!  Contains module signal_dot
!    Module containing common code used in parsing diagnostic dot file.
!
!*******************************************************************************

      MODULE signal_dot
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

      TYPE signal_dot_file
         PRIVATE
         INTEGER :: iou
         INTEGER :: i_line
         CHARACTER (len=6) :: signal_id
      END TYPE

      PUBLIC  :: signal_dot_open, signal_dot_close,
     &           signal_dot_parse_chord, signal_dot_read_keyword
      PRIVATE :: signal_dot_check_status, signal_dot_read_line

!-------------------------------------------------------------------------------
!  Module Parameters
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER :: len_key = 30, c_name_len = 30
      INTEGER, PARAMETER :: line_len = 200

      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Open diagnostic dot file.
!*******************************************************************************
      FUNCTION signal_dot_open(file, signal_id)
      USE safe_open_mod
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER*(*), INTENT(in)     :: file
      CHARACTER (len=*), INTENT(in) :: signal_id
      TYPE(signal_dot_file)         :: signal_dot_open

!-------------------------------------------------------------------------------
!   Local Variables
!-------------------------------------------------------------------------------
      INTEGER :: istat

!-------------------------------------------------------------------------------
!  Local Parameters
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER :: iou_0=32

!-------------------------------------------------------------------------------
!  Start of Executable Code
!-------------------------------------------------------------------------------
!  Open up the file
      signal_dot_open%iou       = iou_0
      signal_dot_open%signal_id = signal_id
      CALL safe_open(signal_dot_open%iou,istat,file,'old','formatted')
      CALL signal_dot_check_status(istat, signal_dot_open,                     &
     &                             "Error opening ",                           &
     &                             " file.")

      END FUNCTION

!*******************************************************************************
!  Close diagnostic dot file.
!*******************************************************************************
      SUBROUTINE signal_dot_close(signal_dot_file_ref)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE(signal_dot_file) :: signal_dot_file_ref

!-------------------------------------------------------------------------------
!  Start of Executable Code
!-------------------------------------------------------------------------------
!  Close the file
      CLOSE(signal_dot_file_ref%iou)

      END SUBROUTINE

!*******************************************************************************
!  Parse a chord from a diagnostic dot file.
!
!  Chords start with the first line as the chord name. The it is followed by one
!  or two lines of 3 coordinates.
!*******************************************************************************
      SUBROUTINE signal_dot_parse_chord(signal_dot_file_ref,                   &
     &                                  coordinate_type,                       &
     &                                  chord_name,                            &
     &                                  xcart_i, xcart_f)
      USE coordinate_utilities
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE(signal_dot_file), INTENT(inout) :: signal_dot_file_ref
      CHARACTER (len=*), INTENT (in) :: coordinate_type
      CHARACTER (len=c_name_len), INTENT(out)          :: chord_name
      REAL(rprec), DIMENSION(3), INTENT(out)           :: xcart_i
      REAL(rprec), DIMENSION(3), OPTIONAL, INTENT(out) :: xcart_f

!-------------------------------------------------------------------------------
!   Local Variables
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(3) :: x_start, x_end
      INTEGER                   :: istat
      CHARACTER (len=line_len)  :: line

!-------------------------------------------------------------------------------
!  Start of Executable Code
!-------------------------------------------------------------------------------
! Read the chord identification.
      line = signal_dot_read_line(signal_dot_file_ref,                         &
     &                            "Error reading ID for ",                     &
     &                            " chord.")
      chord_name = TRIM(line)

! Read the start position.
      line = signal_dot_read_line(signal_dot_file_ref,                         &
     &                            "Error reading start position for ",         &
     &                            " chord.")
! Reread the line, looking for three reals
      READ(line, *, iostat=istat) x_start(1:3)
      CALL signal_dot_check_status(istat, signal_dot_file_ref,                 &
     &                             "Error ",                                   &
     &                             " chord start position.")

      IF (PRESENT(xcart_f)) THEN
! Read the end position.
         line = signal_dot_read_line(signal_dot_file_ref,                      &
     &                               "Error reading end position for ",        &
     &                               " chord.")
! Reread the line, looking for three reals
         READ(line, *, iostat=istat) x_end(1:3)
         CALL signal_dot_check_status(istat, signal_dot_file_ref,              &
     &                                "Error ",                                &
     &                                " chord end position.")
      ENDIF

! Done reading information for this chord. Generate the chord
      IF (TRIM(coordinate_type) .eq. 'RPHiDegZ') THEN
         x_start(2) = x_start(2)*degree ! Convert from degrees to radians.
         xcart_i = cyl_to_cart(x_start)

         IF (PRESENT(xcart_f)) THEN
            x_end(2) = x_end(2)*degree  ! Convert from degrees to radians.
            xcart_f = cyl_to_cart(x_end)
         ENDIF
      ELSE
         xcart_i = x_start
         IF (PRESENT(xcart_f)) THEN
            xcart_f = x_end
         ENDIF
      ENDIF

      END SUBROUTINE

!*******************************************************************************
!  Read a keyword from the diagnostic dot file.
!
!  Read lines until a valid keyword is reached ot the end of the file is
!  reached.
!*******************************************************************************
      FUNCTION signal_dot_read_keyword(signal_dot_file_ref, keywords)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER (len=len_key), DIMENSION(:), INTENT(in) :: keywords
      CHARACTER (len=len_key) :: signal_dot_read_keyword
      TYPE(signal_dot_file), INTENT(inout) :: signal_dot_file_ref

!-------------------------------------------------------------------------------
!   Local Variables
!-------------------------------------------------------------------------------
      INTEGER                  :: i
      CHARACTER (len=line_len) :: line

!-------------------------------------------------------------------------------
!  Start of Executable Code
!-------------------------------------------------------------------------------
! Contine looping until a vaild keyword is read. If the end of the file is
! reached an error will kill the program.
      DO
         line = signal_dot_read_line(signal_dot_file_ref,                      &
     &                               "Failed to read ",                        &
     &                               " keyword before end of file.")
         DO i = 1, SIZE(keywords)
             IF (TRIM(ADJUSTL(line)) .eq. TRIM(keywords(i))) THEN
                signal_dot_read_keyword = keywords(i)
                RETURN
             ENDIF
         END DO
      END DO

      END FUNCTION

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Read line from diagnostic dot file.
!*******************************************************************************
      FUNCTION signal_dot_read_line(signal_dot_file_ref,                       &
     &                              pre_message,                               &
     &                              post_message)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER (len=*), INTENT(in) :: pre_message, post_message
      CHARACTER (len=line_len) :: signal_dot_read_line
      TYPE(signal_dot_file), INTENT(inout) :: signal_dot_file_ref

!-------------------------------------------------------------------------------
!   Local Variables
!-------------------------------------------------------------------------------
      INTEGER :: istat

!-------------------------------------------------------------------------------
!  Start of Executable Code
!-------------------------------------------------------------------------------
      READ (signal_dot_file_ref%iou, '(a)', iostat=istat)                      &
     &   signal_dot_read_line
      signal_dot_file_ref%i_line = signal_dot_file_ref%i_line + 1
      CALL signal_dot_check_status(istat, signal_dot_file_ref,                 &
     &                             pre_message, post_message)

      END FUNCTION

!*******************************************************************************
!  Check and handle errors.
!*******************************************************************************
      SUBROUTINE signal_dot_check_status(istat,                                &
     &                                   signal_dot_file_ref,                  &
     &                                   pre_message,                          &
     &                                   post_message)
      IMPLICIT NONE
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in)           :: istat
      TYPE(signal_dot_file), INTENT(inout) :: signal_dot_file_ref
      CHARACTER (len=*), INTENT(in) :: pre_message, post_message

!-------------------------------------------------------------------------------
!  Start of Executable Code
!-------------------------------------------------------------------------------
      IF (istat .ne. 0) THEN
         WRITE(*,*) pre_message//TRIM(signal_dot_file_ref%signal_id)//         &
     &              post_message//" line ",signal_dot_file_ref%i_line
         STOP
      END IF

      END SUBROUTINE

      END MODULE
