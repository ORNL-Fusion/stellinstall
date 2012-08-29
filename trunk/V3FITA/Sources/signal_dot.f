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

      PUBLIC signal_dot_parse_chord
      PRIVATE signal_dot_check_status

      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Parse a chord from a diagnostic dot file.
!*******************************************************************************
      SUBROUTINE signal_dot_parse_chord(iou, istat,                            &
     &                                  coordinate_type,                       &
     &                                  signal_id,                             &
     &                                  line, i_line,                          &
     &                                  xcart_i, xcart_f,                      &
     &                                  chord_name)
      USE coordinate_utilities
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(inout)                 :: iou, istat
      CHARACTER (len=*), INTENT(in)          :: signal_id
      CHARACTER (len=200), INTENT(inout)     :: line
      INTEGER, INTENT(inout)                 :: i_line
      REAL(rprec), DIMENSION(3), INTENT(out) :: xcart_i, xcart_f
      CHARACTER (len=30), INTENT(out)        :: chord_name

!-------------------------------------------------------------------------------
!  Dummy Arguments
!-------------------------------------------------------------------------------
      CHARACTER (len=*) :: coordinate_type

!-------------------------------------------------------------------------------
!   Local Variables
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(3) :: x_start, x_end

!-------------------------------------------------------------------------------
!  Start of Executable Code
!-------------------------------------------------------------------------------

! Parse beyond the keyword in the first line.
! Read the second line, chord identification.
      READ (iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      CALL signal_dot_check_status(istat, i_line,                              &
     &                             "Error reading ID for ",                    &
     &                             " chord, line ",                            &
     &                             signal_id)

      chord_name = TRIM(line)

! Read the third line, start position.
      READ (iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      CALL signal_dot_check_status(istat, i_line,                              &
     &                             "Error ",                                   &
     &                             " start position line",                     &
     &                             signal_id)

! Reread the line, looking for three reals
      READ(line, *, iostat=istat) x_start(1:3)
      CALL signal_dot_check_status(istat, i_line,                              &
     &                             "Error ",                                   &
     &                             " start position",                          &
     &                             signal_id)

! Read the fourth line, end position.
      READ (iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      CALL signal_dot_check_status(istat, i_line,                              &
     &                             "Error ",                                   &
     &                             " end position line",                       &
     &                             signal_id)

! Reread the line, looking for three reals.
      READ(line, *, iostat=istat) x_end(1:3)
      CALL signal_dot_check_status(istat, i_line,                              &
     &                             "Error ",                                   &
     &                             " end position",                            &
     &                             signal_id)

! Done reading information for this sxrch. Generate the sxrch
      IF (TRIM(coordinate_type) .eq. 'RPHiDegZ') THEN
         x_start(2) = x_start(2)*pi/180. ! Convert to degrees to radians.
         xcart_i = cyl_to_cart(x_start)

         x_end(2) = x_end(2)*pi/180.     ! Convert to degrees to radians.
         xcart_f = cyl_to_cart(x_end)
      ELSE
         xcart_i = x_start
         xcart_f = x_end
      ENDIF

      END SUBROUTINE

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
! Check and handle errors.
!*******************************************************************************
      SUBROUTINE signal_dot_check_status(istat,                                &
     &                                   i_line,                               &
     &                                   pre_message,                          &
     &                                   post_message,                         &
     &                                   signal_id)
      IMPLICIT NONE
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in)           :: istat, i_line
      CHARACTER (len=*), INTENT(in) :: pre_message, post_message,              &
     &                                 signal_id

!-------------------------------------------------------------------------------
!  Start of Executable Code
!-------------------------------------------------------------------------------
      IF (istat .ne. 0) THEN
         WRITE(*,*) pre_message // signal_id // post_message, i_line
         STOP
      END IF

      END SUBROUTINE

      END MODULE
