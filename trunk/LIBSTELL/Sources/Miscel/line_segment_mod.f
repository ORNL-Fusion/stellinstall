!*******************************************************************************
!  File line_segment.f
!
!  Module is part of the libstell.a. This module containes code to create
!  profile constructed of line sigments. These line segments are assumed to be
!  specified such that xx(i) < xx(i + 1)
!
!*******************************************************************************

      MODULE line_segment
      USE stel_kinds
      USE stel_constants

      PUBLIC line_seg, line_seg_int, line_seg_test
      PRIVATE get_indices, y_value, y_value_int, slope, offset

      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Find the linearly interpolated value of y at position x from line sigments
!  specified by yy(xx).
!*******************************************************************************
      SUBROUTINE line_seg(x, y, xx, yy, n, iflag)
      IMPLICIT none

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: n
      INTEGER, INTENT(out) :: iflag

      REAL(rprec), INTENT(out) :: y
      REAL(rprec), INTENT(in)  :: x
      REAL(rprec), DIMENSION(n), INTENT(in)  :: xx, yy
!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER :: ilow, ihigh

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Check for valid line segments.
      IF (n .le. 1) THEN
         iflag = -1 ! Line sigments require at least two points
         RETURN
      ELSE IF (xx(1) .ge. xx(2)) THEN
         iflag = -2 ! Line sigments must be specified in increasing order of x
         RETURN
      END IF

! If x is outside of array, linearly interprolate from the last two points.
      IF (x .lt. xx(1)) THEN
         ilow = 1; ihigh = 2
      ELSE IF (x .gt. xx(n)) THEN
         ilow = n - 1; ihigh = n
      ELSE
         CALL get_indices(x, xx, 1, n, ilow, ihigh)
      END IF

      y = y_value(x, yy, xx, ilow, ihigh)

      END SUBROUTINE

!*******************************************************************************
!  Find the linearly integrated value of y at position x from line sigments
!  specified by yy(xx).
!*******************************************************************************
      SUBROUTINE line_seg_int(x, y, xx, yy, n, iflag)
      IMPLICIT none

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: n
      INTEGER, INTENT(out) :: iflag

      REAL(rprec), INTENT(out) :: y
      REAL(rprec), INTENT(in)  :: x
      REAL(rprec), DIMENSION(n), INTENT(in)  :: xx, yy
!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER :: ilow, ihigh, i0low, i0high, i

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Check for valid line segments.
      IF (n .le. 1) THEN
         iflag = -1 ! Line sigments require at least two points
         RETURN
      ELSE IF (xx(1) .ge. xx(2)) THEN
         iflag = -2 ! Line sigments must be specified in increasing order of x
         RETURN
      END IF

! If x is outside of array, linearly interprolate from the last two points.
      IF (x .lt. xx(1)) THEN
         ilow = 1; ihigh = 2
      ELSE IF (x .gt. xx(n)) THEN
         ilow = n - 1; ihigh = n
      ELSE
         CALL get_indices(x, xx, 1, n, ilow, ihigh)
      END IF

! Determine where zero is.
      IF (zero .lt. xx(1)) THEN
         i0low = 1; i0high = 2
      ELSE IF (zero .gt. xx(n)) THEN
         i0low = n - 1; i0high = n
      ELSE
         CALL get_indices(zero, xx, 1, n, i0low, i0high)
      END IF

! If 0 and x have the same low and high indices, then the integrate directly
! from 0 to x. Other wise Integrate from:
!   0 to xx(i0high); xx(i0high) to xx(ilow); xx(ilow) to x
      IF (ilow .eq. i0low) THEN
         y = y_value_int(zero, x, yy, xx, ilow, ihigh)
      ELSE
         y = y_value_int(zero, xx(i0high), yy, xx, i0low, i0high)
         DO i = i0high, ilow - 1
             y = y + y_value_int(xx(i), xx(i+1), yy, xx, i, i+1)
         END DO
         y = y + y_value_int(xx(ilow), x, yy, xx, ilow, ihigh)
      END IF

      END SUBROUTINE

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Find bounding indices of the xx array
!*******************************************************************************
      RECURSIVE                                                                &
     &SUBROUTINE get_indices(x, xx, lBound, uBound, ilow, ihigh)
      IMPLICIT none

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), INTENT(in)  :: x
      REAL(rprec), DIMENSION(:), INTENT(in)  :: xx
      INTEGER, INTENT(in) :: lBound, uBound ! Lower and upper bounds of the xx array
      INTEGER, INTENT(out) :: ilow, ihigh
!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER :: hBound ! Half way index of the xx array.

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Perform a tree search to find the interval that x falls between. This
! algorthim works by dividing the array into two sub arrays then determining if
! x is in the lower or upper array. Then the algortrim recurses on the sub array
! until n equals 2. This assumes that the array xx is correctly sorted.
      IF (uBound - lBound .eq. 1) THEN
         ilow = lBound; ihigh = uBound
         RETURN
      END IF

      hBound = (uBound + lBound)/2
      IF (x .ge. xx(lBound) .and. x .lt. xx(hBound)) THEN
         CALL get_indices(x, xx, lBound, hBound, ilow, ihigh)
      ELSE
         CALL get_indices(x, xx, hBound, uBound, ilow, ihigh)
      END IF
      END SUBROUTINE

!*******************************************************************************
!  Find m of y = m*x + b
!*******************************************************************************
      FUNCTION slope(yy, xx, ilow, ihigh)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), INTENT(in) :: xx, yy
      INTEGER, INTENT(in) :: ilow, ihigh
      REAL(rprec) :: slope
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      slope = (yy(ihigh) - yy(ilow))/(xx(ihigh) - xx(ilow))

      END FUNCTION

!*******************************************************************************
!  Find b of y = m*x + b
!*******************************************************************************
      FUNCTION offset(yy, xx, ilow, ihigh)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), INTENT(in) :: xx, yy
      INTEGER, INTENT(in) :: ilow, ihigh
      REAL(rprec) :: offset
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      offset = (xx(ihigh)*yy(ilow) - xx(ilow)*yy(ihigh))                       &
     &       / (xx(ihigh) - xx(ilow))

      END FUNCTION

!*******************************************************************************
!  Find y(x)
!*******************************************************************************
      FUNCTION y_value(x, yy, xx, ilow, ihigh)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), INTENT(in) :: xx, yy
      REAL(rprec), INTENT(in) :: x
      INTEGER, INTENT(in) :: ilow, ihigh
      REAL(rprec) :: y_value
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      y_value = slope(yy, xx, ilow, ihigh)*x                                   &
     &        + offset(yy, xx, ilow, ihigh)

      END FUNCTION

!*******************************************************************************
!  Find Int[y(x),x]
!*******************************************************************************
      FUNCTION y_value_int(x0, x1, yy, xx, ilow, ihigh)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), INTENT(in) :: xx, yy
      REAL(rprec), INTENT(in) :: x0, x1
      INTEGER, INTENT(in) :: ilow, ihigh
      REAL(rprec) :: y_value_int
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      y_value_int = slope(yy, xx, ilow, ihigh)/2.0*(x1**2.0 - x0**2.0)         &
     &            + offset(yy, xx, ilow, ihigh)*(x1 - x0)

      END FUNCTION

!*******************************************************************************
!  Unit Test Functions
!*******************************************************************************
      FUNCTION line_seg_test()
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: line_seg_test
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! This is currently just a stub testing code will be placed in here. If all
! tests pass, return true.
      line_seg_test = .true.
      END FUNCTION

      END MODULE line_segment
