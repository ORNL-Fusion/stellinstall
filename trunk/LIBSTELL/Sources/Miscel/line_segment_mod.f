!*******************************************************************************
!  File line_segment.f
!
!  Module is part of the libstell.a. This module containes code to create
!  profile constructed of line sigments.
!
!*******************************************************************************

      MODULE line_segment
      USE stel_kinds

      PUBLIC line_seg, line_seg_int, line_seg_test
      PRIVATE get_indices

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

      IF (n .le. 1) THEN
         iflag = -1 ! Line sigments require at least two points
         RETURN
      END IF

! If x is outside of array, linearly interprolate from the last two points.
      IF (x .lt. MIN(xx(1), xx(n))) THEN
         ilow = 1; ihigh = 2
      ELSE IF (x .gt. MAX(xx(1), xx(n))) THEN
         ilow = n - 1; ihigh = n
      ELSE
         CALL get_indices(x, xx, 1, n, ilow, ihigh)
      END IF

      y = (yy(ihigh) - yy(ilow))/(xx(ihigh) - xx(ilow))*x                      &
     &  + (xx(ihigh)*yy(ilow) - xx(ilow)*yy(ihigh))                            &
     &  /(xx(ihigh) - xx(ilow))
      END SUBROUTINE

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Fine bounding indices of the xx array
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
      IF (abs(uBound - lBound) .eq. 1) THEN
         ilow = lBound; ihigh = uBound
         RETURN
      END IF

      hBound = (uBound + lBound)/2
      IF (x .ge. MIN(xx(lBound), xx(hBound)) .and.                             &
     &    x .lt. MAX(xx(lBound), xx(hBound))) THEN
         CALL get_indices(x, xx, lBound, hBound, ilow, ihigh)
      ELSE
         CALL get_indices(x, xx, hBound, uBound, ilow, ihigh)
      END IF
      END SUBROUTINE

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
      line_seg_test = .true.
      END FUNCTION
      END MODULE line_segment