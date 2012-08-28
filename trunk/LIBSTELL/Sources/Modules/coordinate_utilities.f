!*******************************************************************************
!  File coordinate_utilities.f
!
!  Module is part of the LIBSTELL. This modules containes code to convert from
!  different coordinate systems.
!
!*******************************************************************************

      MODULE coordinate_utilities
      USE stel_kinds
      USE stel_constants

      PUBLIC :: cart_to_cyl, cood_utils_test
      PRIVATE :: check

      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Convert from carteasan coordinates to cylindical coordinates.
!*******************************************************************************
      FUNCTION cart_to_cyl(cart)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(3), INTENT(in) :: cart
      REAL(rprec), DIMENSION(3)             :: cart_to_cyl

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      cart_to_cyl(1) = SQRT(cart(1)**2 + cart(2)**2)
      cart_to_cyl(2) = ATAN2(cart(2), cart(1))
      cart_to_cyl(3) = cart(3)

      END FUNCTION

!*******************************************************************************
!  Unit Test Functions
!*******************************************************************************
!*******************************************************************************
!  Main test function
!*******************************************************************************
      FUNCTION cood_utils_test()
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: cood_utils_test
      REAL(rprec), DIMENSION(3) :: result

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Test cart_to_cyl function. cart(0,0,1) = cyl(0,0,1)
      result = cart_to_cyl((/ 0.0d+0, 0.0d+0, 1.0d+0 /))
      cood_utils_test = check(0.0d+0, result(1), 1, "cart_to_cyl")
      cood_utils_test = check(0.0d+0, result(2), 2, "cart_to_cyl")
      cood_utils_test = check(1.0d+0, result(3), 3, "cart_to_cyl")

! Test cart_to_cyl function. cart(-1,0,0) = cyl(1,0,1)
      result = cart_to_cyl((/ -1.0d+0, 0.0d+0, 0.0d+0 /))
      cood_utils_test = check(1.0d+0, result(1), 4, "cart_to_cyl")
      cood_utils_test = check(pi, result(2), 5, "cart_to_cyl")
      cood_utils_test = check(0.0d+0, result(3), 6, "cart_to_cyl")

! Test cart_to_cyl function. cart(-1,0,-1) = cyl(1,0,1)
      result = cart_to_cyl((/ 0.0d+0, 1.0d+0, -1.0d+0 /))
      cood_utils_test = check(1.0d+0, result(1), 7, "cart_to_cyl")
      cood_utils_test = check(pi/2.0d+0, result(2), 8, "cart_to_cyl")
      cood_utils_test = check(-1.0d+0, result(3), 9, "cart_to_cyl")

      END FUNCTION

!*******************************************************************************
!  Check Test result
!*******************************************************************************
      FUNCTION check(expected, recieved, testNum, name)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: check
      REAL(rprec), INTENT(in) :: expected, recieved
      INTEGER, INTENT(in) :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      check = expected .eq. recieved
      IF (.not.check) THEN
         write(*,*) "coordinate_utilities.f: ", name, " test", testNum,        &
     &              "failed."
         write(*,*) "Expected", expected, "Recieved", recieved
      END IF

      END FUNCTION

      END MODULE
