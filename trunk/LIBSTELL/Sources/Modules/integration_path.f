!*******************************************************************************
!  File integration_path.f
!
!  Module is part of the LIBSTELL. This modules containes code define and
!  integrate along an arbitray path.
!
!*******************************************************************************

      MODULE integration_path
      USE stel_kinds

!*******************************************************************************
!  Object representing a vertex point in cartesian coordinates.
!*******************************************************************************
      TYPE vertex
         REAL(rprec) :: x, y, z
      END TYPE

!*******************************************************************************
!  Object representing a path. A path is defined by two or more vertex points.
!*******************************************************************************
      TYPE path
         TYPE(vertex), DIMENSION(:), POINTER :: vertices
      END TYPE

      INTERFACE check
         MODULE PROCEDURE check_log, check_real, check_int
      END INTERFACE

      PUBLIC :: path_append_vertex, path_destroy, path_test
      PRIVATE :: check, check_real, check_log, check_int

      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Append a vertex to a path.
!
!  If there are no vertexs appended allocate the array pointer for the first
!  vertex. Otherwise resize the array to add a new path.
!
!*******************************************************************************
      SUBROUTINE path_append_vertex(this, x, y, z)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE(path), INTENT(inout) :: this
      REAL(rprec), INTENT(in)   :: x, y, z

!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER :: status
      TYPE(vertex), DIMENSION(:), ALLOCATABLE :: vertices

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      IF (ASSOCIATED(this%vertices)) THEN
         ALLOCATE(vertices(SIZE(this%vertices) + 1), STAT = status)

         vertices(1:UBOUND(this%vertices, 1)) = this%vertices
         DEALLOCATE(this%vertices)
         this%vertices = vertices
      ELSE
         ALLOCATE(this%vertices(1), STAT = status)
 
         this%vertices(1)%x = x
         this%vertices(1)%z = y
         this%vertices(1)%z = z
      ENDIF

      END SUBROUTINE

!*******************************************************************************
!  Destroy a path.
!
!  Deallocate all memory associtaed with a path.
!
!*******************************************************************************
      SUBROUTINE path_destroy(this)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE(path), INTENT(inout) :: this

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      IF (ASSOCIATED(this%vertices)) THEN
         DEALLOCATE(this%vertices)
      ENDIF

      END SUBROUTINE

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************

!*******************************************************************************
!  Unit Test Functions
!*******************************************************************************
!*******************************************************************************
!  Main test function
!*******************************************************************************
      FUNCTION path_test()
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: path_test
      TYPE(path) :: test_path

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Test to make sure the vertices begin in an unallocated state.
      path_test = check(.false., ASSOCIATED(test_path%vertices), 1,            &
     &                  "ASSOCIATED")
      IF (.not.path_test) RETURN

! Test to make sure first vertex is created.
      CALL path_append_vertex(test_path,                                       &
     &                        1.0_rprec, 2.0_rprec, 3.0_rprec)
      path_test = check(.true., ASSOCIATED(test_path%vertices), 1,             &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(1.0_rprec, test_path%vertices(1)%x, 2,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(2.0_rprec, test_path%vertices(1)%y, 3,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(3.0_rprec, test_path%vertices(1)%z, 4,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

! Test to make sure second vertex is appended.
      CALL path_append_vertex(test_path,                                       &
     &                        4.0_rprec, 5.0_rprec, 6.0_rprec)
      path_test = check(.true., ASSOCIATED(test_path%vertices), 5,             &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(2, SIZE(test_path%vertices), 6,                        &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(4.0_rprec, test_path%vertices(2)%x, 7,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(5.0_rprec, test_path%vertices(2)%y, 8,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(6.0_rprec, test_path%vertices(2)%z, 9,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

! Test to make sure path object is destroyed.
      CALL path_destroy(test_path)       
      path_test = check(.false., ASSOCIATED(test_path%vertices), 1,            &
     &                  "path_destroy")
      IF (.not.path_test) RETURN

      END FUNCTION

!*******************************************************************************
!  Check Test result logical
!*******************************************************************************
      FUNCTION check_log(expected, recieved, testNum, name)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: check_log
      LOGICAL, INTENT(in) :: expected, recieved
      INTEGER, INTENT(in) :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      check_log = expected .eqv. recieved
      IF (.not.check_log) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Recieved", recieved
      END IF

      END FUNCTION

!*******************************************************************************
!  Check Test result real
!*******************************************************************************
      FUNCTION check_real(expected, recieved, testNum, name)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: check_real
      REAL(rprec), INTENT(in) :: expected, recieved
      INTEGER, INTENT(in) :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      check_real = expected .eq. recieved
      IF (.not.check_real) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Recieved", recieved
      END IF

      END FUNCTION

!*******************************************************************************
!  Check Test result Integer
!*******************************************************************************
      FUNCTION check_int(expected, recieved, testNum, name)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: check_int
      INTEGER, INTENT(in) :: expected, recieved, testNum
      CHARACTER (LEN=*), INTENT(in) :: name
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      check_int = expected .eq. recieved
      IF (.not.check_int) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Recieved", recieved
      END IF

      END FUNCTION

      END MODULE
