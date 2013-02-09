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
         TYPE(vertex), DIMENSION(:), POINTER :: vertices => NULL()
      END TYPE

      INTERFACE check
         MODULE PROCEDURE check_log, check_real, check_int
      END INTERFACE

      PUBLIC :: path_append_vertex, path_destroy, path_test                    &
     &          path_integrate
      PRIVATE :: check, check_real, check_log, check_int, check_status         &
     &           integrate, test_function

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
! Create a temporary array the same size as current number of vertices. Copy the
! contents of the old array to the new array then delete the old array.
         ALLOCATE(vertices(SIZE(this%vertices)), STAT = status)
         CALL check_status(status)
         vertices = this%vertices
         DEALLOCATE(this%vertices)

! Create a new array with one extra element. Copy the contents of the temporary
! array to the new one the delete the temporary array.
         ALLOCATE(this%vertices(SIZE(vertices) + 1), STAT = status)
         CALL check_status(status)
         this%vertices(1:SIZE(vertices)) = vertices
         DEALLOCATE(vertices)

         this%vertices(SIZE(this%vertices))%x = x
         this%vertices(SIZE(this%vertices))%y = y
         this%vertices(SIZE(this%vertices))%z = z
      ELSE
         ALLOCATE(this%vertices(1), STAT = status)
         CALL check_status(status)

         this%vertices(1)%x = x
         this%vertices(1)%y = y
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
         this%vertices => NULL()
      ENDIF

      END SUBROUTINE

!*******************************************************************************
!  Integrate along the path.
!
!  This integrates the specified function along the path.
!
!*******************************************************************************
      FUNCTION path_integrate(this, integration_function, context)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE(path), INTENT(inout) :: this
      REAL(rprec)               :: path_integrate
      CHARACTER(len=1), INTENT(in) :: context(:)
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart, dx)
         USE stel_kinds
         CHARACTER(len=1), INTENT(in)          :: context(:)
         REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
         REAL(rprec), INTENT(in)               :: dx
         REAL(rprec)               :: integration_function
         END FUNCTION
      END INTERFACE

!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER :: i

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      IF (SIZE(this%vertices) .lt. 2) THEN
         STOP 'Path integrals require a minimum of two vertices.'
      ENDIF

      DO i = 1, SIZE(this%vertices) - 1
         path_integrate = path_integrate +                                     &
     &                    integrate(context,                                   &
     &                              this%vertices(i),                          &
     &                              this%vertices(i + 1),                      &
     &                              integration_function)
      ENDDO

      END FUNCTION

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Integrate along the segment.
!
!  This integrates the specified function along the line specified by two
!  vertices.
!
!*******************************************************************************
      FUNCTION integrate(context, vertex1, vertex2,                            &
     &                   integration_function)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE(vertex), INTENT(in) :: vertex1, vertex2
      REAL(rprec)              :: integrate
      CHARACTER(len=1), INTENT(in) :: context(:)
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart, dx)
         USE stel_kinds
         CHARACTER(len=1), INTENT(in)          :: context(:)
         REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
         REAL(rprec), INTENT(in)               :: dx
         REAL(rprec)               :: integration_function
         END FUNCTION
      END INTERFACE

!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(3) :: xcart, dxcart
      REAL(rprec)               :: dx
      INTEGER                   :: nsteps, i

!-------------------------------------------------------------------------------
!  Local Parameter declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), PARAMETER :: dxOptimal = 0.0025

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Determine the number of integration steps to take by dividing the path length
! by the step length and rounding to the nearest integer.
      dxcart(1) = vertex2%x - vertex1%x
      dxcart(2) = vertex2%y - vertex1%y
      dxcart(3) = vertex2%z - vertex1%z

! Store the entire length in dx temporarly.
      dx = SQRT(DOT_PRODUCT(dxcart, dxcart))
      nsteps = INT(dx/dxOptimal)

! Choose the actual step size.
      dxcart = dxcart/nsteps
      dx = SQRT(DOT_PRODUCT(dxcart, dxcart))

      integrate = 0.0
      xcart = (/ vertex1%x, vertex1%y, vertex1%z /)
      DO i = 1, nsteps
         xcart = xcart + dxcart
         integrate = integrate                                                 &
     &             + integration_function(context, xcart, dxcart, dx)
      ENDDO

      END FUNCTION

!*******************************************************************************
!  Check Allocation Status
!
!  If an array allocation has failed terminate the program.
!
!*******************************************************************************
      SUBROUTINE check_status(status)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: status

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      IF (status .ne. 0) THEN
         STOP 'ALLOCATION error in LIBSTELL integration_path.f'
      ENDIF

      END SUBROUTINE

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
      LOGICAL                       :: path_test
      REAL(rprec)                   :: result
      TYPE(path)                    :: test_path
      CHARACTER(len=1), ALLOCATABLE :: context(:)

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
      path_test = check(1.0_rprec, test_path%vertices(1)%x, 7,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(2.0_rprec, test_path%vertices(1)%y, 8,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(3.0_rprec, test_path%vertices(1)%z, 9,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(4.0_rprec, test_path%vertices(2)%x, 10,                &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(5.0_rprec, test_path%vertices(2)%y, 11,                &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(6.0_rprec, test_path%vertices(2)%z, 12,                &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

! Test to make sure third vertex is appended.
      CALL path_append_vertex(test_path,                                       &
     &                        7.0_rprec, 8.0_rprec, 9.0_rprec)
      path_test = check(.true., ASSOCIATED(test_path%vertices), 13,            &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(3, SIZE(test_path%vertices), 14,                       &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(7.0_rprec, test_path%vertices(3)%x, 15,                &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(8.0_rprec, test_path%vertices(3)%y, 16,                &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(9.0_rprec, test_path%vertices(3)%z, 17,                &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

! Test to make sure path object is destroyed.
      CALL path_destroy(test_path)
      path_test = check(.false., ASSOCIATED(test_path%vertices), 1,            &
     &                  "path_destroy")
      IF (.not.path_test) RETURN

! Test path integration.
      CALL path_append_vertex(test_path,                                       &
     &                        1.0_rprec, 0.0_rprec, 0.0_rprec)
      CALL path_append_vertex(test_path,                                       &
     &                        0.0_rprec, 0.0_rprec, 0.0_rprec)
      result = path_integrate(test_path, test_function, context)
      path_test = check(400.0_rprec, result, 1, "path_integrate")
      CALL path_destroy(test_path)

      END FUNCTION

!*******************************************************************************
!  Check Test result logical
!*******************************************************************************
      FUNCTION check_log(expected, recieved, testNum, name)
      IMPLICIT NONE

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
      IMPLICIT NONE

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
      IMPLICIT NONE

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

!*******************************************************************************
!  Test Integration Function
!*******************************************************************************
      FUNCTION test_function(context, xcart, dxcart, dx)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)               :: test_function

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      test_function = 1.0

      END FUNCTION

      END MODULE
