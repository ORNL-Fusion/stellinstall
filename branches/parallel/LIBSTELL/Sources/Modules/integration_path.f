!*******************************************************************************
!>  @file magnetic.f
!>  @brief Contains module @ref integration_path
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module is part of the LIBSTELL. This modules containes code to define and
!>  integrate along an arbitray path.
!*******************************************************************************

      MODULE integration_path
      USE stel_kinds
#if defined (MPI_OPT)
      USE mpi
#endif

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) vertex
!  2) path
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  A single point in space defined by an z, y, z coordinate. A vertex is
!>  structured as a singly linked list.
!-------------------------------------------------------------------------------
      TYPE vertex
         REAL (rprec)           :: x, y, z
         TYPE (vertex), POINTER :: next => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref vertex types using
!>  @ref path_construct_new or @ref path_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE path_construct
         MODULE PROCEDURE path_construct_new,                                  &
     &                    path_construct_mpi
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for checking the results of the unit tests
!-------------------------------------------------------------------------------
      INTERFACE check
         MODULE PROCEDURE check_log, check_real, check_int
      END INTERFACE

      PUBLIC :: path_construct, path_append_vertex, path_destruct,             &
     &          path_test, path_integrate
      PRIVATE :: check, check_real, check_log, check_int,                      &
     &           integrate, test_function

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a single @ref vertex.
!>
!>  Allocates memory and initializes a @ref vertex object.
!>
!>  @param[in] x X position of the vertex object.
!>  @param[in] y Y position of the vertex object.
!>  @param[in] z Z position of the vertex object.
!>  @returns A pointer to a constructed @ref vertex object.
!-------------------------------------------------------------------------------
      FUNCTION path_construct_new(x, y, z)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER   :: path_construct_new
      REAL (rprec), INTENT(in) :: x
      REAL (rprec), INTENT(in) :: y
      REAL (rprec), INTENT(in) :: z

!  Start of executable code
      ALLOCATE(path_construct_new)

      path_construct_new%x = x
      path_construct_new%y = y
      path_construct_new%z = z

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref vertex object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref vertex object. This gets called
!>  recursively until the last object is reached.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref vertex object.
!-------------------------------------------------------------------------------
      RECURSIVE FUNCTION path_construct_mpi(intra_comm)                        &
     &   RESULT(new_vertex)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER :: new_vertex
      INTEGER, INTENT(in)    :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                :: error
      LOGICAL                :: next

!  Start of executable code
      ALLOCATE(new_vertex)

      CALL MPI_BCAST(next, 1, MPI_LOGICAL, 0, intra_comm, error)
      IF (next) THEN
         new_vertex%next => path_construct_mpi(intra_comm)
      END IF

      CALL MPI_BCAST(new_vertex%x, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(new_vertex%y, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(new_vertex%z, 1, MPI_REAL8, 0, intra_comm, error)

#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref vertex object.
!>
!>  Deallocates memory and uninitializes a @ref vertex object. This recursively
!>  deconstructed the next vertex until the last in the linked list is found.
!>
!>  @param[inout] this A @ref vertex instance.
!-------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE path_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER :: this

!  Start of executable code.
      IF (ASSOCIATED(this%next)) THEN
         CALL path_destruct(this%next)
         this%next => null()
      END IF

      DEALLOCATE(this)
      this => null()

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Append a @ref vertex to a path.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, a new vertex is allocated and appended to the path.
!>  This allow works as a constructer and allocates the first vertex if needed.
!>
!>  @param[inout] this Vertex to append path to.
!>  @param[in]    x    X position of the vertex object.
!>  @param[in]    y    Y position of the vertex object.
!>  @param[in]    z    Z position of the vertex object.
!-------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE path_append_vertex(this, x, y, z)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER   :: this
      REAL (rprec), INTENT(in) :: x
      REAL (rprec), INTENT(in) :: y
      REAL (rprec), INTENT(in) :: z

!  Start of executable code
      IF (ASSOCIATED(this)) THEN
         IF (ASSOCIATED(this%next)) THEN
            CALL path_append_vertex(this%next, x, y, z)
         ELSE
            this%next => path_construct(x, y, z)
         END IF
      ELSE
         this => path_construct(x, y, z)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Integrate along the path.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, a new vertex is allocated and appended to the path.
!>  The integrand is proveded by means of call back function.
!>
!>  @param[inout] this                 Starting vertex to integrate to the end.
!>  @param        integration_function Function pointer that defines the
!>                                     integrand.
!>  @param[in]    context              Generic object that contains data for the
!>                                     integration function.
!>  The total integrated path to the end.
!-------------------------------------------------------------------------------
      RECURSIVE FUNCTION path_integrate(this, integration_function,            &
     &                                  context) RESULT(total)
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: total
      TYPE (vertex), INTENT(inout)  :: this
      CHARACTER (len=1), INTENT(in) :: context(:)
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart, dx)
         USE stel_kinds
         REAL (rprec) :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
         REAL (rprec), INTENT(in)               :: dx
         END FUNCTION
      END INTERFACE

!  Start of executable code
      If (ASSOCIATED(this%next)) THEN
         total = path_integrate(this%next, integration_function,               &
     &                          context)
         total = total + integrate(context, this, this%next,                   &
     &                             integration_function)
      ELSE
         total = 0.0
      END IF

      END FUNCTION

!*******************************************************************************
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize data to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref vertex instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE path_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), INTENT(inout) :: this
      INTEGER, INTENT(in)          :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                      :: error
      LOGICAL                      :: next

!  Start of executable code
      next = ASSOCIATED(this%next)
      CALL MPI_BCAST(next, 1, MPI_LOGICAL, 0, intra_comm, error)
      IF (next) THEN
         CALL path_sync_children(this%next, intra_comm)
      END IF

      CALL MPI_BCAST(this%x, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%y, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%z, 1, MPI_REAL8, 0, intra_comm, error)

#endif

      END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Line integrate between to points.
!>
!>  This divids the straight line path defined by two vertices and line
!>  integrates the function. The integrand is proveded by means of call back
!>  function.
!>
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integrate(context, vertex1, vertex2,                            &
     &                   integration_function)
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                 :: integrate
      CHARACTER(len=1), INTENT(in) :: context(:)
      TYPE (vertex), INTENT(in)    :: vertex1
      TYPE (vertex), INTENT(in)    :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart, dx)
         USE stel_kinds
         REAL (rprec) :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
         REAL (rprec), INTENT(in)               :: dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3) :: xcart
      REAL (rprec), DIMENSION(3) :: dxcart
      REAL (rprec)               :: dx
      INTEGER                    :: i, nsteps

!  local parameters
      REAL(rprec), PARAMETER :: dxOptimal = 0.0025

!  Start of executable code
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
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Path unit test function.
!>
!>  This runs the associated unit tests and returns the result.
!>
!>  @returns True if the tests pass and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION path_test()

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: path_test

!  local variables
      REAL (rprec)                  :: result
      TYPE (vertex), POINTER        :: test_path => null()
      CHARACTER(len=1), ALLOCATABLE :: context(:)

!  Start of executable code
! Test to make sure the vertices begin in an unallocated state.
      path_test = check(.false., ASSOCIATED(test_path), 1, "ASSOCIATED")
      IF (.not.path_test) RETURN

! Test to make sure first vertex is created.
      CALL path_append_vertex(test_path,                                       &
     &                        1.0_rprec, 2.0_rprec, 3.0_rprec)
      path_test = check(.true., ASSOCIATED(test_path), 1,                      &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(1.0_rprec, test_path%x, 2, "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(2.0_rprec, test_path%y, 3, "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(3.0_rprec, test_path%z, 4, "path_append_vertex")
      IF (.not.path_test) RETURN

! Test to make sure second vertex is appended.
      CALL path_append_vertex(test_path,                                       &
     &                        4.0_rprec, 5.0_rprec, 6.0_rprec)
      path_test = check(.true., ASSOCIATED(test_path%next), 5,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(1.0_rprec, test_path%x, 6, "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(2.0_rprec, test_path%y, 7, "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(3.0_rprec, test_path%z, 8, "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(4.0_rprec, test_path%next%x, 9,                        &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(5.0_rprec, test_path%next%y, 10,                       &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(6.0_rprec, test_path%next%z, 11,                       &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

! Test to make sure third vertex is appended.
      CALL path_append_vertex(test_path,                                       &
     &                        7.0_rprec, 8.0_rprec, 9.0_rprec)
      path_test = check(.true., ASSOCIATED(test_path%next%next),               &
     &                  12, "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(7.0_rprec, test_path%next%next%x, 12,                  &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(8.0_rprec, test_path%next%next%y, 13,                  &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(9.0_rprec, test_path%next%next%z, 14,                  &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

! Test to make sure path object is destroyed.
      CALL path_destruct(test_path)
      path_test = check(.false., ASSOCIATED(test_path), 1,                     &
     &                  "path_destruct")
      IF (.not.path_test) RETURN

! Test path integration.
      CALL path_append_vertex(test_path,                                       &
     &                        1.0_rprec, 0.0_rprec, 0.0_rprec)
      CALL path_append_vertex(test_path,                                       &
     &                        0.0_rprec, 0.0_rprec, 0.0_rprec)
      result = path_integrate(test_path, test_function, context)
      path_test = check(400.0_rprec, result, 1, "path_integrate")
      CALL path_destruct(test_path)

      END FUNCTION

!*******************************************************************************
!  CHECK FUNCTIONS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Check a logical value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_log(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_log
      LOGICAL, INTENT(in)           :: expected
      LOGICAL, INTENT(in)           :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check_log = expected .eqv. received
      IF (.not.check_log) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check a real value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_real(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_real
      REAL(rprec), INTENT(in)       :: expected
      REAL(rprec), INTENT(in)       :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check_real = expected .eq. received
      IF (.not.check_real) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check an integer value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_int(expected, received, testNum, name)
      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_int
      INTEGER, INTENT(in)           :: expected
      INTEGER, INTENT(in)           :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check_int = expected .eq. received
      IF (.not.check_int) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Call back function to test the integration.
!>
!>  This defined a function of f(x) = 1.
!>
!>  @param[in] context Generic object that contains data for the integration
!>                     function.
!>  @param[in] xcart   Current point in the integration.
!>  @param[in] dxcart  Vector direction of the current integration path.
!>  @param[in] dx      Scaler length of the current integration path.
!>  @returns One
!-------------------------------------------------------------------------------
      FUNCTION test_function(context, xcart, dxcart, dx)
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: test_function
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: dx

!  Start of executable code
      test_function = 1.0

      END FUNCTION

      END MODULE
