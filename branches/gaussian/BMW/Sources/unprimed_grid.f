!*******************************************************************************
!>  @file unprimed_grid.f
!>  @brief Contains module @ref unprimed_grid.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref primed_grid_class. This contains the
!>  state variables to define the currents and positions of the volumn integral.
!*******************************************************************************
      MODULE unprimed_grid
      USE stel_kinds, ONLY: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) unprimed grid base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a unprimed grid. This is grid the volume integral
!>  will be summed over.
!-------------------------------------------------------------------------------
      TYPE unprimed_grid_class
!>  Vector potential in the R direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_r => null()
!>  Vector potential in the Phi direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_p => null()
!>  Vector potential in the Z direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_z => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref unprimed_grid_class object.
!>
!>  Allocates memory and initializes a @ref unprimed_grid_class object with an
!>  This computes the vector potential and positions on the unprimed grid.
!>
!>  @param[in] mgrid A @ref m_grid_class object.
!>  @param[in] pgrid A @ref primed_grid_class object.
!>  @returns A pointer to a constructed @ref unprimed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION unprimed_grid_construct(mgrid, pgrid)
      USE primed_grid
      USE m_grid

      IMPLICIT NONE

!  Declare Arguments
      TYPE (unprimed_grid_class), POINTER :: unprimed_grid_construct
      TYPE (m_grid_class), INTENT(in)             :: mgrid
      TYPE (primed_grid_class), INTENT(in)        :: pgrid

!  local variables
      REAL (rprec)                                :: start_time
      INTEGER                                     :: i
      INTEGER                                     :: j
      INTEGER                                     :: k
      INTEGER                                     :: num_r
      INTEGER                                     :: num_p
      INTEGER                                     :: num_z
      REAL (rprec)                                :: x
      REAL (rprec)                                :: y
      REAL (rprec)                                :: r
      REAL (rprec)                                :: z
      REAL (rprec)                                :: v
      REAL (rprec)                                :: ax
      REAL (rprec)                                :: ay
      REAL (rprec)                                :: cosv
      REAL (rprec)                                :: sinv
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: rp

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(unprimed_grid_construct)

      num_r = SIZE(mgrid%a_r, 1)
      num_p = SIZE(mgrid%a_p, 3)
      num_z = SIZE(mgrid%a_z, 2)

      ALLOCATE(unprimed_grid_construct%a_r(num_r, num_z, num_p))
      ALLOCATE(unprimed_grid_construct%a_p(num_r, num_z, num_p))
      ALLOCATE(unprimed_grid_construct%a_z(num_r, num_z, num_p))

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i, j, k, v, r, z, x, y, ax, ay, rp, cosv, sinv)

      ALLOCATE(rp(SIZE(pgrid%x, 1),                                            &
     &            SIZE(pgrid%x, 2),                                            &
     &            SIZE(pgrid%x, 3)))

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO k = 1, num_p
         v = (k - 1.0)*pgrid%dv

         WRITE (*,1000) v

         cosv = COS(v)
         sinv = SIN(v)

         DO j = 1, num_z
            z = (j - 1.0)*mgrid%dz + mgrid%zmin

            DO i = 1, num_r
               r = (i - 1.0)*mgrid%dr + mgrid%rmin

               x = r*cosv
               y = r*sinv

               rp = SQRT((pgrid%x - x)**2.0 +                                  &
     &                   (pgrid%y - y)**2.0 +                                  &
     &                   (pgrid%z - z)**2.0)

               ax = SUM(pgrid%j_x/rp)*pgrid%dvol
               ay = SUM(pgrid%j_y/rp)*pgrid%dvol

               unprimed_grid_construct%a_r(i,j,k) =  x/r*ax + y/r*ay           &
     &                                            + mgrid%a_r(i,j,k)
               unprimed_grid_construct%a_p(i,j,k) = -y/r*ax + x/r*ay           &
     &                                            + mgrid%a_p(i,j,k)
               unprimed_grid_construct%a_z(i,j,k) = SUM(pgrid%j_z/rp)          &
     &                                            * pgrid%dvol                 &
     &                                            + mgrid%a_z(i,j,k)
            END DO
         END DO
      END DO
!$OMP END DO

      DEALLOCATE(rp)

!$OMP END PARALLEL

      WRITE (*,1001)

      CALL profiler_set_stop_time('unprimed_grid_construct', start_time)

1000  FORMAT('Phi = ',f4.2)
1001  FORMAT('Unprimed Grid Finished')

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref primed_grid_class object.
!>
!>  Deallocates memory and uninitializes a @ref primed_grid_class object.
!>
!>  @param[inout] this A @ref unprimed_grid_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE unprimed_grid_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (unprimed_grid_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%a_r)) THEN
         DEALLOCATE(this%a_r)
         this%a_r => null()
      END IF

      IF (ASSOCIATED(this%a_p)) THEN
         DEALLOCATE(this%a_p)
         this%a_p => null()
      END IF

      IF (ASSOCIATED(this%a_z)) THEN
         DEALLOCATE(this%a_z)
         this%a_z => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

      END MODULE
