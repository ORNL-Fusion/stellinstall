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

!>  B field in the R direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: b_r => null()
!>  B field in the Phi direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: b_p => null()
!>  B field in the Z direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: b_z => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref unprimed_grid_class object.
!>
!>  Allocates memory and initializes a @ref unprimed_grid_class object with an
!>  This computes the vector potential, magnetic field and positions on the
!>  unprimed grid. Magnetic field components are provided on a truncated grid.
!>  Due to the central differencing, field componetents cannot be provided for
!>  first and last points on the r and z grids.
!>
!>  @param[in] mgrid A @ref m_grid_class object.
!>  @param[in] pgrid A @ref primed_grid_class object.
!>  @returns A pointer to a constructed @ref unprimed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION unprimed_grid_construct(mgrid, pgrid)
      USE primed_grid
      USE m_grid
!$    USE omp_lib

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
      REAL (rprec)                                :: v
      REAL (rprec)                                :: ax
      REAL (rprec)                                :: ay
      REAL (rprec)                                :: ar_p
      REAL (rprec)                                :: ar_m
      REAL (rprec)                                :: ap_p
      REAL (rprec)                                :: ap_m
      REAL (rprec)                                :: az_p
      REAL (rprec)                                :: az_m
      REAL (rprec)                                :: r_p
      REAL (rprec)                                :: r_m
      INTEGER                                     :: k_p
      INTEGER                                     :: k_m
      REAL (rprec)                                :: cosv
      REAL (rprec)                                :: sinv
      REAL (rprec)                                :: done
      REAL (rprec)                                :: total
      REAL (rprec)                                :: current
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: rp
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: zp

!  local parameters
      CHARACTER (len=*), PARAMETER                :: clear_screen =            &
     &   achar(27)//"[2K"//achar(13)

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(unprimed_grid_construct)

      num_r = SIZE(mgrid%a_r, 1)
      num_p = SIZE(mgrid%a_p, 3)
      num_z = SIZE(mgrid%a_z, 2)

      ALLOCATE(unprimed_grid_construct%a_r(num_r, num_z, num_p))
      ALLOCATE(unprimed_grid_construct%a_p(num_r, num_z, num_p))
      ALLOCATE(unprimed_grid_construct%a_z(num_r, num_z, num_p))

      ALLOCATE(unprimed_grid_construct%b_r(num_r, num_z, num_p))
      ALLOCATE(unprimed_grid_construct%b_p(num_r, num_z, num_p))
      ALLOCATE(unprimed_grid_construct%b_z(num_r, num_z, num_p))

      total = num_p*num_z*num_r
!$    total = total/OMP_GET_MAX_THREADS()
      current = 0.0

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i, j, k, v, r, zp, x, y, ax, ay, rp, cosv, sinv, ar_p,
!$OMP&         ar_m, ap_p, ap_m, az_p, az_m, r_p, r_m, k_p, k_m, done,
!$OMP&         current)

      ALLOCATE(rp(SIZE(pgrid%x, 1),                                            &
     &            SIZE(pgrid%x, 2),                                            &
     &            SIZE(pgrid%x, 3)))

      ALLOCATE(zp(SIZE(pgrid%x, 1),                                            &
     &            SIZE(pgrid%x, 2),                                            &
     &            SIZE(pgrid%x, 3)))

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO k = 1, num_p
         v = (k - 1.0)*pgrid%dv

         cosv = COS(v)
         sinv = SIN(v)

         DO j = 1, num_z
            zp = (pgrid%z - ((j - 1.0)*mgrid%dz + mgrid%zmin))**2.0

            DO i = 1, num_r

!$             IF (OMP_GET_THREAD_NUM() .eq. 0) THEN
                  current = current + 1.0
                  done = 100*current/total
                  SELECT CASE (MOD(INT(current),4))
                     CASE (0)
                        WRITE (*,1000,ADVANCE='NO') clear_screen, done
                     CASE (1)
                        WRITE (*,1001,ADVANCE='NO') clear_screen, done
                     CASE (2)
                        WRITE (*,1002,ADVANCE='NO') clear_screen, done
                     CASE (3)
                        WRITE (*,1003,ADVANCE='NO') clear_screen, done
                  END SELECT
!$             END IF

               r = (i - 1.0)*mgrid%dr + mgrid%rmin

               x = r*cosv
               y = r*sinv

               rp = SQRT((pgrid%x - x)**2.0 +                                  &
     &                   (pgrid%y - y)**2.0 +                                  &
     &                   zp)

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
      DEALLOCATE(zp)

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO k = 1, num_p
         IF (k .eq. 1) THEN
            k_p = k + 1
            k_m = num_p
         ELSE IF (k .eq. num_p) THEN
            k_p = 1
            k_m = k - 1
         ELSE
            k_p = k + 1
            k_m = k - 1
         END IF

         DO j = 2, num_z - 1

            DO i = 2, num_r - 1
               r_p = i*mgrid%dr + mgrid%rmin
               r_m = (i - 2.0)*mgrid%dr + mgrid%rmin
               r = (i - 1.0)*mgrid%dr + mgrid%rmin

               ar_p = unprimed_grid_construct%a_r(i + 1,j,k)
               ar_m = unprimed_grid_construct%a_r(i - 1,j,k)

               az_p = unprimed_grid_construct%a_z(i,j + 1,k)
               az_m = unprimed_grid_construct%a_z(i,j - 1,k)

               ap_p = unprimed_grid_construct%a_p(i,j,k_p)
               ap_m = unprimed_grid_construct%a_p(i,j,k_m)

               unprimed_grid_construct%b_r(i,j,k) =                            & ! 1/rdazdp - dapdz
     &            (az_p - az_m)/(2.0*pgrid%dv*r) -                             &
     &            (ap_p - ap_m)/(2.0*mgrid%dz)

               unprimed_grid_construct%b_p(i,j,k) =                            & ! dardz - dazdr
     &            (ar_p - ar_m)/(2.0*mgrid%dz) -                               &
     &            (az_p - az_m)/(2.0*mgrid%dr)

               unprimed_grid_construct%b_z(i,j,k) =                            & ! 1/r(drapdr - dardp)
     &            (r_p*ap_p - r_m*ap_m)/(2.0*mgrid%dr*r) -                     &
     &            (ar_p - ar_m)/(2.0*pgrid%dv*r)

            END DO
         END DO
      END DO
!$OMP END DO

!$OMP WORKSHARE
      unprimed_grid_construct%b_r(1,:,:) = 0.0
      unprimed_grid_construct%b_p(1,:,:) = 0.0
      unprimed_grid_construct%b_z(1,:,:) = 0.0
      unprimed_grid_construct%b_r(num_r,:,:) = 0.0
      unprimed_grid_construct%b_p(num_r,:,:) = 0.0
      unprimed_grid_construct%b_z(num_r,:,:) = 0.0

      unprimed_grid_construct%b_r(:,1,:) = 0.0
      unprimed_grid_construct%b_p(:,1,:) = 0.0
      unprimed_grid_construct%b_z(:,1,:) = 0.0
      unprimed_grid_construct%b_r(:,num_z,:) = 0.0
      unprimed_grid_construct%b_p(:,num_z,:) = 0.0
      unprimed_grid_construct%b_z(:,num_z,:) = 0.0
!$OMP END WORKSHARE

!$OMP END PARALLEL

      WRITE (*,1004) clear_screen

      CALL profiler_set_stop_time('unprimed_grid_construct', start_time)

1000  FORMAT(a,' | ',f6.2,' % Finished')
1001  FORMAT(a,' / ',f6.2,' % Finished')
1002  FORMAT(a,' - ',f6.2,' % Finished')
1003  FORMAT(a,' \ ',f6.2,' % Finished')
1004  FORMAT(a,'Unprimed Grid Finished')

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

      IF (ASSOCIATED(this%b_r)) THEN
         DEALLOCATE(this%b_r)
         this%b_r => null()
      END IF

      IF (ASSOCIATED(this%b_p)) THEN
         DEALLOCATE(this%b_p)
         this%b_p => null()
      END IF

      IF (ASSOCIATED(this%b_z)) THEN
         DEALLOCATE(this%b_z)
         this%b_z => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

      END MODULE
