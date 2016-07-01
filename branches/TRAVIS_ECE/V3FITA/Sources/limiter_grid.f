!*******************************************************************************
!>  @file limiter_grid.f
!>  @brief Contains module @ref limiter
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref limiter_grid_class.
!>  @par Super Class:
!>  @ref geometric
!*******************************************************************************

      MODULE limiter_grid
      USE stel_kinds
      USE mpi_inc
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  limiter_grid module parameters
!*******************************************************************************
!>  NETCDF r grid start.
      CHARACTER (len=*), PARAMETER :: nc_r0         = 'r0'
!>  NETCDF r grid spacing.
      CHARACTER (len=*), PARAMETER :: nc_dr         = 'dr'
!>  NETCDF z grid start.
      CHARACTER (len=*), PARAMETER :: nc_z0         = 'z0'
!>  NETCDF z grid spacing.
      CHARACTER (len=*), PARAMETER :: nc_dz         = 'dz'
!>  NETCDF grid phi angles.
      CHARACTER (len=*), PARAMETER :: nc_phi_angles = 'phi_angles'
!>  NETCDF grid iso surfaces.
      CHARACTER (len=*), PARAMETER :: nc_iso_grids  = 'iso_grids'

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) limiter_grid base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a limiter signal.
!>  @par Super class:
!>  @ref geometric
!-------------------------------------------------------------------------------
      TYPE limiter_grid_class
!>  Limiter grid starting r.
         REAL (rprec)                            :: r0
!>  Limiter grid r grid spacing.
         REAL (rprec)                            :: dr
!>  Limiter grid starting z.
         REAL (rprec)                            :: z0
!>  Limiter grid z grid spacing.
         REAL (rprec)                            :: dz
!>  Limiter phi planes the limiters are defined at.
         REAL (rprec), DIMENSION(:), POINTER     :: phi => null()
!>  Limiter iso grid.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: grid => null()
      END TYPE limiter_grid_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a limiter_grid_class.
!>
!>  Allocates memory and initializes a @ref limiter_grid_class object from a
!>  netcdf file.
!>
!>  @param[in] lgrid_file
!>  @returns A pointer to a constructed @ref limiter_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION limiter_grid_construct(lgrid_file)
      USE ezcdf
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_grid_class), POINTER :: limiter_grid_construct
      CHARACTER (len=*), INTENT(in)      :: lgrid_file

!  local variables
      INTEGER                            :: lgrid_iou
      INTEGER                            :: status
      INTEGER, DIMENSION(3)              :: dim_lengths
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(limiter_grid_construct)

      CALL cdf_open(lgrid_iou, TRIM(lgrid_file), 'r', status)
      CALL assert_eq(0, status, 'limiter_grid_construct: ' //                  &
     &               'failed to open ', TRIM(lgrid_file))

      CALL cdf_read(lgrid_iou, nc_r0, limiter_grid_construct%r0)
      CALL cdf_read(lgrid_iou, nc_dr, limiter_grid_construct%dr)
      CALL cdf_read(lgrid_iou, nc_z0, limiter_grid_construct%z0)
      CALL cdf_read(lgrid_iou, nc_dz, limiter_grid_construct%dz)

      CALL cdf_inquire(lgrid_iou, nc_phi_angles, dim_lengths)
      ALLOCATE(limiter_grid_construct%phi(dim_lengths(1)))
      CALL cdf_read(lgrid_iou, nc_phi_angles,                                  &
     &              limiter_grid_construct%phi)

      CALL cdf_inquire(lgrid_iou, nc_iso_grids, dim_lengths)
      ALLOCATE(limiter_grid_construct%grid(dim_lengths(1),                     &
     &                                     dim_lengths(2),                     &
     &                                     dim_lengths(3)))
      CALL cdf_read(lgrid_iou, nc_iso_grids,                                   &
     &              limiter_grid_construct%grid)

      CALL cdf_close(lgrid_iou)

      CALL profiler_set_stop_time('limiter_grid_construct', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref limiter_grid_class object.
!>
!>  Deallocates memory and uninitializes a @ref limiter_grid_class object.
!>
!>  @param[inout] this A @ref limiter_grid_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_grid_destruct(this)

!  Declare Arguments
      TYPE (limiter_grid_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%phi)) THEN
         DEALLOCATE(this%phi)
         this%phi => null()
      END IF

      IF (ASSOCIATED(this%grid)) THEN
         DEALLOCATE(this%grid)
         this%grid => null()
      END IF

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Find the number of radial points in the grid.
!>
!>  Find the number of radial points in the grid by measuing the size of the
!>  appropiate dimension.
!>
!>  @param[in] this A @ref limiter_grid_class instance.
!>  @returns Number of radial grid points.
!-------------------------------------------------------------------------------
       FUNCTION limiter_grid_get_num_r(this)

       IMPLICIT NONE

!  Declare Arguments
       INTEGER                               :: limiter_grid_get_num_r
       TYPE (limiter_grid_class), INTENT(in) :: this

!  local variables
       REAL (rprec)                          :: start_time

!  Start of executable code
       start_time = profiler_get_start_time()

       limiter_grid_get_num_r = SIZE(this%grid, 1)

       CALL profiler_set_stop_time('limiter_grid_get_num_r', start_time)

       END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Find the number of vertical points in the grid.
!>
!>  Find the number of vertical points in the grid by measuing the size of the
!>  appropiate dimension.
!>
!>  @param[in] this A @ref limiter_grid_class instance.
!>  @returns Number of radial grid points.
!-------------------------------------------------------------------------------
       FUNCTION limiter_grid_get_num_z(this)

       IMPLICIT NONE

!  Declare Arguments
       INTEGER                               :: limiter_grid_get_num_z
       TYPE (limiter_grid_class), INTENT(in) :: this

!  local variables
       REAL (rprec)                          :: start_time

!  Start of executable code
       start_time = profiler_get_start_time()

       limiter_grid_get_num_z = SIZE(this%grid, 2)

       CALL profiler_set_stop_time('limiter_grid_get_num_z', start_time)

       END FUNCTION

       END MODULE
