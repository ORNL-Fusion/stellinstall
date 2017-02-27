!*******************************************************************************
!>  @file m_grid.f
!>  @brief Contains module @ref m_grid.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref m_grid_class. This contains the
!>  state variables to define the vacuum vector potential.
!*******************************************************************************
      MODULE m_grid
      USE stel_kinds, ONLY: rprec
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) primed grid base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a primed grid. This is grid the volume integral will
!>  be summed over.
!-------------------------------------------------------------------------------
      TYPE m_grid_class
!>  Minimum R position.
         REAL (rprec)                            :: rmin
!>  Maximum R position.
         REAL (rprec)                            :: rmax
!>  Minimum Z position.
         REAL (rprec)                            :: zmin
!>  Maximum Z position.
         REAL (rprec)                            :: zmax

!>  Radial grid size.
         REAL (rprec)                            :: dr
!>  Vertical grid size.
         REAL (rprec)                            :: dz

!>  Number of field periods.
         INTEGER                                 :: nfp

!>  Current density in the X direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_r => null()
!>  Current density in the Y direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_p => null()
!>  Current density in the Z direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_z => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object with an
!>  This computes the currents and positions on the primed grid.
!>
!>  @param[in] mgrid_file_name File name for vacuum fields.
!>  @returns A pointer to a constructed @ref m_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION m_grid_construct(mgrid_file_name)
      USE ezcdf
      USE read_wout_mod, Only: extcur
      IMPLICIT NONE

!  Declare Arguments
      TYPE (m_grid_class), POINTER                :: m_grid_construct
      CHARACTER (len=*), INTENT(in)               :: mgrid_file_name

!  local variables
      REAL (rprec)                                :: start_time
      INTEGER                                     :: mgrid_ncid
      INTEGER                                     :: i
      INTEGER                                     :: status
      INTEGER                                     :: varid
      INTEGER                                     :: num_r
      INTEGER                                     :: num_p
      INTEGER                                     :: num_z
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp_buffer
      CHARACTER (len=6)                           :: temp_string

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(m_grid_construct)

      CALL cdf_open(mgrid_ncid, TRIM(mgrid_file_name), 'r', status)

      CALL cdf_read(mgrid_ncid, 'ir', num_r)
      CALL cdf_read(mgrid_ncid, 'kp', num_p)
      CALL cdf_read(mgrid_ncid, 'jz', num_z)

      CALL cdf_read(mgrid_ncid, 'nfp', m_grid_construct%nfp)

      CALL cdf_read(mgrid_ncid, 'rmax', m_grid_construct%rmax)
      CALL cdf_read(mgrid_ncid, 'zmax', m_grid_construct%zmax)
      CALL cdf_read(mgrid_ncid, 'rmin', m_grid_construct%rmin)
      CALL cdf_read(mgrid_ncid, 'zmin', m_grid_construct%zmin)

      m_grid_construct%dr = (m_grid_construct%rmax -                           &
     &                       m_grid_construct%rmin)/(num_r - 1.0)
      m_grid_construct%dz = (m_grid_construct%zmax -                           &
     &                       m_grid_construct%zmin)/(num_z - 1.0)

      ALLOCATE(m_grid_construct%a_r(num_r, num_z, num_p))
      ALLOCATE(m_grid_construct%a_p(num_r, num_z, num_p))
      ALLOCATE(m_grid_construct%a_z(num_r, num_z, num_p))

      ALLOCATE(temp_buffer(num_r, num_z, num_p))

      DO i = 1, SIZE(extcur)
         WRITE (temp_string, 1000) i
         CALL cdf_read(mgrid_ncid, temp_string, temp_buffer)
         m_grid_construct%a_r = m_grid_construct%a_r                           &
     &                        + temp_buffer*extcur(i)

         WRITE (temp_string, 1001) i
         CALL cdf_read(mgrid_ncid, temp_string, temp_buffer)
         m_grid_construct%a_p = m_grid_construct%a_p                           &
     &                        + temp_buffer*extcur(i)

         WRITE (temp_string, 1002) i
         CALL cdf_read(mgrid_ncid, temp_string, temp_buffer)
         m_grid_construct%a_z = m_grid_construct%a_z                           &
     &                        + temp_buffer*extcur(i)
      END DO

      DEALLOCATE(temp_buffer)

      CALL cdf_close(mgrid_ncid)

      WRITE (*,2000)

      CALL profiler_set_stop_time('m_grid_construct', start_time)

1000  FORMAT('ar_',i3.3)
1001  FORMAT('ap_',i3.3)
1002  FORMAT('az_',i3.3)

2000  FORMAT('M Grid Ready')

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref m_grid_class object.
!>
!>  Deallocates memory and uninitializes a @ref m_grid_class object.
!>
!>  @param[inout] this A @ref m_grid_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE m_grid_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (m_grid_class), POINTER :: this

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
