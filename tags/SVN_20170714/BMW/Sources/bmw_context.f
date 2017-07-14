!*******************************************************************************
!>  @file bmw_context.f
!>  @brief Contains module @ref bmw_context.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref bmw_context_class. This contains the
!>  state variables needed by BMW.
!*******************************************************************************
      MODULE bmw_context
      USE m_grid
      USE primed_grid
      USE unprimed_grid

      IMPLICIT NONE

!*******************************************************************************
!  bmw context module parameters
!*******************************************************************************
!>  Version number.
      INTEGER, PARAMETER :: series = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) bmw context base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a bmw context. This contains all memory needed to
!>  operate bmw.
!-------------------------------------------------------------------------------
      TYPE bmw_context_class
!>  Vacuum grid object.
         TYPE (m_grid_class), POINTER        :: m_grid => null()
!>  Primed grid object.
         TYPE (primed_grid_class), POINTER   :: p_grid => null()
!>  Unprimed grid object.
         TYPE (unprimed_grid_class), POINTER :: up_grid => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref bmw_context_class object.
!>
!>  Allocates memory and initializes a @ref bmw_context_class object.
!>
!>  @param[in] mgrid_file_name Path and name of the mgrid file.
!>  @param[in] wout_file_name  Path and name of the wout file.
!>  @param[in] flags           Option flags.
!>  @returns A pointer to a constructed @ref bmw_context_class object.
!-------------------------------------------------------------------------------
      FUNCTION bmw_context_construct(mgrid_file_name, wout_file_name,          &
     &                               flags)
      USE read_wout_mod, Only: read_wout_file

      IMPLICIT NONE

!  Declare Arguments
      TYPE (bmw_context_class), POINTER :: bmw_context_construct
      CHARACTER (len=*), INTENT(in)     :: mgrid_file_name
      CHARACTER (len=*), INTENT(in)     :: wout_file_name
      INTEGER, INTENT(in)               :: flags

!  local variables
      REAL (rprec)                      :: start_time
      INTEGER                           :: status
      INTEGER                           :: num_p

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL read_wout_file(wout_file_name, status)
      IF (status .ne. 0) THEN
         WRITE (*,1000) wout_file_name
         CALL EXIT(status)
      END IF

      ALLOCATE(bmw_context_construct)

      bmw_context_construct%m_grid => m_grid_construct(mgrid_file_name)

      num_p = SIZE(bmw_context_construct%m_grid%a_p, 3)

      bmw_context_construct%p_grid =>                                          &
     &   primed_grid_construct(num_p*bmw_context_construct%m_grid%nfp,         &
     &                         flags)
      bmw_context_construct%up_grid =>                                         &
     &   unprimed_grid_construct(bmw_context_construct%m_grid,                 &                                        &
     &                           bmw_context_construct%p_grid)

      CALL profiler_set_stop_time('bmw_context_construct', start_time)

1000  FORMAT(a,' is an invalid wout file.')

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref bmw_context_class object.
!>
!>  Deallocates memory and uninitializes a @ref bmw_context_class object.
!>
!>  @param[inout] this A @ref bmw_context_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE bmw_context_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (bmw_context_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%m_grid)) THEN
         CALL m_grid_destruct(this%m_grid)
         this%m_grid => null()
      END IF

      IF (ASSOCIATED(this%p_grid)) THEN
         CALL primed_grid_destruct(this%p_grid)
         this%p_grid => null()
      END IF

      IF (ASSOCIATED(this%up_grid)) THEN
         CALL unprimed_grid_destruct(this%up_grid)
         this%up_grid => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!>  @page result_file_bmw BMW Result File
!>
!>  @tableofcontents
!>  @section result_file_bmw_intro_sec Introduction
!>  This page documents the contents of a result NetCDF.
!>
!>  @section result_file_vmec_var_sec Variables
!>  @header{Variable(Dimensions), Description, Code Reference}
!>  @begin_table
!>     @item{series,                          Version number.,                    bmw_context::series}
!>     @item{nfp,                             Number of field periods,            m_grid_class::m_grid_class::nfp}
!>     @item{rmin,                            Minimum R point.,                   m_grid_class::m_grid_class::rmin}
!>     @item{rmax,                            Maximum R point.,                   m_grid_class::m_grid_class::rmax}
!>     @item{zmin,                            Minimum Z point.,                   m_grid_class::m_grid_class::rmin}
!>     @item{zmax,                            Maximum Z point.,                   m_grid_class::m_grid_class::zmax}
!>     @item{ar_grid (num_p\, num_z\, num_r), Vector potential in R direction.,   unprimed_grid::unprimed_grid_class::ar_grid}
!>     @item{ap_grid (num_p\, num_z\, num_r), Vector potential in Phi direction., unprimed_grid::unprimed_grid_class::ap_grid}
!>     @item{az_grid (num_p\, num_z\, num_r), Vector potential in Z direction.,   unprimed_grid::unprimed_grid_class::az_grid}
!>  @end_table
!-------------------------------------------------------------------------------
!>  @brief Write NetCDF based result file
!>
!>  Defines dimensions and variables for the BMW result file. Flat internal
!>  arrays need to be reshaped to the dimension.
!>
!>  @param[in] this             A @ref bmw_context_class instance.
!>  @param[in] result_file_name NetCDF file id of the result file.
!-------------------------------------------------------------------------------
      SUBROUTINE bmw_context_write(this, result_file_name)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (bmw_context_class), INTENT(in)        :: this
      CHARACTER (len=*), INTENT(in)               :: result_file_name

!  local variables
      REAL (rprec)                                :: start_time
      INTEGER                                     :: status
      INTEGER                                     :: result_iou

!  local parameters
      CHARACTER (len=*), DIMENSION(3), PARAMETER  ::                           &
     &   dims = (/ 'r  ','z  ','phi' /)

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL cdf_open(result_iou, TRIM(result_file_name), 'w', status)

      CALL cdf_define(result_iou, 'series', series)

      CALL cdf_define(result_iou, 'nfp', this%m_grid%nfp)

      CALL cdf_define(result_iou, 'rmin', this%m_grid%rmin)
      CALL cdf_define(result_iou, 'rmax', this%m_grid%rmax)
      CALL cdf_define(result_iou, 'zmin', this%m_grid%zmin)
      CALL cdf_define(result_iou, 'zmax', this%m_grid%zmax)

      CALL cdf_define(result_iou, 'ar_grid', this%up_grid%a_r,                 &
     &                dimname=dims)
      CALL cdf_define(result_iou, 'ap_grid', this%up_grid%a_p,                 &
     &                dimname=dims)
      CALL cdf_define(result_iou, 'az_grid', this%up_grid%a_z,                 &
     &                dimname=dims)

      CALL cdf_define(result_iou, 'br_grid', this%up_grid%b_r,                 &
     &                dimname=dims)
      CALL cdf_define(result_iou, 'bp_grid', this%up_grid%b_p,                 &
     &                dimname=dims)
      CALL cdf_define(result_iou, 'bz_grid', this%up_grid%b_z,                 &
     &                dimname=dims)

      CALL cdf_write(result_iou, 'series', series)

      CALL cdf_write(result_iou, 'nfp', this%m_grid%nfp)

      CALL cdf_write(result_iou, 'rmin', this%m_grid%rmin)
      CALL cdf_write(result_iou, 'rmax', this%m_grid%rmax)
      CALL cdf_write(result_iou, 'zmin', this%m_grid%zmin)
      CALL cdf_write(result_iou, 'zmax', this%m_grid%zmax)

      CALL cdf_write(result_iou, 'ar_grid', this%up_grid%a_r)
      CALL cdf_write(result_iou, 'ap_grid', this%up_grid%a_p)
      CALL cdf_write(result_iou, 'az_grid', this%up_grid%a_z)

      CALL cdf_write(result_iou, 'br_grid', this%up_grid%b_r)
      CALL cdf_write(result_iou, 'bp_grid', this%up_grid%b_p)
      CALL cdf_write(result_iou, 'bz_grid', this%up_grid%b_z)

      CALL cdf_close(result_iou)

      CALL profiler_set_stop_time('bmw_context_write', start_time)

      END SUBROUTINE

      END MODULE
