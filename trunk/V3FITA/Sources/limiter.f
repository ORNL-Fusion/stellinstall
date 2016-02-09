!*******************************************************************************
!>  @file limiter.f
!>  @brief Contains module @ref limiter
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref limiter_class.
!>  @par Super Class:
!>  @ref geometric
!>
!>  @see limiter_iso_T
!*******************************************************************************

      MODULE limiter

      USE stel_kinds, only: rprec, dp
      USE limiter_iso_T
      USE limiter_grid

      IMPLICIT NONE
!*******************************************************************************
!  geometric module parameters
!*******************************************************************************
!>  Type descriptor for limiter no type.
      INTEGER, PARAMETER :: limiter_no_type      = -1
!>  Type descriptor for limiter type defined by an iso function.
      INTEGER, PARAMETER :: limiter_iso_type     = 0
!>  Type descriptor for limiter type defined by polygons.
!>  @note This limiter type has not be implemented yet.
      INTEGER, PARAMETER :: limiter_polygon_type = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) limiter base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a limiter signal.
!>  @par Super class:
!>  @ref geometric
!-------------------------------------------------------------------------------
      TYPE limiter_class
!>  Type descirptor of the limiter type.
!>  @par Possible values are:
!>  * @ref limiter_no_type
!>  * @ref limiter_iso_type
!>  * @ref limiter_polygon_type
         INTEGER                            :: type = limiter_no_type
!>  If true, the reconstruction tries to place the edge on the limiter.
!>  Otherwise, the plasma edge is kept inside the limiter.
         LOGICAL                            :: on_edge = .false.
!>  Iso surface function discribing the limiter. Not used for polygon limiters.
         TYPE (limiter_iso),       POINTER  :: iso => null()
!>  Grid function discribing the limiter. Not used for iso limiters.
         TYPE (limiter_grid_class), POINTER :: grid => null()
      END TYPE limiter_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref limiter_class types either using
!>  @ref limiter_construct_iso, @ref limiter_construct_polygon or
!>  @ref limiter_construct_mpi
!>
!>  @note
!>  @ref limiter_construct_polygon has not been implemented yet.
!-------------------------------------------------------------------------------
      INTERFACE limiter_construct
         MODULE PROCEDURE limiter_construct_iso,                               &
     &                    limiter_construct_polygon,                           &
     &                    limiter_construct_mpi
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a limiter containing an iso function.
!>
!>  Allocates memory and initializes a @ref limiter_class object.
!>
!>  @param[in] iso         A contour function to define the limiter.
!>  @param[in] on_edge     Specifies if the edge should touch the limiter or
!>                         not. If true, the reconstruction tries to have the
!>                         limiter touch the limiter. If false, the edge will
!>                         only just fall inside the limiter.
!>  @returns A pointer to a constructed @ref limiter_class object.
!-------------------------------------------------------------------------------
      FUNCTION limiter_construct_iso(iso, on_edge)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), POINTER           :: limiter_construct_iso
      TYPE (limiter_iso), POINTER             :: iso
      LOGICAL, INTENT(in)                     :: on_edge

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(limiter_construct_iso)

      limiter_construct_iso%type = limiter_iso_type
      limiter_construct_iso%on_edge = on_edge
      limiter_construct_iso%iso => iso

      CALL profiler_set_stop_time('limiter_construct_iso', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a limiter containing a polygon.
!>
!>  Allocates memory and initializes a @ref limiter_class object.
!>
!>  @param[in] grid    Grid of limiter contours.
!>  @param[in] on_edge Specifies if the edge should touch the limiter or not. If
!>                     true, the reconstruction tries to have the limiter touch
!>                     the limiter. If false, the edge will only just fall
!>                     inside the limiter.
!>  @returns A pointer to a constructed @ref limiter_class object.
!>  @note Polygon limiters have not be implemented yet.
!-------------------------------------------------------------------------------
      FUNCTION limiter_construct_polygon(grid, on_edge)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), POINTER      :: limiter_construct_polygon
      TYPE (limiter_grid_class), POINTER :: grid
      LOGICAL, INTENT(in)                :: on_edge

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(limiter_construct_polygon)

      limiter_construct_polygon%type = limiter_polygon_type
      limiter_construct_polygon%on_edge = on_edge

      limiter_construct_polygon%grid => grid

      CALL profiler_set_stop_time('limiter_construct_polygon',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref limiter_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref limiter_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref limiter_class object.
!-------------------------------------------------------------------------------
      FUNCTION limiter_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), POINTER :: limiter_construct_mpi
      INTEGER, INTENT(in)           :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                       :: error

!  Start of executable code
      ALLOCATE(limiter_construct_mpi)

      CALL MPI_BCAST(limiter_construct_mpi%type, 1, MPI_INTEGER, 0,            &
     &               intra_comm, error)
      CALL MPI_BCAST(limiter_construct_mpi%on_edge, 1, MPI_LOGICAL, 0,         &
     &               intra_comm, error)

      SELECT CASE (limiter_construct_mpi%type)

         CASE (limiter_iso_type)
            limiter_construct_mpi%iso =>                                       &
     &         limiter_iso_construct(intra_comm)

         CASE (limiter_polygon_type)
            limiter_construct_mpi%grid =>                                      &
     &         limiter_grid_construct(intra_comm)

      END SELECT

#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref limiter_class object.
!>
!>  Deallocates memory and uninitializes a @ref limiter_class object.
!> 
!>  @param[inout] this A @ref limiter_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%iso)) THEN
         CALL limiter_iso_destruct(this%iso)
         this%iso => null()
      END IF

      IF (ASSOCIATED(this%grid)) THEN
         CALL limiter_grid_destruct(this%grid)
         this%grid => null()
      END IF

      this%type = limiter_no_type
      this%on_edge = .false.

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  If the limiter is defined by an iso function returns the value of
!>  @ref limiter_get_modeled_iso_signal. Otherwise returns zero.
!>
!>  @param[in]    this    A @ref limiter_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   If the limiter type is defined by an iso function,
!>                        sigma value of @ref limiter_get_modeled_iso_signal.
!>  @returns If the limiter type is defined by an iso function, return the value
!>           of @ref limiter_get_modeled_iso_signal. Otherwise returns zero.
!-------------------------------------------------------------------------------
      FUNCTION limiter_get_modeled_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: limiter_get_modeled_signal
      TYPE (limiter_class), INTENT(in)        :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (this%type)

         CASE (limiter_iso_type)
            limiter_get_modeled_signal =                                       &
     &         limiter_get_modeled_iso_signal(this, a_model, sigma)

         CASE (limiter_polygon_type)
            limiter_get_modeled_signal =                                       &
     &         limiter_get_modeled_polygon_signal(this, a_model, sigma)

      END SELECT

      CALL profiler_set_stop_time('limiter_get_modeled_signal',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the limiter type.
!>
!>  Returns a description of the limiter type for use when writting output files.
!>
!>  @param[in] this A @ref limiter_class instance.
!>  @returns A string describing the limiter type.
!-------------------------------------------------------------------------------
      FUNCTION limiter_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: limiter_get_signal_type
      TYPE (limiter_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (this%type)

         CASE (limiter_iso_type)
            limiter_get_signal_type = 'edge_lim iso'

         CASE (limiter_polygon_type)
            limiter_get_signal_type = 'edge_lim poly'

      END SELECT

      CALL profiler_set_stop_time('limiter_get_signal_type', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  Returns a description of the array indices for use when writting output
!>  files.
!>
!>  @param[in] this A @ref limiter_class instance.
!>  @returns A string describing the model and model sigma array indices.
!>  @note Polygon limiters have not be implemented yet. This coding may change
!>  when those are supported.
!-------------------------------------------------------------------------------
      FUNCTION limiter_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   limiter_get_header
      TYPE (limiter_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      limiter_get_header(1) = 'r (m)'
      limiter_get_header(2) = 'phi (rad)'
      limiter_get_header(3) = 'z (m)'
      limiter_get_header(4) = 'model_sig(1)'
      limiter_get_header(5) = 'model_sig(2)'
      limiter_get_header(6) = 'model_sig(3)'
      limiter_get_header(7) = 'model_sig(4)'

      CALL profiler_set_stop_time('limiter_get_header', start_time)

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
!>  @param[inout] this       A @ref limiter_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                 :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                             :: error
      REAL (rprec)                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL MPI_BCAST(this%type, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%on_edge, 1, MPI_LOGICAL, 0, intra_comm, error)

      SELECT CASE (this%type)

         CASE (limiter_iso_type)
            CALL limiter_iso_sync_children(this%iso, intra_comm)

         CASE (limiter_polygon_type)
            CALL limiter_grid_sync_children(this%grid, intra_comm)

      END SELECT

      CALL profiler_set_stop_time('limiter_sync_children', start_time)
#endif

      END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal from an iso function.
!>
!>  If @ref limiter_class::on_edge or the plasma edge extends beyond the
!>  limiter, returns the R, Phi and Z position of the maximum of the iso
!>  function value. Otherwise returns zero.
!>  @par
!>  The R and Z positions for a fixed toroidal plane are retrived from the the
!>  @ref equilibrium method @ref equilibrium_get_plasma_edge
!>
!>  @param[in]    this    A @ref limiter_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   Returns the value of modeled sigma.
!>  @returns If @ref limiter_class::on_edge or the plasma edge extends beyond
!>  the limiter, returns the R, Phi and Z position of the maximum of the iso
!>  function value. Otherwise returns zero.
!-------------------------------------------------------------------------------
      FUNCTION limiter_get_modeled_iso_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: limiter_get_modeled_iso_signal
      TYPE (limiter_class), INTENT(in)        :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  local variables
      INTEGER                                 :: phi_index, theta_index
      INTEGER                                 :: num_theta
!  The r and z arrays will be allocated by the equilibrium.
      REAL(rprec), DIMENSION(:), POINTER      :: r, z
      REAL(rprec), DIMENSION(3)               :: rphiz_at_max
      REAL(rprec)                             :: fval_max, fval
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      r => null()
      z => null()

      limiter_get_modeled_iso_signal = 0.0
      rphiz_at_max = 0.0
      fval_max = -1.0E10

!  Loop over each phi plane. r, z are allocated by equilibrium_get_plasma_edge.
!  There is no need to allocate and deallocate r and z for each phi angle.
      DO phi_index = 1, SIZE(this%iso%vgrid)
         num_theta =                                                           &
     &      equilibrium_get_plasma_edge(a_model%equilibrium,                   &
     &                                  this%iso%vgrid(phi_index),             &
     &                                  r, z)

!  Loop over each theta angle
         DO theta_index = 1, num_theta
            fval = limiter_iso_value(this%iso,                                 &
     &                               (/ r(theta_index),                        &
     &                                  this%iso%vgrid(phi_index),             &
     &                                  z(theta_index) /))

            IF (fval .gt. fval_max) THEN
               fval_max = fval
               rphiz_at_max = (/ r(theta_index),                               &
     &                           this%iso%vgrid(phi_index),                    &
     &                           z(theta_index) /)
            END IF
         END DO
      END DO

      IF (this%on_edge) THEN
         limiter_get_modeled_iso_signal(1) = fval_max
      ELSE
         limiter_get_modeled_iso_signal(1) = MAX(fval_max, 0.0)
      END IF
      limiter_get_modeled_iso_signal(2:4) = rphiz_at_max

!  Dealloctae r and z since they are no longer needed.
      IF (ASSOCIATED(r)) THEN
         DEALLOCATE(r)
      END IF
      IF (ASSOCIATED(z)) THEN
         DEALLOCATE(z)
      END IF
      r => null()
      z => null()

      sigma = 0.0

      CALL profiler_set_stop_time('limiter_get_modeled_iso_signal',            &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal from an iso grid.
!>
!>  If @ref limiter_class::on_edge or the plasma edge extends beyond the
!>  limiter, returns the R, Phi and Z position of the maximum of the iso
!>  grid value. Otherwise returns zero.
!>  @par
!>  The R and Z positions for a fixed toroidal plane are retrived from the the
!>  @ref equilibrium method @ref equilibrium_get_plasma_edge
!>
!>  @param[in]    this    A @ref limiter_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   Returns the value of modeled sigma.
!>  @returns If @ref limiter_class::on_edge or the plasma edge extends beyond
!>  the limiter, returns the R, Phi and Z position of the maximum of the iso
!>  grid value. Otherwise returns zero.
!-------------------------------------------------------------------------------
      FUNCTION limiter_get_modeled_polygon_signal(this, a_model, sigma)
      USE model
      USE bivariate

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: limiter_get_modeled_polygon_signal
      TYPE (limiter_class), INTENT(in)        :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  local variables
      TYPE (bivariate_type), POINTER          :: bivariate_object
      INTEGER                                 :: phi_index, index
      INTEGER                                 :: num_theta
!  The r and z arrays will be allocated by the equilibrium.
      REAL (rprec), DIMENSION(:), POINTER     :: r, z
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: rgrid, zgrid
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: fval
      INTEGER                                 :: num_r, num_z
      REAL (rprec), DIMENSION(3)              :: rphiz_at_max
      REAL(rprec)                             :: fval_max
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      r => null()
      z => null()

!  Get the r and z axis for the limiter grid.
      num_r = limiter_grid_get_num_r(this%grid)
      num_z = limiter_grid_get_num_z(this%grid)

      ALLOCATE(rgrid(num_r))
      DO index = 1, num_r
         rgrid(index) = (index - 1)*this%grid%dr + this%grid%r0
      END DO

      ALLOCATE(zgrid(num_z))
      DO index = 1, num_z
         zgrid(index) = (index - 1)*this%grid%dz + this%grid%z0
      END DO

      limiter_get_modeled_polygon_signal = 0.0
      rphiz_at_max = 0.0
      fval_max = -1.0E10

!  Loop over each phi plane. r, z are allocated by equilibrium_get_plasma_edge.
!  There is no need to allocate and deallocate r and z for each phi angle.
      DO phi_index = 1, SIZE(this%grid%phi)
         num_theta =                                                           &
     &      equilibrium_get_plasma_edge(a_model%equilibrium,                   &
     &                                  this%grid%phi(phi_index),              &
     &                                  r, z)

         ALLOCATE(fval(num_theta))
         bivariate_object => bivariate_construct(1, num_theta)

         CALL bivariate_set_grids(bivariate_object, r, z, rgrid, zgrid)
         CALL bivariate_get_4pt(bivariate_object,                              &
     &                          this%grid%grid(:,:,phi_index), fval)

!  Find the largest distance away from the
         index = MAXLOC(fval, 1)
         IF (fval(index) .gt. fval_max) THEN
            fval_max = fval(index)
            rphiz_at_max = (/ r(index), this%grid%phi(phi_index),              &
     &                        z(index) /)
         END IF

         CALL bivariate_destruct(bivariate_object)
         DEALLOCATE(fval)
      END DO

      IF (this%on_edge) THEN
         limiter_get_modeled_polygon_signal(1) = fval_max
      ELSE
         limiter_get_modeled_polygon_signal(1) = MAX(fval_max, 0.0)
      END IF
      limiter_get_modeled_polygon_signal(2:4) = rphiz_at_max

!  Dealloctae r and z since they are no longer needed.
      IF (ASSOCIATED(r)) THEN
         DEALLOCATE(r)
      END IF
      IF (ASSOCIATED(z)) THEN
         DEALLOCATE(z)
      END IF
      r => null()
      z => null()

      DEALLOCATE(rgrid)
      DEALLOCATE(zgrid)

      sigma = 0.0

      CALL profiler_set_stop_time('limiter_get_modeled_polygon_signal',        &
     &                            start_time)

      END FUNCTION

      END MODULE
