!*******************************************************************************
!>  @file surface_context.f
!>  @brief Contains module @ref surface_context.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref surface_context_class. This module
!>  contains all the code necessary define a surface for cariddi to read a
!>  virtual surface from.
!*******************************************************************************

      MODULE surface_context
      USE stel_kinds, ONLY: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) surface context class.
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a surface context. This contains all memory needed
!>  to operate surface.
!-------------------------------------------------------------------------------
      TYPE :: surface_context_class
!>  Array of radial points on the current surface.
         REAL (rprec), DIMENSION(:), POINTER :: x
!>  Array of phi points on the current surface.
         REAL (rprec), DIMENSION(:), POINTER :: y
!>  Array of z points on the current surface.
         REAL (rprec), DIMENSION(:), POINTER :: z
!>  Normal component in r direction.
         REAL (rprec), DIMENSION(:), POINTER :: x_normal
!>  Normal component in phi direction.
         REAL (rprec), DIMENSION(:), POINTER :: y_normal
!>  Normal component in z direction.
         REAL (rprec), DIMENSION(:), POINTER :: z_normal
!>  Virtual Surface current in the r direction.
         REAL (rprec), DIMENSION(:), POINTER :: k_x
!>  Virtual Surface current in the phi direction.
         REAL (rprec), DIMENSION(:), POINTER :: k_y
!>  Virtual Surface current in the z direction.
         REAL (rprec), DIMENSION(:), POINTER :: k_z
      CONTAINS
         PROCEDURE, PASS :: set_current => surface_context_set_current
         PROCEDURE, PASS :: write => surface_context_write
         FINAL :: surface_context_destruct
      END TYPE

!-------------------------------------------------------------------------------
!>  Interface to constructors.
!-------------------------------------------------------------------------------
      INTERFACE surface_context_class
         MODULE PROCEDURE surface_context_construct
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref surface_context_class object.
!>
!>  Allocates memory and initializes a @ref surface_context_class object.
!>
!>  @param[in] wout_file_name    Path and name of the wout file.
!>  @param[in] surface_file_name Path defining the surface to compute fields to.
!>  @param[in] parallel          @ref bmw_parallel_context_class object
!>                               instance.
!>  @param[in] io_unit           Unit number to write messages to.
!>  @returns A pointer to a constructed @ref bmw_context_class object.
!-------------------------------------------------------------------------------
      FUNCTION surface_context_construct(wout_file_name,                       &
     &                                   surface_file_name, parallel,          &
     &                                   io_unit)
      USE bmw_parallel_context
      USE read_wout_mod, ONLY: read_wout_file, extcur
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      CLASS (surface_context_class), POINTER ::                                &
     &   surface_context_construct
      CHARACTER (len=*), INTENT(in)          :: wout_file_name
      CHARACTER (len=*), INTENT(in)          :: surface_file_name
      CLASS (bmw_parallel_context_class), INTENT(in) :: parallel
      INTEGER, INTENT(in)                    :: io_unit

!  local variables
      REAL (rprec)                           :: start_time
      INTEGER                                :: status
      INTEGER                                :: iou
      INTEGER                                :: size
      INTEGER                                :: i
      INTEGER                                :: ncid
      INTEGER                                :: size_dim
      INTEGER                                :: xid
      INTEGER                                :: yid
      INTEGER                                :: zid
      INTEGER                                :: xnid
      INTEGER                                :: ynid
      INTEGER                                :: znid

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(surface_context_construct)

      CALL read_wout_file(wout_file_name, status)
      IF (status .ne. 0) THEN
         IF (parallel%offset .eq. 0) THEN
            WRITE (io_unit,1000) TRIM(wout_file_name)
         END IF
         CALL bmw_parallel_context_abort(status)
      END IF
      IF (.not.ALLOCATED(extcur)) THEN
         IF (parallel%offset .eq. 0) THEN
            WRITE (io_unit,1001) TRIM(wout_file_name)
         END IF
         CALL bmw_parallel_context_abort(status)
      END IF

!  Parse the surface file.
      status = NF_OPEN(TRIM(surface_file_name), NF_NOWRITE, ncid)

      status = NF_INQ_DIMID(ncid, 'size', size_dim)
      status = NF_INQ_DIMLEN(ncid, size_dim, size)

      ALLOCATE(surface_context_construct%x(size))
      ALLOCATE(surface_context_construct%y(size))
      ALLOCATE(surface_context_construct%z(size))

      ALLOCATE(surface_context_construct%x_normal(size))
      ALLOCATE(surface_context_construct%y_normal(size))
      ALLOCATE(surface_context_construct%z_normal(size))

      ALLOCATE(surface_context_construct%k_x(size))
      ALLOCATE(surface_context_construct%k_y(size))
      ALLOCATE(surface_context_construct%k_z(size))

      status = NF_INQ_VARID(ncid, 'x', xid)
      status = NF_INQ_VARID(ncid, 'y', yid)
      status = NF_INQ_VARID(ncid, 'z', zid)
      status = NF_INQ_VARID(ncid, 'x_norm', xnid)
      status = NF_INQ_VARID(ncid, 'y_norm', ynid)
      status = NF_INQ_VARID(ncid, 'z_norm', znid)

      status = NF_GET_VARA_DOUBLE(ncid, xid, 1, size,                          &
     &                            surface_context_construct%x)
      status = NF_GET_VARA_DOUBLE(ncid, yid, 1, size,                          &
     &                            surface_context_construct%y)
      status = NF_GET_VARA_DOUBLE(ncid, zid, 1, size,                          &
     &                            surface_context_construct%z)
      status = NF_GET_VARA_DOUBLE(ncid, xnid, 1, size,                         &
     &                            surface_context_construct%x_normal)
      status = NF_GET_VARA_DOUBLE(ncid, ynid, 1, size,                         &
     &                            surface_context_construct%y_normal)
      status = NF_GET_VARA_DOUBLE(ncid, znid, 1, size,                         &
     &                            surface_context_construct%z_normal)

      CALL profiler_set_stop_time('surface_context_construct',                 &
     &                            start_time)

1000  FORMAT(a,' is an invalid wout file.')
1001  FORMAT(a,' is not a free boundary wout file.')

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
      SUBROUTINE surface_context_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (surface_context_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%x)) THEN
         DEALLOCATE(this%x)
         this%x => null()
      END IF

      IF (ASSOCIATED(this%y)) THEN
         DEALLOCATE(this%y)
         this%y => null()
      END IF

      IF (ASSOCIATED(this%z)) THEN
         DEALLOCATE(this%z)
         this%z => null()
      END IF

      IF (ASSOCIATED(this%x_normal)) THEN
         DEALLOCATE(this%x_normal)
         this%x_normal => null()
      END IF

      IF (ASSOCIATED(this%y_normal)) THEN
         DEALLOCATE(this%y_normal)
         this%y_normal => null()
      END IF

      IF (ASSOCIATED(this%z_normal)) THEN
         DEALLOCATE(this%z_normal)
         this%z_normal => null()
      END IF

      IF (ASSOCIATED(this%k_x)) THEN
         DEALLOCATE(this%k_x)
         this%k_x => null()
      END IF

      IF (ASSOCIATED(this%k_y)) THEN
         DEALLOCATE(this%k_y)
         this%k_y => null()
      END IF

      IF (ASSOCIATED(this%k_z)) THEN
         DEALLOCATE(this%k_z)
         this%k_z => null()
      END IF

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Set the virtual current on the surface.
!>
!>  This uses the biot-savart to set a vitual current on an arbitrary
!>  surface. The vitual current is taken from the jump condition assuming the
!>  surface is a perfect conductor. At each point in the surface, the plasma
!>  current in volumn integrated.
!>
!>    b = mu0/4pi Int Int Int (Jsupu*esubu + Jsupv*esubv) X
!>                            (x - xp)/(x - xp)*signJ*J*da*du*dv             (1)
!>
!>  Where the primed (p) variables are the values on the last closed flux
!>  surface grid.
!>
!>  @param[inout] this     A @ref surface_context_class instance.
!>  @param[in]    io_unit  Unit number to write messages to.
!>  @param[in]    parallel @ref bmw_parallel_context_class object instance.
!>  @returns A pointer to a constructed @ref bmw_context_class object.
!-------------------------------------------------------------------------------
      SUBROUTINE surface_context_set_current(this, io_unit, parallel)
      USE bmw_parallel_context
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: ns, rmnc, rmns, zmnc, zmns, xm, xn,             &
     &                         bsubumnc, bsubvmnc, bsubumns, bsubvmns,         &
     &                         xm_nyq, xn_nyq, mnmax, mnmax_nyq, lasym,        &
     &                         isigng
      USE bmw_state_flags, ONLY: clear_screen, progress,                       &
     &                           bmw_state_flags_off
      USE, INTRINSIC :: iso_fortran_env, Only : output_unit
      USE primed_grid
!$    USE omp_lib

!  Declare Arguments
      CLASS (surface_context_class), INTENT(inout)   :: this
      INTEGER, INTENT(in)                            :: io_unit
      CLASS (bmw_parallel_context_class), INTENT(in) :: parallel

!  local variables
      REAL (rprec)                                   :: start_time
      REAL (rprec)                                   :: total
      REAL (rprec)                                   :: done
      REAL (rprec)                                   :: current
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE    :: lx
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE    :: ly
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE    :: lz
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE    :: l
      REAL (rprec)                                   :: bx
      REAL (rprec)                                   :: by
      REAL (rprec)                                   :: bz
      INTEGER                                        :: i
      CLASS (primed_grid_class), POINTER             :: p_grid

!  Start of executable code
      start_time = profiler_get_start_time()

      p_grid => primed_grid_class(101, bmw_state_flags_off, '',                &
     &                            parallel, io_unit)

      total = parallel%end(SIZE(this%x))
      total = CEILING(total/parallel%num_threads)
      current = 0.0

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i, lx, ly, lz, l, bx, by, bz)

      ALLOCATE(lx(SIZE(p_grid%x, 1),                                           &
     &            SIZE(p_grid%x, 2),                                           &
     &            SIZE(p_grid%x, 3)))
      ALLOCATE(ly(SIZE(p_grid%x, 1),                                           &
     &            SIZE(p_grid%x, 2),                                           &
     &            SIZE(p_grid%x, 3)))
      ALLOCATE(lz(SIZE(p_grid%x, 1),                                           &
     &            SIZE(p_grid%x, 2),                                           &
     &            SIZE(p_grid%x, 3)))
      ALLOCATE(l(SIZE(p_grid%x, 1),                                            &
     &           SIZE(p_grid%x, 2),                                            &
     &           SIZE(p_grid%x, 3)))

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = parallel%start(SIZE(this%x)), parallel%end(SIZE(this%x))
         lx = this%x(i) - p_grid%x
         ly = this%y(i) - p_grid%y
         lz = this%z(i) - p_grid%z

         l = SQRT(lx*lx + ly*ly + lz*lz)**3.0

         bx = SUM((p_grid%j_y*lz - p_grid%j_z*ly)/l)*p_grid%dvol
         by = SUM((p_grid%j_z*lx - p_grid%j_x*lz)/l)*p_grid%dvol
         bz = SUM((p_grid%j_x*ly - p_grid%j_y*lx)/l)*p_grid%dvol

!  NOTE: The proper VC current for the coupling surface requires a minus sign
!        here (-n_hat X B_in). However, the field produced by this current is
!        meant to cancel out the field due to currents enclosed so B_vc = -B_im.
!        For the coupling ignore the minus so the currents can be used without
!        adding additonal signs.
         this%k_x(i) = (this%y_normal(i)*bz - this%z_normal(i)*by)/mu0
         this%k_y(i) = (this%z_normal(i)*bx - this%x_normal(i)*bz)/mu0
         this%k_z(i) = (this%x_normal(i)*by - this%y_normal(i)*bx)/mu0

!  Progress Reporting
         IF (parallel%offset .eq. 0) THEN
!$          IF (OMP_GET_THREAD_NUM() .eq. 0) THEN
               current = current + 1.0
               done = 100.0*current/total

               WRITE (io_unit,1000,ADVANCE='NO')                               &
     &            clear_screen, progress(MOD(INT(current),4)), done

               IF (io_unit .ne. output_unit) THEN
                  BACKSPACE (io_unit)
               END IF
!$          END IF
         END IF
      END DO
!$OMP END DO

      DEALLOCATE(lx)
      DEALLOCATE(ly)
      DEALLOCATE(lz)
      DEALLOCATE(l)

!$OMP END PARALLEL

      IF (parallel%stride .gt. 1) THEN
         CALL parallel%reduce(this%k_x)
         CALL parallel%reduce(this%k_y)
         CALL parallel%reduce(this%k_z)
      END IF

      DEALLOCATE(p_grid)
      p_grid => null()

      IF (parallel%offset .eq. 0) THEN
         WRITE (io_unit,1001) clear_screen
      END IF

      CALL profiler_set_stop_time('surface_context_set_current',               &
     &                            start_time)

1000  FORMAT(a,a,f6.2,' % Finished')
1001  FORMAT(a,'Coupling Grid Finished')

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
      SUBROUTINE surface_context_write(this, output_file, io_unit,             &
     &                                 parallel)
      USE bmw_parallel_context
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      CLASS (surface_context_class), INTENT(inout)   :: this
      CHARACTER (len=*), INTENT(in)                  :: output_file
      INTEGER, INTENT(in)                            :: io_unit
      CLASS (bmw_parallel_context_class), INTENT(in) :: parallel

!  local variables
      INTEGER                                        :: iou
      INTEGER                                        :: i
      INTEGER                                        :: status
      REAL (rprec)                                   :: start_time
      INTEGER                                        :: ncid
      INTEGER                                        :: size_dim
      INTEGER                                        :: kxid
      INTEGER                                        :: kyid
      INTEGER                                        :: kzid

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (parallel%offset .eq. 0) THEN
         status = NF_OPEN(TRIM(output_file), NF_WRITE, ncid)

         status = NF_INQ_DIMID(ncid, 'size', size_dim)

         status = NF_INQ_VARID(ncid, 'k_x', kxid)
         IF (status .ne. NF_NOERR) THEN
            status = NF_REDEF(ncid)
            status = NF_DEF_VAR(ncid, 'k_x', NF_DOUBLE, 1, size_dim,           &
     &                          kxid)
         END IF
         status = NF_INQ_VARID(ncid, 'k_y', kyid)
         IF (status .ne. NF_NOERR) THEN
            status = NF_REDEF(ncid)
            status = NF_DEF_VAR(ncid, 'k_y', NF_DOUBLE, 1, size_dim,           &
     &                          kyid)
         END IF
         status = NF_INQ_VARID(ncid, 'k_z', kzid)
         IF (status .ne. NF_NOERR) THEN
            status = NF_REDEF(ncid)
            status = NF_DEF_VAR(ncid, 'k_z', NF_DOUBLE, 1, size_dim,           &
     &                          kzid)
         END IF

         status = NF_ENDDEF(ncid)

         status = NF_PUT_VARA_DOUBLE(ncid, kxid, 1, SIZE(this%x),              &
     &                               this%k_x)
         status = NF_PUT_VARA_DOUBLE(ncid, kyid, 1, SIZE(this%x),              &
     &                               this%k_y)
         status = NF_PUT_VARA_DOUBLE(ncid, kzid, 1, SIZE(this%x),              &
     &                               this%k_z)

         status = NF_CLOSE(ncid)
      END IF

      CALL profiler_set_stop_time('surface_context_write', start_time)

      END SUBROUTINE

      END MODULE
