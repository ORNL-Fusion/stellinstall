!*******************************************************************************
!>  @file limiter_iso_T.f
!>  @brief Contains module @ref limiter_iso_T
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref limiter_iso.
!*******************************************************************************
      MODULE limiter_iso_T

      USE stel_kinds, only : rprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
      USE v3_utilities
#if defined (MPI_OPT)
      USE mpi
#endif
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) extcurz base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a limiter_iso function.
!>
!>  1) Limiter iso: limiter_iso
!>  Contains the data to define a scalar function of position, so that the
!>  iso-contour f=0 of the function corresponds to a geometric limit to the
!>  plasma. Also contains data to specify the minimum number poloidal points on
!>  the s=1 surface to use, and the toroidal planes on which to s=1 surface will
!>  be evaluated
!>
!>  Second attempt:
!>
!>  e = SUM_over_i(0,4)_j(0,4) [arz(i,j) (r - rc)^i (z - zc)^j]
!>  f = e / |grad(e)|
!>
!>  - @note
!>  1. For now, axisymmetric
!>  2. Easy to put in circles, ellipses, and planes
!>  3. Function f is approximately distance.
!-------------------------------------------------------------------------------
      TYPE limiter_iso
!>  Coefficients for the iso function.
         REAL(rprec), DIMENSION(0:4,0:4)    :: arz
!>  R offset for function.
         REAL(rprec)                        :: rc
!>  Z offset for function.
         REAL(rprec)                        :: zc
!>  Minimum number of poloidal angles for s=1 surface.
         INTEGER                            :: numin
!>  Values of toroidal angle at which to compute s=1 surface.
         REAL(rprec), DIMENSION(:), POINTER :: vgrid => null()
      END TYPE limiter_iso

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref limiter_iso types using
!>  @ref limiter_iso_construct_new or @ref limiter_iso_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE limiter_iso_construct
         MODULE PROCEDURE limiter_iso_construct_new,                           &
     &                    limiter_iso_construct_mpi
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a limiter iso function.
!>
!>  Allocates memory and initializes a @ref limiter_iso object.
!>
!>  @param[in] arz   Coefficients for the iso function.
!>  @param[in] rc    R offset for function.
!>  @param[in] zc    Z offset for function.
!>  @param[in] numin Minimum number of poloidal angles for s=1 surface.
!>  @param[in] vgrid Values of toroidal angle at which to compute s=1 surface.
!>  @returns A pointer to a constructed @ref limiter_iso object.
!-------------------------------------------------------------------------------
      FUNCTION limiter_iso_construct_new(arz, rc, zc, numin, vgrid)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (limiter_iso), POINTER :: limiter_iso_construct_new
      REAL(rprec), DIMENSION(0:4,0:4), INTENT(in) :: arz
      REAL(rprec), INTENT(in)                     :: rc
      REAL(rprec), INTENT(in)                     :: zc
      INTEGER, INTENT(in)                         :: numin
      REAL(rprec), DIMENSION(:), INTENT(in)       :: vgrid

!  local variables
      REAL (rprec)                                :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  First, destroy this
      ALLOCATE(limiter_iso_construct_new)
      
      limiter_iso_construct_new%arz = arz
      limiter_iso_construct_new%rc = rc
      limiter_iso_construct_new%zc = zc
      limiter_iso_construct_new%numin = numin

      ALLOCATE(limiter_iso_construct_new%vgrid(SIZE(vgrid)))
      limiter_iso_construct_new%vgrid = vgrid

      CALL profiler_set_stop_time('limiter_iso_construct_new',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref limiter_iso object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref limiter_iso object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref limiter_iso object.
!-------------------------------------------------------------------------------
      FUNCTION limiter_iso_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_iso), POINTER :: limiter_iso_construct_mpi
      INTEGER, INTENT(in)         :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                     :: error
      INTEGER                     :: num_v

!  Start of executable code
      ALLOCATE(limiter_iso_construct_mpi)

      CALL MPI_BCAST(limiter_iso_construct_mpi%arz, 25, MPI_REAL8, 0,          &
     &               intra_comm, error)
      CALL MPI_BCAST(limiter_iso_construct_mpi%rc, 1, MPI_REAL8, 0,            &
     &               intra_comm, error)
      CALL MPI_BCAST(limiter_iso_construct_mpi%zc, 1, MPI_REAL8, 0,            &
     &               intra_comm, error)
      CALL MPI_BCAST(limiter_iso_construct_mpi%numin, 1, MPI_INTEGER, 0,       &
     &               intra_comm, error)

      CALL MPI_BCAST(num_v, 1, MPI_INTEGER, 0, intra_comm, error)
      ALLOCATE(limiter_iso_construct_mpi%vgrid(num_v))
      CALL MPI_BCAST(limiter_iso_construct_mpi%vgrid, num_v, MPI_REAL8,        &
     &               0, intra_comm, error)

#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref limiter_iso object.
!>
!>  Deallocates memory and uninitializes a @ref limiter_iso object.
!>
!>  @param[inout] this A @ref limiter_iso instance.
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_iso_destruct(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (limiter_iso), POINTER :: this

!  Start of executable code

!  Get rid of all components
      IF (ASSOCIATED(this%vgrid)) THEN
         DEALLOCATE(this%vgrid)
         this%vgrid => null()
      ENDIF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Computes the iso value at an r, phi, z position.
!>
!>  Limiter iso function is computed as
!>
!>  e = SUM_over_i(0,4)_j(0,4) [arz(i,j) (r - rc)^i (z - zc)^j]
!>  f = e / |grad(e)|
!>
!>
!>  @param[in] this    A @ref limiter_iso instance.
!>  @param[in] rpz_arg The r, phi, z point to evaluate at.
!>  @returns The value of iso function at the r, phi, z point.
!-------------------------------------------------------------------------------
      FUNCTION limiter_iso_value(this, rpz_arg)
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: limiter_iso_value
      TYPE (limiter_iso), INTENT(in)         :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: rpz_arg

!  local variables
      REAL(rprec)                            :: e, er, ez, grad_e
      INTEGER                                :: i, j
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Compute function
      e = zero
      DO i = 0, 4
         DO j = 0, 4
            e = e + this%arz(i,j)                                              &
     &        * ((rpz_arg(1) - this%rc)**i)                                    &
     &        * ((rpz_arg(3) - this%zc)**j)
         END DO
      END DO

!  Compute de/dr
      er = zero
      DO i = 1, 4
         DO j = 0, 4
            er = er + this%arz(i,j)*i                                          &
     &         * ((rpz_arg(1) - this%rc)**(i-1))                               &
     &         * ((rpz_arg(3) - this%zc)**j)
         END DO
      END DO

!  Compute de/dr
      ez = zero
      DO i = 0, 4
         DO j = 1, 4
         ez = ez + this % arz(i,j)*j                                           &
     &      * ((rpz_arg(1) - this%rc)**i)                                      &
     &      * ((rpz_arg(3) - this%zc)**(j - 1))
         END DO
      END DO

!  Avoid divide by zero errors.
      grad_e = MAX(1.E-12, SQRT(er*er + ez*ez))
      
      limiter_iso_value = e/grad_e

      CALL profiler_set_stop_time('limiter_iso_value', start_time)

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
!>  @param[inout] this       A @ref limiter_iso instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_iso_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_iso), INTENT(inout) :: this
      INTEGER, INTENT(in)               :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                           :: error
      INTEGER                           :: num_v
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL MPI_BCAST(this%arz, 25, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%rc, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%zc, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%numin, 1, MPI_INTEGER, 0, intra_comm, error)

      num_v = SIZE(this%vgrid)
      CALL MPI_BCAST(num_v, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%vgrid, num_v, MPI_REAL8, 0, intra_comm, error)

      CALL profiler_set_stop_time('limiter_iso_value', start_time)
#endif

      END SUBROUTINE

      END MODULE
