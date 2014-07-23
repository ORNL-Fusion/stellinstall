!*******************************************************************************
!>  @file extcurz.f
!>  Contains module @ref extcurz.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Represents a signal obtained by integrating around the magnetic field to get
!>  the current enclosed in the Z direction using Ampere's Law. Defines the base
!>  class of the type @ref extcurz_class.
!>  @par Super Class:
!>  @ref diagnostic
!*******************************************************************************

      MODULE extcurz

      USE stel_kinds, only: rprec
#if defined (MPI_OPT)
      USE mpi
#endif

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) extcurz base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a external Z currents enclosed by an inegration
!>  loop.
!>  @par Super Class:
!>  @ref diagnostic
!-------------------------------------------------------------------------------
      TYPE extcurz_class
         REAL (rprec) :: r
         REAL (rprec) :: theta
      END TYPE extcurz_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref extcurz_class types using
!>  @ref extcurz_construct_new or @ref extcurz_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE extcurz_construct
         MODULE PROCEDURE extcurz_construct_new,                               &
     &                    extcurz_construct_mpi
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref extcurz_class object measureing temperature.
!>
!>  Allocates memory and initializes a @ref extcurz_class object.
!>
!>  @param[in] r     The radial position to integrate at.
!>  @param[in] theta The poloidal angle to integrate at.
!>  @returns A pointer to a constructed @ref thomson_class object.
!-------------------------------------------------------------------------------
      FUNCTION extcurz_construct_new(r, theta)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (extcurz_class), POINTER  :: extcurz_construct_new
      REAL(rprec), INTENT(in)        :: r
      REAL(rprec), INTENT(in)        :: theta

!  Start of executable code
      ALLOCATE(extcurz_construct_new)
      extcurz_construct_new%r = r
      extcurz_construct_new%theta = theta

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref extcurz_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref extcurz_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref extcurz_class object.
!-------------------------------------------------------------------------------
      FUNCTION extcurz_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (extcurz_class), POINTER :: extcurz_construct_mpi
      INTEGER, INTENT(in)           :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                       :: error

!  Start of executable code
      ALLOCATE(extcurz_construct_mpi)

      CALL MPI_BCAST(extcurz_construct_mpi%r, 1, MPI_REAL8, 0,                 &
     &               intra_comm, error)
      CALL MPI_BCAST(extcurz_construct_mpi%theta, 1, MPI_REAL8, 0,             &
     &               intra_comm, error)

#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref extcurz_class object.
!>
!>  Deallocates memory and uninitializes a @ref extcurz_class object.
!>
!>  @param[inout] this A @ref extcurz_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE extcurz_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (extcurz_class), POINTER :: this

!  Start of executable code
      this%r = 0.0
      this%theta = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  Calculates the inclused Z current from Ampere's Law from an loop integration
!>  of the magnetic field. The loop integrated magnetic field is provided by
!>  @ref equilibrium_get_Int_B_dphi
!>
!>  @param[in]    this    A @ref extcurz_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION extcurz_get_modeled_signal(this, a_model, sigma)
      USE model
      USE stel_constants, only: mu0

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: extcurz_get_modeled_signal
      TYPE (extcurz_class), INTENT(in)        :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Start of executable code
      extcurz_get_modeled_signal = 0.0

      extcurz_get_modeled_signal(1) =                                          &
     &   equilibrium_get_Int_B_dphi(a_model%equilibrium, this%r,               &
     &                              this%theta)/mu0

      sigma = 0

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the extcurz type.
!>
!>  Returns a description of the extcurz type for use when writting output files.
!>
!>  @param[in] this A @ref extcurz_class instance.
!>  @returns A string describing the extcurz type.
!-------------------------------------------------------------------------------
      FUNCTION extcurz_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: extcurz_get_signal_type
      TYPE (extcurz_class), INTENT(in) :: this

!  Start of executable code
      extcurz_get_signal_type = 'extcur'

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
!>  @param[inout] this       A @ref extcurz_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE extcurz_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (extcurz_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                 :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                            :: error

!  Start of executable code
      CALL MPI_BCAST(this%r, 1, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%theta, 1, MPI_REAL8, 0, intra_comm, error)

#endif

      END SUBROUTINE

      END MODULE
