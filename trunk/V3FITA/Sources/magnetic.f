!*******************************************************************************
!  File magnetic.f
!  Contains module magnetic
!  Defines the base class of the type magnetic:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE magnetic

      USE stel_kinds, only: rprec
      USE mddc_T, only: mddc_mrf, mddc_mrf_destroy

      IMPLICIT NONE

!*******************************************************************************
!  magnetic module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: magnetic_no_type = -1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) magnetic base class
!
!*******************************************************************************

      TYPE magnetic_class
         TYPE (mddc_mrf), POINTER :: responce_function => null()
      END TYPE magnetic_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct(responce_function)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER       :: magnetic_construct
      TYPE (mddc_mrf), POINTER, INTENT(in) :: responce_function

!  Start of executable code
      ALLOCATE(magnetic_construct)

      magnetic_construct%responce_function => responce_function

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE magnetic_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%responce_function)) THEN
         CALL mddc_mrf_destroy(this%responce_function)
         this%responce_function => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION magnetic_get_modeled_signal(this, a_model, sigma)
      USE model
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: magnetic_get_modeled_signal
      TYPE (magnetic_class), INTENT(in)       :: this
      TYPE (model_class), INTENT(in)          :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Start of executable code
!  Calculate the plasma only signal from the responce function.
      magnetic_get_modeled_signal(2) =                                         &
     &   equilibrium_get_plasma_mag_signal(a_model%equilibrium,                &
     &      this%responce_function%a_r,                                        &
     &      this%responce_function%a_f,                                        &
     &      this%responce_function%a_z,                                        &
     &      this%responce_function%rmax,                                       &
     &      this%responce_function%rmin,                                       &
     &      this%responce_function%ir,                                         &
     &      this%responce_function%zmax,                                       &
     &      this%responce_function%zmin,                                       &
     &      this%responce_function%jz,                                         &
     &      this%responce_function%lstell_sym)

!  Calculate the induced signal from external currents.
      magnetic_get_modeled_signal(3) =                                         &
     &   equilibrium_get_ext_mag_signal(a_model%equilibrium,                   &
     &      this%responce_function%rdiag_coilg_1,                              &
     &      this%responce_function%extcur_mg,                                  &
     &      this%responce_function%n_field_cg)

!  Add plasma signals and induced signal to create total signal.
      magnetic_get_modeled_signal(1) =                                         &
     &   SUM(magnetic_get_modeled_signal(2:3))

      sigma = 0.0

      END FUNCTION

      END MODULE
