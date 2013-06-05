!*******************************************************************************
!  File equilibrium.f
!  Contains module equilibrium
!  Defines the base class of the type equilibrium:
!
!  Super Class: None
!  Sub Classes: vmec
!
!*******************************************************************************

      MODULE equilibrium
      USE vmec_equilibrium

      IMPLICIT NONE

!*******************************************************************************
!  equilibrium module parameters
!*******************************************************************************

      INTEGER, PARAMETER :: equilibrium_no_type = -1
      INTEGER, PARAMETER :: equilibrium_vmec_type = 0

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) equilibrium base class
!
!*******************************************************************************

      TYPE equilibrium_class
         INTEGER                    :: type = equilibrium_no_type
         LOGICAL                    :: not_converged
         TYPE (vmec_class), POINTER :: vmec => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE equilibrium_construct
         MODULE PROCEDURE equilibrium_construct_vmec
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct an signal equilibrium a vmec equilibrium.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_construct_vmec(vmec_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), POINTER :: equilibrium_construct_vmec
      TYPE (vmec_class), POINTER, INTENT(in) :: vmec_object

!  Start of executable code
      ALLOCATE(equilibrium_construct_vmec)

      equilibrium_construct_vmec%type = equilibrium_vmec_type
      equilibrium_construct_vmec%vmec => vmec_object
      equilibrium_construct_vmec%not_converged = .true.

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE equilibrium_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%vmec)) THEN
         CALL vmec_destruct(this%vmec)
         this%vmec => null()
      END IF

      this%type = equilibrium_no_type
      this%not_converged = .true.

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
      SUBROUTINE equilibrium_set_param(this, id, i_index, j_index,             &
     &                                 value)

      IMPLICIT NONE

!  Declare Sruments
      TYPE (equilibrium_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: id
      INTEGER, INTENT(in)                     :: i_index
      INTEGER, INTENT(in)                     :: j_index
      REAL (rprec), INTENT(in)                :: value

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            this%not_converged = this%not_converged .or.                       &
     &         vmec_set_param(this%vmec, id, i_index, j_index, value)

      END SELECT

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the id number for a equilibrium parameter
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_param_id(this, param_name)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: equilibrium_get_param_id
      TYPE (equilibrium_class), INTENT(in)         :: this
      CHARACTER (len=data_name_length), INTENT(in) :: param_name

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_param_id =                                         &
     &         vmec_get_param_id(this%vmec, param_name)

         CASE DEFAULT
            equilibrium_get_param_id = data_no_id

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!   Get Scaler Real value for an id number
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_scaler_value(this, id)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_scaler_value
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_scaler_value =                                     &
     &         vmec_get_scaler_value(this%vmec, id)

         CASE DEFAULT
            equilibrium_get_scaler_value = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get equilibrium parameter name for an id number
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_param_name(this, id)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=data_name_length) :: equilibrium_get_param_name
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_param_name =                                       &
     &         vmec_get_param_name(this%vmec, id)

         CASE DEFAULT
            equilibrium_get_param_name = ''

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get equilibrium electron density at position
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_ne(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_ne
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_ne = vmec_get_ne(this%vmec, x_cart)

         CASE DEFAULT
            equilibrium_get_ne = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get equilibrium electron temperature at position
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_te(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_te
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_te = vmec_get_te(this%vmec, x_cart)

         CASE DEFAULT
            equilibrium_get_te = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get equilibrium soft x-ray emissivity at position
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_sxrem(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_sxrem
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_sxrem = vmec_get_sxrem(this%vmec, x_cart)

         CASE DEFAULT
            equilibrium_get_sxrem = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get equilibrium plasma pressure at position
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_p(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_p
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_p = vmec_get_p(this%vmec, x_cart)

         CASE DEFAULT
            equilibrium_get_p = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get equilibrium b field at position
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_B_vec(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3) :: equilibrium_get_B_vec
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_B_vec = vmec_get_B_vec(this%vmec, x_cart)

         CASE DEFAULT
            equilibrium_get_B_vec = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get equilibrium line integrated b at position
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_Int_B_dphi(this, r, theta)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_Int_B_dphi
      TYPE (equilibrium_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)             :: r
      REAL (rprec), INTENT(in)             :: theta

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_Int_B_dphi =                                       &
     &         vmec_get_Int_B_dphi(this%vmec, r, theta)

         CASE DEFAULT
            equilibrium_get_Int_B_dphi = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get plasma edge at phi position
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_plasma_edge(this, phi, r, z)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: equilibrium_get_plasma_edge
      TYPE (equilibrium_class), INTENT(in)             :: this
      REAL (rprec), INTENT (in)                        :: phi
      REAL (rprec), DIMENSION(:), POINTER, INTENT(out) :: r
      REAL (rprec), DIMENSION(:), POINTER, INTENT(out) :: z

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_plasma_edge =                                      &
     &         vmec_get_plasma_edge(this%vmec, phi, r, z)

         CASE DEFAULT
            equilibrium_get_plasma_edge = 0
            r => null()
            z => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get plasma only magnetic signal.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_plasma_mag_signal(this,                         &
     &                                           r_fun, phi_fun, z_fun,        &
     &                                           r_max, r_min, num_r,          &
     &                                           z_max, z_min, num_z,          &
     &                                           is_stell_sym)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_plasma_mag_signal
      TYPE (equilibrium_class), INTENT(in)       :: this
      REAL (rprec), DIMENSION(:,:,:), INTENT(in) :: r_fun
      REAL (rprec), DIMENSION(:,:,:), INTENT(in) :: phi_fun
      REAL (rprec), DIMENSION(:,:,:), INTENT(in) :: z_fun
      REAL (rprec), INTENT(in)                   :: r_max
      REAL (rprec), INTENT(in)                   :: r_min
      INTEGER, INTENT(in)                        :: num_r
      REAL (rprec), INTENT(in)                   :: z_max
      REAL (rprec), INTENT(in)                   :: z_min
      INTEGER, INTENT(in)                        :: num_z
      LOGICAL, INTENT(in)                        :: is_stell_sym

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_plasma_mag_signal =                                &
     &         vmec_get_plasma_mag_signal(this%vmec,                           &
     &                                    r_fun, phi_fun, z_fun,               &
     &                                    r_max, r_min, num_r,                 &
     &                                    z_max, z_min, num_z,                 &
     &                                    is_stell_sym)

         CASE DEFAULT
            equilibrium_get_plasma_mag_signal = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get external current magnetic signal.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_ext_mag_signal(this, diag_coils,                &
     &                                        current_scales,                  &
     &                                        num_coil_groups)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_ext_mag_signal
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(:), INTENT(in) :: diag_coils
      REAL (rprec), DIMENSION(:), INTENT(in) :: current_scales
      INTEGER, INTENT(in)                    :: num_coil_groups

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_ext_mag_signal =                                   &
     &         vmec_get_ext_mag_signal(this%vmec, diag_coils,                  &
     &                                 current_scales,                         &
     &                                 num_coil_groups)

         CASE DEFAULT
            equilibrium_get_ext_mag_signal = 0.0

      END SELECT

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Query variable id for scaler type
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_scaler_value(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_scaler_value
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_scaler_value =                                      &
     &         vmec_is_scaler_value(this%vmec, id)

         CASE DEFAULT
            equilibrium_is_scaler_value = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 1D array
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_1d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_1d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_1d_array = vmec_is_1d_array(this%vmec, id)

         CASE DEFAULT
            equilibrium_is_1d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 2D array
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_2d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_2d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_2d_array = vmec_is_2d_array(this%vmec, id)

         CASE DEFAULT
            equilibrium_is_2d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 3D array
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_3d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_3d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_3d_array = vmec_is_3d_array(this%vmec, id)

         CASE DEFAULT
            equilibrium_is_3d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 4D array
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_4d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_4d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_4d_array = vmec_is_4d_array(this%vmec, id)

         CASE DEFAULT
            equilibrium_is_4d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id is a valid reconstruction parameter
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_recon_param(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_recon_param
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_recon_param =                                       &
     &         vmec_is_recon_param(this%vmec, id)

         CASE DEFAULT
            equilibrium_is_recon_param = .false.

      END SELECT

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
      FUNCTION equilibrium_converge(this, num_iter)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                 :: equilibrium_converge
      TYPE (equilibrium_class), INTENT(inout) :: this
      INTEGER, INTENT(inout)                  :: num_iter

!  Start of executable code
      IF (this%not_converged) THEN
         SELECT CASE (this%type)

            CASE (equilibrium_vmec_type)
               equilibrium_converge = vmec_converge(this%vmec, num_iter)
               this%not_converged = .not.equilibrium_converge

            CASE DEFAULT
               equilibrium_converge = .false.
               this%not_converged = .true.

         END SELECT
      ELSE
         equilibrium_converge = .true.
      ENDIF

      END FUNCTION

      END MODULE
