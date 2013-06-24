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
!  Construct an signal equilibrium with a vmec equilibrium.
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
!-------------------------------------------------------------------------------
!  Set reconstruction parameter
!-------------------------------------------------------------------------------
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

!  Set param needs to come first so that function gets executed. Otherwise, when
!  this%not_converged was true, the expression will evaluate to false all the
!  time and not bother calling vmec_set_param.
         CASE (equilibrium_vmec_type)
            this%not_converged = vmec_set_param(this%vmec, id,                 &
     &                                          i_index, j_index,              &
     &                                          value)                         &
     &                         .or. this%not_converged

      END SELECT

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Set magnetic cache
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_set_magnetic_cache(this, total_tor_planes,        &
     &                                          tor_planes_used)

      IMPLICIT NONE

!  Declare Sruments
      TYPE (equilibrium_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: total_tor_planes
      INTEGER, INTENT(in)                     :: tor_planes_used

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_set_magnetic_cache(this%vmec, total_tor_planes,          &
     &                                   tor_planes_used)

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
      TYPE (equilibrium_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)        :: param_name

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
!  Get the value of a parameter for an id number
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_param_value(this, id, i_index, j_index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_param_value
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id
      INTEGER, INTENT(in)                  :: i_index
      INTEGER, INTENT(in)                  :: j_index

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_param_value =                                      &
     &         vmec_get_param_value(this%vmec, id, i_index, j_index)

         CASE DEFAULT
            equilibrium_get_param_value = 0.0

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
!  Get volume magnetic volume integration radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_rgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_rgrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_rgrid =>                           &
     &         vmec_get_magnetic_volume_rgrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_rgrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get volume magnetic volume integration z grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_zgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_zgrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_zgrid =>                           &
     &         vmec_get_magnetic_volume_zgrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_zgrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get volume magnetic volume integration jr grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_jrgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_jrgrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_jrgrid =>                          &
     &         vmec_get_magnetic_volume_jrgrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_jrgrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get volume magnetic volume integration jphi grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_jphigrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_jphigrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_jphigrid =>                        &
     &         vmec_get_magnetic_volume_jphigrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_jphigrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get volume magnetic volume integration jz grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_jzgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_jzgrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_jzgrid =>                          &
     &         vmec_get_magnetic_volume_jzgrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_jzgrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get volume integration element.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_volume_int_element(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_volume_int_element
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

        CASE (equilibrium_vmec_type)
            equilibrium_get_volume_int_element =                               &
     &         vmec_get_volume_int_element(this%vmec)

         CASE DEFAULT
            equilibrium_get_volume_int_element = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get external current.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_ext_currents(this, scale_currents)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:), POINTER  ::                                  &
     &   equilibrium_get_ext_currents
      TYPE (equilibrium_class), INTENT(in) :: this
      LOGICAL, INTENT(out)                 :: scale_currents

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_ext_currents =>                                    &
     &         vmec_get_ext_currents(this%vmec, scale_currents)

         CASE DEFAULT
            scale_currents = .false.
            equilibrium_get_ext_currents => null()

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
!-------------------------------------------------------------------------------
!  Converge the equilibrium.
!-------------------------------------------------------------------------------
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

      WRITE (*,*) 'Equilibrium converged ', equilibrium_converge,              &
     &            num_iter

      END FUNCTION

!-------------------------------------------------------------------------------
!  Save the internal state of the equilibrium.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_save_state(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(inout) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_save_state(this%vmec)

      END SELECT

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Reset the internal state of the equilibrium.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_reset_state(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(inout) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_reset_state(this%vmec)

      END SELECT

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Write the equilibrium out to disk.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_write(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: iou

!  Start of executable code
      SELECT CASE (this%type)

      CASE (equilibrium_vmec_type)
         CALL vmec_write(this%vmec, iou)

      END SELECT

      END SUBROUTINE

      END MODULE
