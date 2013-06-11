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
      USE mddc_T, only: mddc_mrf, mddc_mrf_destroy, mddc_construct

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

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE magnetic_construct
      MODULE PROCEDURE magnetic_construct_netcdf,                              &
     &                 magnetic_construct_responce
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic with a responce function.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct_responce(responce_function)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER :: magnetic_construct_responce
      TYPE (mddc_mrf), POINTER, INTENT(in) :: responce_function

!  Start of executable code
      ALLOCATE(magnetic_construct_responce)

      magnetic_construct_responce%responce_function => responce_function

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a magnetic from a netcdf file .
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct_netcdf(mdsig_iou)
      USE data_parameters
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER  :: magnetic_construct_netcdf
      INTEGER, INTENT(in)             :: mdsig_iou

!  local variables
      TYPE (mddc_mrf), POINTER        :: responce_function

      CHARACTER (len=data_name_length) :: code_name
      CHARACTER (len=data_name_length) :: code_version
      CHARACTER (len=data_name_length) :: date_run
      CHARACTER (len=data_name_length) :: field_coils_id
      INTEGER                          :: n_field_cg
      REAL (rprec), DIMENSION(:), POINTER :: rdiag_coilg_1
      REAL (rprec), DIMENSION(:), POINTER :: extcur_mg
      INTEGER                          :: kp
      REAL (rprec)                     :: rmin
      REAL (rprec)                     :: rmax
      REAL (rprec)                     :: zmin
      REAL (rprec)                     :: zmax
      INTEGER                          :: n_field_periods
      LOGICAL                          :: lstell_sym
      INTEGER, DIMENSION(3)            :: dim_lengths
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_r
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_f
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_z

!  Start of executable code
      ALLOCATE(responce_function)

      CALL cdf_read(mdsig_iou,                                                 &
     &              'mddc_mrf_code_name', code_name)
      CALL cdf_read(mdsig_iou,                                                 &
     &              'mddc_mrf_code_version', code_version)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_date_run', date_run)
      CALL cdf_read(mdsig_iou,                                                 &
     &              'mddc_mrf_field_coils_id', field_coils_id)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_n_field_cg', n_field_cg)

      ALLOCATE(rdiag_coilg_1(n_field_cg))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rdiag_coilg_1', rdiag_coilg_1)

      ALLOCATE(extcur_mg(n_field_cg))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_extcur_mg', extcur_mg)

      CALL cdf_read(mdsig_iou, 'mddc_mrf_kp', kp)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rmin', rmin)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rmax', rmax)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_zmin', zmin)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_zmax', zmax)
      CALL cdf_read(mdsig_iou,                                                 &
     &              'mddc_mrf_n_field_periods', n_field_periods)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_lstell_sym', lstell_sym)

      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_r', dim_lengths)
      ALLOCATE(a_r(dim_lengths(1),dim_lengths(2),dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_r', a_r)

      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_f', dim_lengths)
      ALLOCATE(a_f(dim_lengths(1),dim_lengths(2),dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_f', a_f)

      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_z', dim_lengths)
      ALLOCATE(a_z(dim_lengths(1),dim_lengths(2),dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_z', a_z)

      CALL mddc_construct(responce_function, code_name, code_version,          &
     &                    date_run, field_coils_id, rdiag_coilg_1,             &
     &                    extcur_mg, kp, rmin, rmax, zmin, zmax,               &
     &                    n_field_periods, lstell_sym, a_r, a_f, a_z)

      magnetic_construct_netcdf => magnetic_construct(responce_function)

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
