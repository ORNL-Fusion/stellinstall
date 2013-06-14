!*******************************************************************************
!  File v3fit_input.f
!  Contains module v3fit_input
!
!  This file contains all the variables and maximum sizes of the inputs for a
!  v3fit namelist input file.
!
!*******************************************************************************

      MODULE v3fit_input
      USE stel_kinds
      USE data_parameters
      USE signal
      USE pprofile_T

!*******************************************************************************
!  v3fit input module parameters
!*******************************************************************************

      INTEGER, PARAMETER :: v3fit_max_diagnostics = 1000
      INTEGER, PARAMETER :: v3fit_max_geometrics = 1000
      INTEGER, PARAMETER :: v3fit_max_combinations = 1000

      INTEGER, PARAMETER :: v3fit_max_signals = v3fit_max_diagnostics          &
     &                                        + v3fit_max_geometrics           &
     &                                        + v3fit_max_combinations

!  Maximum number of reconstruction parameters
      INTEGER, PARAMETER :: v3fit_max_parameters = 100

      INTEGER, PARAMETER :: v3fit_work_size = 100

      INTEGER, PARAMETER :: v3fit_max_spec_size = 100

!  Maximum number of phi positions for limiters
      INTEGER, PARAMETER :: v3fit_max_lif_size = 100

!  Maximum number of signals a combination signal can refer to.
      INTEGER, PARAMETER :: v3fit_max_cos_size = 100

!  Reconstruction step name length
      INTEGER, PARAMETER :: v3fit_step_name_length = 4

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) v3fit_namelist_class
!  2) v3fit_context_class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Derived type to hold all data from the namelist input file. Containing all
!  all this in a single place will allow memory to be freed after the
!  reconstruct is setup.
!-------------------------------------------------------------------------------
!  Namelist v3fit_main_nli definition
!
!          File name variables
!  main_nli_filename     File name for main namelist input
!  mdsig_list_filename   File name for list of MDSIG files
!  vmec_nli_filename     File name for VMEC namelist input
!  sxrch_dot_filename	 file holding sxr chord information
!  thscte_dot_filename	 file holding Thomson Scattering information
!  ipch_dot_filenam      file holding interferometry-polarimetry information 
!
!       Array allocation sizes, numbers to allocate
!  na_d_desc             d_desc array, diagnostic_desc
!  na_g_desc             g_desc array, geometric_desc
!  na_c_desc             c_desc array, coosig_desc
!  na_s_desc             s_desc array, signal_desc
!  na_rparam             rparam array, recon_param
!
!       Task specification, work variables
!  my_task               Character: to specify task in MAIN select case
!  i_work                Array of integers, work space
!  r_work                Array of reals, work space
!  c_work                Array of character variables, work space
!
!      Model profile specification
!  pp_ne_ptype           model electron density profile, parameterized profile type
!  pp_ne_b               array of b_coefficients for the electron density profile
!  pp_ne_as              array of as_coefficients electron density splines
!  pp_ne_af              array of af_coefficients electron density splines
!  pp_sxrem_ptype        model soft x-ray emissivity profile, parameterized profile type
!  pp_sxrem_b            array of b_coefficients for the soft x-ray emissivity profile
!  pp_sxrem_as           array of as_coefficients soft x-ray emissivity splines
!  pp_sxrem_af           array of af_coefficients soft x-ray emissivity splines
!  pp_te_ptype           model electron temperature profile, parameterized profile type
!  pp_te_b               array of b_coefficients for the e-temperature profile
!  pp_te_as              array of as_coefficients e-temperature splines
!  pp_te_af              array of af_coefficients e-temperature splines
!
!      Model Specification Variables
!  model_ne_type     Character variable to specify how the Soft x-ray
!                       emissivity will be calculated by the model.
!      'none'             Not part of the model.
!      'pp_ne'            Calculated from the parameterized profile pp_ne
!      'pp_te_vmec_p'     Calculated from the pp_te and the VMEC pressure
!  model_sxrem_type     Character variable to specify how the Soft x-ray
!                       emissivity will be calculated by the model.
!      'none'             Not part of the model.
!      'pp_sxrem'         Calculated from the parameterized profile pp_sxrem
!      'pp_ne_vmec_p'     Calculated from the pp_ne and the VMEC pressure
!  model_te_type        Character variable to specify how the e-temperature
!                       will be calculated by the model.
!      'none'             Not part of the model.
!      'pp_te'            Calculated from the parameterized profile pp_te
!      'pp_ne_vmec_p'     Calculated from the pp_ne and the VMEC pressure
!  ne_pp_unit           Scaling/conversion factor for the density.
!                       ne * ne_pp_unit is assumed to have units of m ** -3.
!  ne_min               Minimum electron density, m ** -3.
!  e_pressure_fraction  Electron pressure fraction of the total pressure
!                       Used when model_te_type = 'pp_ne_vmec_p'
!
!        Reconstruction constraints
!  n_rc                  Number of reconstruction constraints
!  rc_type               Array of reconstruction constraint types
!  rc_index              Array of reconstruction constraint indices
!  rc_value              Array of reconstruction constraint values
!
!       Reconstruction parameter specification
!  n_rp                  Number of reconstruction parameters
!  rp_type               Array of reconstruction parameter types
!  rp_index              Array of reconstruction parameter indices
!  rp_index2             Array of reconstruction parameter 2nd indices
!  rp_vrnc               Array of reconstruction parameter variances
!
!        Signal Data specification
!  n_sdata_o             Number of signal_data observations
!  iw_sdo_verbose        Integer to control write out of signal_data
!                           observations. (-1, no write, otherwise, verbosity)
!  sdo_data_a            Array of data for observed signals
!  sdo_sigma_a           Array of sigmas for observed signals
!  sdo_weight_a          Array of weight for observed signals
!
! These two are used by the extcurz diagnostic signal.
!  extcurz_s0            S position to integrate about.
!  extcurz_u0            U position to integrate about.
!
!        Signal Data sigma and weight specification
!  sdo_s_spec_imin       Array of minimum index values, for specifying sigmas
!  sdo_s_spec_imax       Array of maximum index values, for specifying sigmas
!  sdo_s_spec_floor      Array of floor values, for specifying sigmas
!  sdo_s_spec_fraction   Array of fractional values, for specifying sigmas
!      floor > 0, sigma = max(floor,fraction * sdo)
!      floor < 0, sigma = sqrt(floor ** 2 + (fraction * sdo) ** 2 )
!  sdo_w_spec_imin       Array of minimum index values, for specifying weights
!  sdo_w_spec_imax       Array of maximum index values, for specifying weights
!  sdo_w_spec_weight     Array of weight values, for specifying weights
!
!      Geometric information
!  r_major_radius        For circular torus geometric
!  a_minor_radius        For circular torus geometric
!
!     lif - Limiter_Iso Function - all declared in v3f_global
!       first dimension of arrays is dimensioned na_lif
!
!  n_lif                 number of limiter_iso functions
!  n_phi_lif             array of number of phi values (:)
!  lif_arz               array of r-z polynomial coefficients, (:,0:4,0:4)
!  lif_rc                array of r offset values (:)
!  lif_zc                array of z offset values (:)
!  lif_sigma             array of sigma values (:)
!  lif_phi_degree        array of phi values (:,na_phis_lif)
!  lif_on_edge           array of logical values (:)
!                           true  - edge of plasma is on limiter
!                           false - edge of plasma is within limiter
!
!  use_boundary          Toggle use of the fixed boundary coeffients as a
!                        reconstruction signal.
!
!     COOSIG - Combination Of Other SIGnals NLI Variables
!  n_coosig         Number of new signals of signal_type coosig
!  n_sig_coosig     Array of number of signals in a combinaton.
!  coosig_indices   2d Integer Array of indices (to list of signals)
!                   for combinations. First index is for which new signal,
!                   second index is for terms in the combination
!  coosig_coeff     2d real array of coefficients for the combination
!                   (a's in the type)
!  coosig_type      Character array of combination types
!                     'sum'  - linear combination
!                     'max'  - maximum of wieghted signals
!                     'min'  - minimum of weighted signals
!  coosig_name      Character array of coosig names
!  coosig_units     Character array of coosig units
!
!    Reconstruction step control NLI variables
!  nrstep         max number of reconstruction steps to perform
!  dg2_stop       Stopping criterion on change in g^2
!  cut_svd        cutoff value for relative Singular Values
!  cut_eff        cutoff value for expected step efficiency
!  cut_marg_eff   cutoff value for expected marginal step efficiency
!  cut_delta_a    cutoff value for expected step size
!  cut_dg2        cutoff value for expected change in g^2
!  astep_max      maximum allowable normalized step size
!  step_type      character specification of reconstruction step type
!                   'sl' - straight line
!                   'seg' - segmented
!                   'lm' - Levenberg-Marquardt
!
!    VMEC zero out xcdot control
!  l_zero_xcdot   reset the xcdot array to zero when converging vmec.
!-------------------------------------------------------------------------------

      CHARACTER (len=data_filename_length) :: main_nli_filename = ''

!  Equilibrium input files
      CHARACTER (len=data_filename_length) :: vmec_nli_filename = ''

!  Signal input files
      CHARACTER (len=data_filename_length) :: mdsig_list_filename = ''
      CHARACTER (len=data_filename_length) :: sxrch_dot_filename = ''
      CHARACTER (len=data_filename_length) :: ipch_dot_filename = ''
      CHARACTER (len=data_filename_length) :: thscte_dot_filename = ''

!  Allay allocation sizes. These should be set to the maxmium size.
      INTEGER :: na_d_desc = v3fit_max_diagnostics
      INTEGER :: na_g_desc = v3fit_max_geometrics
      INTEGER :: na_c_desc = v3fit_max_combinations
      INTEGER :: na_s_desc = v3fit_max_signals
      INTEGER :: na_rparam = v3fit_max_parameters

!  Task specification and work variables
      CHARACTER (len=data_name_length)         :: my_task = ''
!  FIXME: The following arrays are depricated and will be removed in the future.
      INTEGER, DIMENSION(v3fit_work_size)      :: i_work = 0
      REAL (rprec), DIMENSION(v3fit_work_size) :: r_work = 0.0
      CHARACTER (len=data_name_length), DIMENSION(v3fit_work_size) ::          &
     &   c_work = ''

!  Model profile specification
!  Electron Denisty
      CHARACTER (len=p_type_len)           :: pp_ne_ptype = 'none'
      REAL (rprec), DIMENSION(ilb_b:iub_b) :: pp_ne_b = 0.0
      REAL (rprec), DIMENSION(iub_asf)     :: pp_ne_as = 0.0
      REAL (rprec), DIMENSION(iub_asf)     :: pp_ne_af = 0.0

!  Soft x-ray emissvity
      CHARACTER (len=p_type_len)           :: pp_sxrem_ptype = 'none'
      REAL (rprec), DIMENSION(ilb_b:iub_b) :: pp_sxrem_b = 0.0
      REAL (rprec), DIMENSION(iub_asf)     :: pp_sxrem_as = 0.0
      REAL (rprec), DIMENSION(iub_asf)     :: pp_sxrem_af = 0.0

!  Electron temperature
      CHARACTER (len=p_type_len)           :: pp_te_ptype = 'none'
      REAL (rprec), DIMENSION(ilb_b:iub_b) :: pp_te_b = 0.0
      REAL (rprec), DIMENSION(iub_asf)     :: pp_te_as = 0.0
      REAL (rprec), DIMENSION(iub_asf)     :: pp_te_af = 0.0

!  Model Specification Variables
!  FIXME: New coding contains a model_ne_type that old code did not have. Make
!  the default value pp_ne for compatibility. Eventially this shoud be set back
!  to none and users should set model_ne_type in the namelist input.
      CHARACTER (len=data_name_length) :: model_ne_type = 'pp_ne'
      CHARACTER (len=data_name_length) :: model_sxrem_type = 'none'
      CHARACTER (len=data_name_length) :: model_te_type = 'none'
      REAL (rprec)                     :: ne_pp_unit = 0.0
      REAL (rprec)                     :: ne_min = 0.0
      REAL (rprec)                     :: e_pressure_fraction = 0.0

!  Reconstruction constraints
      INTEGER                            :: n_rc = 0
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_parameters) :: rc_type = ''
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters) :: rc_index = 0
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_parameters) :: rc_value = 0.0

!  Reconstruction parameters
      INTEGER                            :: n_rp = 0
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_parameters) :: rp_type = ''
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters) :: rp_index = 0
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters) :: rp_index2 = 0
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_parameters) :: rp_vrnc = 0.0

!  Signal Data specification
      INTEGER                         :: n_sdata_o = 0
      INTEGER                         :: iw_sdo_verbose = -1
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_signals) :: sdo_data_a = 0.0
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_signals) :: sdo_sigma_a = 0.0
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_signals) :: sdo_weight_a = 1.0

!  External current signal specification
!  FIXME: The label of s0 and u0 are vmec specific. These should be changed to a
!  radial coordinate and poloidal angle.
      REAL (rprec) :: extcurz_s0 = -1.0
      REAL (rprec) :: extcurz_u0 = 0.0

!  Signal data sigma specification
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   sdo_s_spec_imin = 0
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   sdo_s_spec_imax = 0
      REAL (rprec), DIMENSION(v3fit_max_spec_size) ::                          &
     &   sdo_s_spec_floor = 0.0
      REAL (rprec), DIMENSION(v3fit_max_spec_size) ::                          &
     &   sdo_s_spec_fraction = 0.0

!  Signal data weight specification
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   sdo_w_spec_imin = 0
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   sdo_w_spec_imax = 0
      REAL (rprec), DIMENSION(v3fit_max_spec_size) ::                          &
     &   sdo_w_spec_weight = 0.0

!  Geometric information
      REAL (rprec) :: r_major_radius
      REAL (rprec) :: a_minor_radius

!  Limiter iso function specification.
      INTEGER                                       :: n_lif = 0
      INTEGER, DIMENSION(v3fit_max_geometrics)      :: n_phi_lif = 0
      REAL (rprec), DIMENSION(v3fit_max_geometrics,0:4,0:4) ::                 &
     &   lif_arz = 0.0
      REAL (rprec), DIMENSION(v3fit_max_geometrics) ::                         &
     &   lif_rc = 0.0
      REAL (rprec), DIMENSION(v3fit_max_geometrics) ::                         &
     &   lif_zc = 0.0
      REAL (rprec), DIMENSION(v3fit_max_geometrics) ::                         &
     &   lif_sigma = 0.001_dp
      REAL (rprec), DIMENSION(v3fit_max_geometrics,v3fit_max_lif_size)         &
     &   :: lif_phi_degree = 0.0
      LOGICAL, DIMENSION(v3fit_max_geometrics)      ::                         &
     &   lif_on_edge = .false.

!  Use external boundary as signal
!  FIXME: The implentation is VMEC specific. This will eventually be removed and
!  replaced by priors when they are implemented.
      LOGICAL :: use_boundary

!  Combinations of other signals specification.
      INTEGER :: n_coosig = 0
      INTEGER, DIMENSION(v3fit_max_cos_size) :: n_sig_coosig = 0
      INTEGER, DIMENSION(v3fit_max_combinations,v3fit_max_cos_size) ::         &
     &   coosig_indices = 0
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_combinations,v3fit_max_cos_size) ::               &
     &   coosig_coeff = 0.0
      CHARACTER(len=combination_type_length),                                  &
     &   DIMENSION(v3fit_max_combinations) :: coosig_type = 'sum'
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_combinations) :: coosig_name = ''
      CHARACTER (len=signal_units_length),                                     &
     &   DIMENSION(v3fit_max_combinations) :: coosig_units = ''

!  Reconstruction step control specification.
      INTEGER                                :: nrstep = 0
      REAL (rprec)                           :: dg2_stop = 0.0
      REAL (rprec)                           :: cut_svd = 0.0
      REAL (rprec)                           :: cut_eff = 0.0
      REAL (rprec)                           :: cut_marg_eff = 0.0
      REAL (rprec)                           :: cut_delta_a = 0.0
      REAL (rprec)                           :: cut_dg2 = 0.0
      REAL (rprec)                           :: astep_max = 0.0
      CHARACTER (len=v3fit_step_name_length) :: step_type = 'sl'

!  Zero out vmec xcdot array
      LOGICAL :: l_zero_xcdot = .true.

!  Declare Namelist
      NAMELIST/v3fit_main_nli/                                                 &
     &   main_nli_filename,                                                    &
!  Equilibrium input files
     &   vmec_nli_filename,                                                    &
!  Signal input files
     &   mdsig_list_filename, sxrch_dot_filename, thscte_dot_filename,         &
     &   ipch_dot_filename,                                                    &
!  Allay allocation sizes
     &   na_d_desc, na_g_desc, na_c_desc, na_s_desc, na_rparam,                &
!  Task specification and work variables
     &   my_task,                                                              &
!  FIXME: The following arrays are depricated and will be removed in the future.
     &   i_work, c_work, r_work,                                               &
!  Model profile specification
!  Electron Denisty
     &   pp_ne_ptype, pp_ne_b, pp_ne_as, pp_ne_af,                             &
!  Soft x-ray emissvity
     &   pp_sxrem_ptype, pp_sxrem_b, pp_sxrem_as, pp_sxrem_af,                 &
!  Electron temperature
     &   pp_te_ptype, pp_te_b, pp_te_as, pp_te_af,                             &
!  Model Specification Variables
     &   model_ne_type, model_sxrem_type, model_te_type, ne_pp_unit,           &
     &   ne_min, e_pressure_fraction,                                          &
!  Reconstruction constraints
     &   n_rc, rc_type, rc_index, rc_value,                                    &
!  Reconstruction parameters
     &   n_rp, rp_type, rp_index, rp_index2, rp_vrnc,                          &
!  Signal Data specification
     &   n_sdata_o, iw_sdo_verbose, sdo_data_a, sdo_sigma_a,                   &
!  External current signal specification
!  FIXME: The label of s0 and u0 are vmec specific. These should be changed to a
!  radial coordinate and poloidal angle.
     &   extcurz_s0, extcurz_u0,                                               &
!  Signal data sigma specification
     &   sdo_s_spec_imin, sdo_s_spec_imax, sdo_s_spec_floor,                   &
     &   sdo_s_spec_fraction,                                                  &
!  Signal data weight specification
     &   sdo_w_spec_imin, sdo_w_spec_imax, sdo_w_spec_weight,                  &
!  Geometric information
     &   r_major_radius, a_minor_radius,                                       &
!  Limiter iso function specification.
     &   n_lif, n_phi_lif, lif_arz, lif_rc, lif_zc, lif_sigma,                 &
     &   lif_phi_degree, lif_on_edge,                                          &
!  Use external boundary as signal
!  FIXME: The implentation is VMEC specific. This will eventually be removed and
!  replaced by priors when they are implemented.
     &   use_boundary,                                                         &
!  Combinations of other signals specification.
     &   n_coosig, n_sig_coosig, coosig_indices, coosig_coeff,                 &
     &   coosig_type, coosig_name, coosig_units,                               &
!  Reconstruction step control specification.
     &   nrstep, dg2_stop, cut_svd, cut_eff, cut_marg_eff, cut_delta_a,        &
     &   cut_dg2, astep_max, step_type,                                        &
!  Zero out vmec xcdot array
     &   l_zero_xcdot

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE v3fit_input_read_namelist(namelist_file)
      USE safe_open_mod
      USE v3_utilities, only: err_fatal

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=data_filename_length), INTENT(in) :: namelist_file

!  local variables
      INTEGER :: iou_mnli, status, j, i, imin, imax

!  Start of executable code
      CALL safe_open(iou_mnli, status, TRIM(namelist_file),                    &
     &               'old', 'formatted')
      CALL assert_eq(0, status, 'v3fit_input_read_namelist' //                 &
     &   ': Safe_open of ' // TRIM(namelist_file) // ' failed')

!  Read the namelist input file.
      READ (iou_mnli, nml=v3fit_main_nli)
      CLOSE (iou_mnli, iostat=status)
      CALL assert_eq(0, status, 'v3fit_input_read_namelist' //                 &
     &   ': Error losing ' // TRIM(namelist_file) // ' failed')

!  Set weights from weight specification.
      DO j = 1, SIZE(sdo_w_spec_imin)
         IF (sdo_w_spec_imin(j) .le. 0) EXIT
         imin = MAX(MIN(sdo_w_spec_imin(j), n_sdata_o), 1)
         imax = MIN(MAX(sdo_w_spec_imax(j), imin), n_sdata_o)
         DO i = imin, imax
            sdo_weight_a(i) = sdo_w_spec_weight(j)
         END DO
      END DO

!  Set sigmas from the sigma specification.
      DO j = 1,SIZE(sdo_s_spec_imin)
         IF (sdo_s_spec_imin(j) .le. 0) EXIT
         imin = MAX(MIN(sdo_s_spec_imin(j), n_sdata_o), 1)
         imax = MIN(MAX(sdo_s_spec_imax(j), imin), n_sdata_o)
         DO i = imin, imax
            IF (sdo_sigma_a(i) .eq. 0) THEN
               IF (sdo_s_spec_floor(j) .ge. 0) THEN
                  sdo_sigma_a(i) = MAX(sdo_s_spec_floor(j),                    &
     &               ABS(sdo_s_spec_fraction(j)*sdo_data_a(i)))
               ELSE
                  sdo_sigma_a(i) = SQRT(sdo_s_spec_floor(j)**2.0 +             &
     &               (sdo_s_spec_fraction(j)*sdo_data_a(i))**2.0)
               END IF
            END IF
         END DO
      END DO

      IF (main_nli_filename .eq. '') THEN
         main_nli_filename = namelist_file
      END IF

      IF ((r_major_radius .gt. 0.0) .and.                                      &
     &    (a_minor_radius .gt. 0.0)) THEN
         IF (n_lif + 1 .le. v3fit_max_geometrics) THEN
            n_lif = n_lif + 1 ! Add this as a new limiter.
            lif_arz(n_lif,:,:) = 0.0
            lif_arz(n_lif,0,0) = -a_minor_radius**2.0
            lif_arz(n_lif,2,0) = 1.0
            lif_arz(n_lif,0,2) = 1.0
            lif_rc(n_lif) = r_major_radius
            lif_zc(n_lif) = 0.0
            lif_sigma(n_lif) = 0.001_rprec
            n_phi_lif(n_lif) = v3fit_max_lif_size
            DO i = 1, v3fit_max_lif_size
               lif_phi_degree(n_lif,i) = (i - 1)*360.0                         &
     &                                 / v3fit_max_lif_size
            END DO
            lif_on_edge(n_lif) = .true.
         ELSE
            CALL err_fatal('Adding limiter by specifying ' //                  &
     &                     'major/minor radii causes n_lif to ' //             &
     &                     'exceed v3fit_max_geometrics')
         END IF
      END IF

!  Search for the last index greater than 0 to determine the number of signals
!  in a combination.
!  FIXME: Remove this when n_sig_coosig is mandatory.
      DO i = 1, n_coosig
         DO j = 1, v3fit_max_cos_size
            IF (coosig_indices(i,j) .eq. 0) THEN
               n_sig_coosig(i) = j - 1
               EXIT ! Exit out of the j loop.
            END IF
         END DO
      END DO

      WRITE (*,*) ' *** V3FIT namelist input read from ' //                    &
     &            TRIM(namelist_file)

      END SUBROUTINE

      END MODULE
