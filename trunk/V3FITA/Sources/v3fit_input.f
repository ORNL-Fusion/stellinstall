!-------------------------------------------------------------------------------
!>  @page namelist_sec Namelist v3fit_main_nli definition
!>
!>  The V3FIT namelist input file is defined in the @ref v3fit_input module.
!>
!>  @section filename_sec Filename variables
!>  * main_nli_filename     File name for main namelist input
!>                          @ref v3fit_input::main_nli_filename
!>  * vmec_nli_filename     File name for VMEC namelist input
!>                          @ref v3fit_input::vmec_nli_filename
!>  * mdsig_list_filename   File name for list of MDSIG files
!>                          @ref v3fit_input::mdsig_list_filename
!>  * sxrch_dot_filename	file holding sxr chord information
!>                          @ref v3fit_input::sxrch_dot_filename
!>  * thscte_dot_filename	file holding Thomson Scattering information
!>                          @ref v3fit_input::thscte_dot_filename
!>  * ipch_dot_filename     file holding interferometry-polarimetry
!>                          information
!>                          @ref v3fit_input::ipch_dot_filename
!>
!>  @section size_sec Array allocation sizes, numbers to allocate
!>  * na_d_desc             d_desc array, diagnostic_desc (deprecated)
!>                          @ref v3fit_input::na_d_desc
!>  * na_g_desc             g_desc array, geometric_desc (deprecated)
!>                          @ref v3fit_input::na_g_desc
!>  * na_c_desc             c_desc array, coosig_desc (deprecated)
!>                          @ref v3fit_input::na_c_desc
!>  * na_s_desc             s_desc array, signal_desc
!>                          @ref v3fit_input::na_s_desc
!>  * na_rparam             rparam array, recon_param (deprecated)
!>                          @ref v3fit_input::na_rparam
!>
!>  @section work_sec Task specification, work variables
!>  * my_task               Character: to specify task in MAIN select case
!>                          @ref v3fit_input::my_task
!>  * i_work                Array of integers, work space (deprecated)
!>                          @ref v3fit_input::i_work
!>  * r_work                Array of reals, work space (deprecated)
!>                          @ref v3fit_input::r_work
!>  * c_work                Array of character variables, work space
!>                          (deprecated)
!>                          @ref v3fit_input::c_work
!>
!>  @section model_sec Model profile specification
!>  * pp_ne_ptype           model electron density profile, parameterized
!>                          profile type
!>                          @ref v3fit_input::pp_ne_ptype
!>  * pp_ne_b               array of b_coefficients for the electron density
!>                          profile
!>                          @ref v3fit_input::pp_ne_b
!>  * pp_ne_as              array of as_coefficients electron density splines
!>                          @ref v3fit_input::pp_ne_as
!>  * pp_ne_af              array of af_coefficients electron density splines
!>                          @ref v3fit_input::pp_ne_af
!>  * pp_sxrem_ptype        model soft x-ray emissivity profile, parameterized
!>                          profile type
!>                          @ref v3fit_input::pp_sxrem_ptype
!>  * pp_sxrem_b            array of b_coefficients for the soft x-ray
!>                          emissivity profile
!>                          @ref v3fit_input::pp_sxrem_b
!>  * pp_sxrem_as           array of as_coefficients soft x-ray emissivity
!>                          splines
!>                          @ref v3fit_input::pp_sxrem_as
!>  * pp_sxrem_af           array of af_coefficients soft x-ray emissivity
!>                          splines
!>                          @ref v3fit_input::pp_sxrem_af
!>  * pp_te_ptype           model electron temperature profile, parameterized
!>                          profile type
!>                          @ref v3fit_input::pp_te_ptype
!>  * pp_te_b               array of b_coefficients for the e-temperature
!>                          profile
!>                          @ref v3fit_input::pp_te_b
!>  * pp_te_as              array of as_coefficients e-temperature splines
!>                          @ref v3fit_input::pp_te_as
!>  * pp_te_af              array of af_coefficients e-temperature splines
!>                          @ref v3fit_input::pp_te_af
!>
!>  @section model_spec_sec Model Specification Variables
!>  * model_ne_type         Character variable to specify how the electron
!>                          density will be calculated by the model.
!>                          @ref v3fit_input::model_ne_type
!>      -# 'none'             Not part of the model.
!>      -# 'pp_ne'            Calculated from the parameterized profile pp_ne
!>      -# 'pp_te_vmec_p'     Calculated from the pp_te and the VMEC pressure
!>  * model_te_type      Character variable to specify how the Soft x-ray
!>                          emissivity will be calculated by the model.
!>                          @ref v3fit_input::model_te_type
!>      -# 'none'             Not part of the model.
!>      -# 'pp_sxrem'         Calculated from the parameterized profile pp_sxrem
!>  * model_te_type         Character variable to specify how the
!>                          e-temperature will be calculated by the model.
!>                          @ref v3fit_input::model_te_type
!>      -# 'none'             Not part of the model.
!>      -# 'pp_te'            Calculated from the parameterized profile pp_te
!>      -# 'pp_ne_vmec_p'     Calculated from the pp_ne and the VMEC pressure
!>  * ne_pp_unit            Scaling/conversion factor for the density.
!>                          ne * ne_pp_unit is assumed to have units of
!>                          m^-3.
!>                          @ref v3fit_input::ne_pp_unit
!>  * ne_min                Minimum electron density, m^-3.
!>                          @ref v3fit_input::ne_min
!>  * e_pressure_fraction   Electron pressure fraction of the total
!>                          pressure. Used when model_te_type = 'pp_ne_vmec_p'
!>                          or model_ne_type = 'pp_te_vmec_p'
!>                          @ref v3fit_input::e_pressure_fraction
!>
!>  @section recon_con_sec Reconstruction constraints
!>  * n_rc                  Number of reconstruction constraints
!>                          @ref v3fit_input::n_rc
!>  * rc_type               Array of reconstruction constraint types
!>                          @ref v3fit_input::rc_type
!>  * rc_index              Array of reconstruction constraint indices
!>                          @ref v3fit_input::rc_index
!>  * rc_value              Array of reconstruction constraint values
!>                          @ref v3fit_input::rc_value
!>
!>  @section recon_param_sec Reconstruction parameter specification
!>  * n_rp                  Number of reconstruction parameters
!>                          @ref v3fit_input::n_rp
!>  * rp_type               Array of reconstruction parameter types
!>                          @ref v3fit_input::rp_type
!>  * rp_index              Array of reconstruction parameter indices
!>                          @ref v3fit_input::rp_index
!>  * rp_index2             Array of reconstruction parameter 2nd indices
!>                          @ref v3fit_input::rp_index2
!>  * rp_vrnc               Array of reconstruction parameter variances
!>                          @ref v3fit_input::rp_vrnc
!>  * rp_range_type         2D Array of reconstruction parameter range types
!>                          @ref v3fit_input::rp_range_type
!>      -# 'infinity'          Unbounded
!>      -# 'value'             Bounded by the value of rp_range_value
!>      -# \<other\>           Bounded by a equilibrium parameter.
!>  * rp_range_value        2D Array of values used if rp_range_type = 'value'
!>                          @ref v3fit_input::rp_range_value
!>  * rp_range_index        2D Array of values of indicies used if
!>                          rp_range_type = \<other\>
!>                          @ref v3fit_input::rp_range_index
!>
!>  @section signal_spec_sec Signal Data specification
!>  * n_sdata_o             Number of signal_data observations
!>                          @ref v3fit_input::n_sdata_o
!>  * iw_sdo_verbose        Integer to control write out of signal_data
!>                          observations. (-1, no write, otherwise, verbosity)
!>                          @ref v3fit_input::iw_sdo_verbose
!>  * sdo_data_a            Array of data for observed signals
!>                          @ref v3fit_input::sdo_data_a
!>  * sdo_sigma_a           Array of sigmas for observed signals
!>                          @ref v3fit_input::sdo_sigma_a
!>  * sdo_weight_a          Array of weight for observed signals
!>                          @ref v3fit_input::sdo_weight_a
!>        These two are used by the extcurz diagnostic signal.
!>  * extcurz_s0            S position to integrate about.
!>                          @ref v3fit_input::extcurz_s0
!>  * extcurz_u0            U position to integrate about.
!>                          @ref v3fit_input::extcurz_u0
!>        Signal Data sigma and weight specification
!>  * sdo_s_spec_imin       Array of minimum index values, for specifying
!>                          sigmas
!>                          @ref v3fit_input::sdo_s_spec_imin
!>  * sdo_s_spec_imax       Array of maximum index values, for specifying
!>                          sigmas
!>                          @ref v3fit_input::sdo_s_spec_imax
!>  * sdo_s_spec_floor      Array of floor values, for specifying sigmas
!>                          @ref v3fit_input::sdo_s_spec_floor
!>  * sdo_s_spec_fraction   Array of fractional values, for specifying sigmas
!>                          @ref v3fit_input::sdo_s_spec_fraction
!>      - floor \> 0, sigma = max(floor,fraction * sdo)
!>      - floor \< 0, sigma = sqrt(floor ** 2 + (fraction * sdo) ** 2)
!>  * sdo_w_spec_imin       Array of minimum index values, for specifying
!>                          weights
!>                          @ref v3fit_input::sdo_w_spec_imin
!>  * sdo_w_spec_imax       Array of maximum index values, for specifying
!>                          weights
!>                          @ref v3fit_input::sdo_w_spec_imax
!>  * sdo_w_spec_weight     Array of weight values, for specifying weights
!>                          @ref v3fit_input::sdo_w_spec_weight
!>
!>  @section geometric_sec Geometric information
!>  * r_major_radius        For circular torus geometric
!>                          @ref v3fit_input::r_major_radius
!>  * a_minor_radius        For circular torus geometric
!>                          @ref v3fit_input::a_minor_radius
!>  @subsection lif_sec Limiter_Iso Function.
!>
!>  First dimension of arrays is dimensioned na_lif
!>  * n_lif                 number of limiter_iso functions
!>                          @ref v3fit_input::n_lif
!>  * n_phi_lif             array of number of phi values (:)
!>                          @ref v3fit_input::n_phi_lif
!>  * lif_arz               array of r-z polynomial coefficients, (:,0:4,0:4)
!>                          @ref v3fit_input::lif_arz
!>  * lif_rc                array of r offset values (:)
!>                          @ref v3fit_input::lif_rc
!>  * lif_zc                array of z offset values (:)
!>                          @ref v3fit_input::lif_zc
!>  * lif_sigma             array of sigma values (:)
!>                          @ref v3fit_input::lif_sigma
!>  * lif_phi_degree        array of phi values (:,na_phis_lif)
!>                          @ref v3fit_input::lif_phi_degree
!>  * lif_on_edge           array of logical values (:)
!>                          @ref v3fit_input::lif_on_edge
!>      - true  - edge of plasma is on limiter
!>      - false - edge of plasma is within limiter
!>
!>  * use_boundary          Toggle use of the fixed boundary coeffients as a
!>                          reconstruction signal. (deprecated)
!>                          @ref v3fit_input::lif_on_edge
!>
!>  @section coosig_sec COOSIG - Combination Of Other SIGnals NLI Variables
!>  * n_coosig              Number of new signals of signal_type coosig
!>                          @ref v3fit_input::n_coosig
!>  * n_sig_coosig          Array of number of signals in a combinaton.
!>                          @ref v3fit_input::n_sig_coosig
!>  * coosig_indices        2d Integer Array of indices (to list of signals)
!>                          for combinations. First index is for which new
!>                          signal, second index is for terms in the combination
!>                          @ref v3fit_input::coosig_indices
!>  * coosig_coeff          2d real array of coefficients for the combination
!>                          (a's in the type)
!>                          @ref v3fit_input::coosig_coeff
!>  * coosig_type           Character array of combination types
!>                          @ref v3fit_input::coosig_type
!>      -# 'sum'  - linear combination
!>      -# 'max'  - maximum of wieghted signals
!>      -# 'min'  - minimum of weighted signals
!>  * coosig_name           Character array of coosig names
!>                          @ref v3fit_input::coosig_name
!>  * coosig_units          Character array of coosig units
!>                          @ref v3fit_input::coosig_units
!>
!>  @section recon_sec Reconstruction step control NLI variables
!>  * nrstep                max number of reconstruction steps to perform
!>                          @ref v3fit_input::nrstep
!>  * dg2_stop              Stopping criterion on change in g^2
!>                          @ref v3fit_input::dg2_stop
!>  * cut_svd               cutoff value for relative Singular Values
!>                          @ref v3fit_input::cut_svd
!>  * cut_eff               cutoff value for expected step efficiency
!>                          @ref v3fit_input::cut_eff
!>  * cut_marg_eff          cutoff value for expected marginal step efficiency
!>                          @ref v3fit_input::coosig_units
!>  * cut_delta_a           cutoff value for expected step size
!>                          @ref v3fit_input::cut_delta_a
!>  * cut_dg2               cutoff value for expected change in g^2
!>                          @ref v3fit_input::cut_dg2
!>  * astep_max             maximum allowable normalized step size
!>                          @ref v3fit_input::astep_max
!>  * step_type             character specification of reconstruction step type
!>                          @ref v3fit_input::step_type
!>      -# 'sl'  - straight line
!>      -# 'seg' - segmented
!>      -# 'lm'  - Levenberg-Marquardt
!>
!>  @section xcdot_sec VMEC zero out xcdot control
!>  * l_zero_xcdot          reset the xcdot array to zero when converging vmec.
!>                          @ref v3fit_input::l_zero_xcdot
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file v3fit_input.f
!>  @brief Contains module @ref v3fit_input
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains all the variables and maximum sizes of the inputs for a
!>  v3fit namelist input file. The module contained within does not represent an
!>  object instance. Instead all variables are contained in a global context.
!>  This is required due to limitations of FORTRAN 95 and namelist inputs. Any
!>  information needed by a V3FIT task should be copied to a @ref v3fit_context
!>  object. All non-parameters are inputs to the namelist input file.
!>  
!>  @ref namelist_sec "Namelist v3fit_main_nli definition"
!>
!>  @note Some of the references are missing here. This is due to a bug in
!>  Doxygen when variable decalarations span multiple lines.
!*******************************************************************************
      MODULE v3fit_input
      USE stel_kinds
      USE data_parameters
      USE signal
      USE pprofile_T

!*******************************************************************************
!  v3fit input module parameters
!*******************************************************************************
!>  Maximum number of diagnostic signal.
      INTEGER, PARAMETER :: v3fit_max_diagnostics = 1000
!>  Maximum number of geometric signal.
      INTEGER, PARAMETER :: v3fit_max_geometrics = 1000
!>  Maximum number of combination signal.
      INTEGER, PARAMETER :: v3fit_max_combinations = 1000
!>  Maximum number of total signal.
      INTEGER, PARAMETER :: v3fit_max_signals = v3fit_max_diagnostics          &
     &                                        + v3fit_max_geometrics           &
     &                                        + v3fit_max_combinations

!>  Maximum number of reconstruction parameters.
      INTEGER, PARAMETER :: v3fit_max_parameters = 100

!>  @deprecated
!>  Maximum size of the work arrays.
      INTEGER, PARAMETER :: v3fit_work_size = 100

!>  Maximum size of the easy specification arrays.
!>  @see sdo_s_spec_imin
!>  @see sdo_s_spec_imax
!>  @see sdo_s_spec_floor
!>  @see sdo_s_spec_fraction
!>  @see sdo_w_spec_imin
!>  @see sdo_w_spec_imax
!>  @see sdo_w_spec_weight
      INTEGER, PARAMETER :: v3fit_max_spec_size = 100

!>  Maximum number of phi positions for limiters.
      INTEGER, PARAMETER :: v3fit_max_lif_size = 100

!>  Maximum number of signals a combination signal can refer to.
      INTEGER, PARAMETER :: v3fit_max_cos_size = 100

!>  Reconstruction step name length.
      INTEGER, PARAMETER :: v3fit_step_name_length = 4

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) v3fit_namelist_class
!  2) v3fit_context_class
!
!*******************************************************************************
!>  File name for main namelist input
      CHARACTER (len=data_filename_length) :: main_nli_filename = ''

!  Equilibrium input files
!>  File name for VMEC namelist input
      CHARACTER (len=data_filename_length) :: vmec_nli_filename = ''

!  Signal input files
!>  File name for list of MDSIG files
      CHARACTER (len=data_filename_length) :: mdsig_list_filename = ''
!>  File holding sxr chord information
      CHARACTER (len=data_filename_length) :: sxrch_dot_filename = ''
!>  File holding Thomson Scattering information
      CHARACTER (len=data_filename_length) :: ipch_dot_filename = ''
!>  File holding interferometry-polarimetry information
      CHARACTER (len=data_filename_length) :: thscte_dot_filename = ''

!  Array allocation sizes. These should be set to the maxmium size.
!>  @deprecated
!>  Maximum number of diagnostic to create.
      INTEGER :: na_d_desc = v3fit_max_diagnostics
!>  @deprecated
!>  Maximum number of geometrics to create.
      INTEGER :: na_g_desc = v3fit_max_geometrics
!>  @deprecated
!>  Maximum number of combinations to create.
      INTEGER :: na_c_desc = v3fit_max_combinations
!>  Maximum number of signals to create.
!>  @note this needs to be greater than or equal to total number of signals
!>  used. The signal arrays are resized to free up memory if this value is
!>  larger than the number of signals constructed. If this value is equal to the
!>  number of signals, the resize step can be avoided.
      INTEGER :: na_s_desc = v3fit_max_signals
!>  @deprecated
!>  Maximum number of reconstruction parameters.
      INTEGER :: na_rparam = v3fit_max_parameters

!  Task specification and work variables
!>  The v3fit task.
!>  @see v3fit
      CHARACTER (len=data_name_length)         :: my_task = ''
!>  @deprecated
!>  Work array for integer inputs.
      INTEGER, DIMENSION(v3fit_work_size)      :: i_work = 0
!>  @deprecated
!>  Work array for real inputs.
      REAL (rprec), DIMENSION(v3fit_work_size) :: r_work = 0.0
!>  @deprecated
!>  Work array for string inputs.
      CHARACTER (len=data_name_length), DIMENSION(v3fit_work_size) ::          &
     &   c_work = ''

!  Model profile specification
!  Electron Denisty
!>  Model electron density profile, parameterized profile type.
!>  @see pprofile_T
      CHARACTER (len=p_type_len)           :: pp_ne_ptype = 'none'
!>  Array of b_coefficients for the electron density profile
!>  @see pprofile_T
      REAL (rprec), DIMENSION(ilb_b:iub_b) :: pp_ne_b = 0.0
!>  Array of as_coefficients electron density splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(iub_asf)     :: pp_ne_as = 0.0
!>  Array of af_coefficients electron density splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(iub_asf)     :: pp_ne_af = 0.0

!  Soft x-ray emissvity
!>  Model soft x-ray emissvity profile, parameterized profile type.
!>  @see pprofile_T
      CHARACTER (len=p_type_len)           :: pp_sxrem_ptype = 'none'
!>  Array of b_coefficients for the soft x-ray emissvity profile profile
!>  @see pprofile_T
      REAL (rprec), DIMENSION(ilb_b:iub_b) :: pp_sxrem_b = 0.0
!>  Array of as_coefficients soft x-ray emissvity profile splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(iub_asf)     :: pp_sxrem_as = 0.0
!>  Array of af_coefficients soft x-ray emissvity profile splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(iub_asf)     :: pp_sxrem_af = 0.0

!  Electron temperature
!>  Model electron temperature profile, parameterized profile type.
!>  @see pprofile_T
      CHARACTER (len=p_type_len)           :: pp_te_ptype = 'none'
!>  Array of b_coefficients for the electron temperature profile
!>  @see pprofile_T
      REAL (rprec), DIMENSION(ilb_b:iub_b) :: pp_te_b = 0.0
!>  Array of as_coefficients electron temperature splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(iub_asf)     :: pp_te_as = 0.0
!>  Array of af_coefficients electron temperature splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(iub_asf)     :: pp_te_af = 0.0

!  Model Specification Variables
!>  Specify how electron density is computed by the model.
!>  @see model
!>  @todo FIXME: New coding contains a model_ne_type that old code did not have.
!>  Make the default value pp_ne for compatibility. Eventially this shoud be set
!>  back to none and users should set model_ne_type in the namelist input.
      CHARACTER (len=data_name_length) :: model_ne_type = 'pp_ne'
!>  Specify how soft x-ray emissivity is computed by the model.
!>  @see model
      CHARACTER (len=data_name_length) :: model_sxrem_type = 'none'
!>  Specify how electron temperature is computed by the model.
!>  @see model
      CHARACTER (len=data_name_length) :: model_te_type = 'none'
!>  Normalization factor for the electron density.
      REAL (rprec)                     :: ne_pp_unit = 0.0
!>  Minimum electron density.
      REAL (rprec)                     :: ne_min = 0.0
!>  Specifies the fraction of the pressure constributed by the electrons.
      REAL (rprec)                     :: e_pressure_fraction = 0.0

!  Reconstruction constraints
!>  Not implemented.
      INTEGER                            :: n_rc = 0
!>  Not implemented.
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_parameters) :: rc_type = ''
!>  Not implemented.
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters) :: rc_index = 0
!>  Not implemented.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_parameters) :: rc_value = 0.0

!  Reconstruction parameters
!>  Number of reconstruction parameters.
      INTEGER                                :: n_rp = 0
!>  Names of reconstruction parameters.
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_parameters)     :: rp_type = ''
!>  Ith index of reconstruction parameters.
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters)     :: rp_index = 0
!>  Jth index of reconstruction parameters.
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters)     :: rp_index2 = 0
!>  Maximum reconstruction increment.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_parameters)     :: rp_vrnc = 0.0
!>  Parameter range types. The first index is the lower range and second index
!>  the upper range.
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_parameters,2) :: rp_range_type = 'infinity'
!>  Parameter range value. The first index is the lower range and second index
!>  the upper range.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_parameters,2)   :: rp_range_value = 0.0
!>  Parameter range indices. The first index is the lower range and second index
!>  the upper range.
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters,2,2) :: rp_range_index = 0

!  Signal data
!>  Number of observed signals.
      INTEGER                         :: n_sdata_o = 0
!>  Not implemented.
      INTEGER                         :: iw_sdo_verbose = -1
!>  Observed signals.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_signals) :: sdo_data_a = 0.0
!>  Observed sigmas.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_signals) :: sdo_sigma_a = 0.0
!>  Observed weights.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_signals) :: sdo_weight_a = 1.0

!  Signal data sigma specification
!>  Sigma specification lower index.
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   sdo_s_spec_imin = 0
!>  Sigma specification upper index.
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   sdo_s_spec_imax = 0
!>  Sigma specification floor.
      REAL (rprec), DIMENSION(v3fit_max_spec_size) ::                          &
     &   sdo_s_spec_floor = 0.0
!>  Sigma specification fraction.
      REAL (rprec), DIMENSION(v3fit_max_spec_size) ::                          &
     &   sdo_s_spec_fraction = 0.0

!  Signal data weight specification
!>  Weight specification lower index.
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   sdo_w_spec_imin = 0
!>  Weight specification upper index.
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   sdo_w_spec_imax = 0
!>  Weight specification weight.
      REAL (rprec), DIMENSION(v3fit_max_spec_size) ::                          &
     &   sdo_w_spec_weight = 0.0

!  External current signal specification
!>  External Z current radial position.
!>  @todo FIXME: The label of s0 and u0 are vmec specific. These should be
!>  changed to a radial coordinate and poloidal angle.
      REAL (rprec) :: extcurz_s0 = -1.0
!>  External Z current poloidal angle.
      REAL (rprec) :: extcurz_u0 = 0.0

!  Geometric information
!>  Geometric limiter major radius.
      REAL (rprec) :: r_major_radius
!>  Geometric limiter minor radius.
      REAL (rprec) :: a_minor_radius

!  Limiter iso function specification.
!>  Number of specified limiter iso functions.
      INTEGER                                       :: n_lif = 0
!>  Number of phi planes for an iso functions.
      INTEGER, DIMENSION(v3fit_max_geometrics)      :: n_phi_lif = 0
!>  Array of r-z polynomial coefficients.
      REAL (rprec), DIMENSION(v3fit_max_geometrics,0:4,0:4) ::                 &
     &   lif_arz = 0.0
!>  Array of r offset values.
      REAL (rprec), DIMENSION(v3fit_max_geometrics) ::                         &
     &   lif_rc = 0.0
!>  Array of z offset values.
      REAL (rprec), DIMENSION(v3fit_max_geometrics) ::                         &
     &   lif_zc = 0.0
!>  Array of sigma values.
      REAL (rprec), DIMENSION(v3fit_max_geometrics) ::                         &
     &   lif_sigma = 0.001_dp
!>  Array of phi values. Number of phi values need to be match @ref n_phi_lif
      REAL (rprec), DIMENSION(v3fit_max_geometrics,v3fit_max_lif_size)         &
     &   :: lif_phi_degree = 0.0
!>  True specifies that the edge hits the limiter. False specifies that the edge
!>  is inside the limiter.
      LOGICAL, DIMENSION(v3fit_max_geometrics)      ::                         &
     &   lif_on_edge = .false.

!  Use external boundary as signal
!>  @deprecated The implentation is VMEC specific. This will eventually be
!>  removed and replaced by priors when they are implemented.
      LOGICAL :: use_boundary

!  Combinations of other signals specification.
!>  Number of combination signals.
      INTEGER :: n_coosig = 0
!>  Number of signals in a combination signal.
      INTEGER, DIMENSION(v3fit_max_cos_size) :: n_sig_coosig = 0
!>  Indices of the number signals to combine. Number if indices needs to match
!>  @ref n_sig_coosig.
      INTEGER, DIMENSION(v3fit_max_combinations,v3fit_max_cos_size) ::         &
     &   coosig_indices = 0
!>  Coefficients of the signals to combine. Number if indices needs to match
!>  @ref n_sig_coosig.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_combinations,v3fit_max_cos_size) ::               &
     &   coosig_coeff = 0.0
!>  Combination type. @see signal
      CHARACTER(len=combination_type_length),                                  &
     &   DIMENSION(v3fit_max_combinations) :: coosig_type = 'sum'
!>  Name of the combination.
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_combinations) :: coosig_name = ''
!>  Units of the combination.
      CHARACTER (len=signal_units_length),                                     &
     &   DIMENSION(v3fit_max_combinations) :: coosig_units = ''

!  Reconstruction step control specification.
!>  Max number of reconstruction steps to perform.
      INTEGER                                :: nrstep = 0
!>  Stopping criterion on change in g^2.
      REAL (rprec)                           :: dg2_stop = 0.0
!>  Cutoff value for relative Singular Values.
      REAL (rprec)                           :: cut_svd = 0.0
!>  Cutoff value for expected step efficiency.
      REAL (rprec)                           :: cut_eff = 0.0
!>  Cutoff value for expected marginal step efficiency.
      REAL (rprec)                           :: cut_marg_eff = 0.0
!>  Cutoff value for expected step size.
      REAL (rprec)                           :: cut_delta_a = 0.0
!>  Cutoff value for expected change in g^2.
      REAL (rprec)                           :: cut_dg2 = 0.0
!>  Maximum allowable normalized step size.
      REAL (rprec)                           :: astep_max = 0.0
!>  Reconstruction step type.
!>  Possible values are:
!>    * 'sl'  straight line
!>    * 'lm'  Levenberg-Marquardt
!>    * 'seg' segmented
!>  @see reconstruction
      CHARACTER (len=v3fit_step_name_length) :: step_type = 'sl'

!>  Zero out vmec xcdot array
!>  @todo FIXME: VMEC specific. Not implemented.
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
!>  @deprecated The following arrays will be removed in the future. Users should
!>  use the @ref namelist_sec::recon_sec "reconstruction step control
!>  specifications" instead.
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
     &   n_rp, rp_type, rp_index, rp_index2, rp_vrnc, rp_range_type,           &
     &   rp_range_value, rp_range_index,                                       &
!  Signal Data specification
     &   n_sdata_o, iw_sdo_verbose, sdo_data_a, sdo_sigma_a,                   &
!  External current signal specification
!>  @todo FIXME: The label of s0 and u0 are vmec specific. These should be
!>  changed to a radial coordinate and poloidal angle.
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
!>  @deprecated The implementation is VMEC specific. This will eventually be
!>  removed and will replaced by priors when they are implemented.
     &   use_boundary,                                                         &
!  Combinations of other signals specification.
     &   n_coosig, n_sig_coosig, coosig_indices, coosig_coeff,                 &
     &   coosig_type, coosig_name, coosig_units,                               &
!  Reconstruction step control specification.
     &   nrstep, dg2_stop, cut_svd, cut_eff, cut_marg_eff, cut_delta_a,        &
     &   cut_dg2, astep_max, step_type,                                        &
!>  @todo FIXME: The implementation is VMEC specific.
     &   l_zero_xcdot

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Reads the namelist input file.
!>
!>  Reads the namelist input file. Replaces the sigma and weight values with the
!>  the computed specification. If using the old, major and minor radius, add a
!>  new limiter to the limiter specification. Determines the number of
!>  combination signals if the user has not specified @ref n_sig_coosig.
!>
!>  @param[in] namelist_file The file name of the namelist input file.
!-------------------------------------------------------------------------------
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
!>  @todo FIXME: Remove this when n_sig_coosig is mandatory.
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
