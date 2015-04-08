!-------------------------------------------------------------------------------
!  The @header, @table_section, @table_subsection, @item and @end_table commands
!  are custom defined commands in Doxygen.in. They are defined under ALIASES.
!  For the page created here, the 80 column limit is exceeded. Arguments of
!  aliases are separated by ','. If you intended ',' to be a string you must use
!  an escaped comma '\,'.
!
!>  @page namelist_sec Namelist v3fit_main_nli definition
!>
!>  @tableofcontents
!>  @section namelist_intro_sec Introduction
!>  This page documents the contents of a namelist input file. V3FIT namelist
!>  variables are defined in the @fixed_width{v3fit_main_nli} common block.
!>
!>  @section namelist_var_sec Namelist Variables
!>  @header{Input variable, Description, Code Reference}
!>
!>  @table_section{filename_sec, Filename variables}
!>     @item{main_nli_filename,           File name for main namelist input.,                                                                                         v3fit_input::main_nli_filename}
!>     @table_subsection{eq_filename_sec, Equilibrium Filename variables}
!>        @item{vmec_nli_filename,        File name for VMEC namelist input.,                                                                                         v3fit_input::vmec_nli_filename}
!>        @item{vmec_wout_input,          File name for the VMEC wout file to start from.,                                                                            v3fit_input::vmec_wout_input}
!>        @item{vacuum_nli_filename,      File name for a vacuum namelist input. The file specification is described in the @ref vacuum_namelist_sec.,                v3fit_input::vacuum_nli_filename}
!>     @table_subsection{signal_filename_sec, Singnal Filename variables}
!>        @item{mdsig_list_filename,      File name for list of MDSIG files.,                                                                                         v3fit_input::mdsig_list_filename}
!>        @item{sxrch_dot_filename,       File holding soft x-ray chord information. The file specification is described in the @ref sxrem_dot_sec.,                  v3fit_input::sxrch_dot_filename}
!>        @item{thscte_dot_filename,      File holding Thomson scattering information. The file specification is described in the @ref thomson_dot_sec.,              v3fit_input::thscte_dot_filename}
!>        @item{ipch_dot_filename,        File holding interferometry-polarimetry chord information. The file specification is described in the @ref intpol_dot_sec., v3fit_input::ipch_dot_filename}
!>        @item{mse_dot_filename,         File holding motional stark effect information. The file specification is described in the @ref mse_ratio_dot_sec.,         v3fit_input::mse_dot_filename}
!>        @item{sxrem_ratio_dot_filename, File holding soft x-ray emissivity ratio information. The file specification is described in the @ref sxrem_ratio_dot_sec., v3fit_input::sxrem_ratio_dot_filename}
!>  @end_table
!>
!>  @table_section{size_sec, Array allocation sizes\, numbers to allocate}
!>     @item{na_s_desc, s_desc array\, signal_desc, v3fit_input::na_s_desc}
!>  @end_table
!>
!>  @table_section{work_sec, Task specification\, work variables}
!>     @item{my_task, Character: to specify task in MAIN select case. Possible
!>                    tasks are:
!>                    -# @fixed_width{'equilibrium'} Solve the equilibrium.
!>                    -# @fixed_width{'v3post'}      Solve the equilibrium and compute modeled signals.
!>                    -# @fixed_width{'reconstruct'} Reconstruct the equilibirum.
!>                    -# @fixed_width{'units_tests'} Run the internal unit tests., v3fit_input::my_task}
!>  @end_table
!>
!>  @table_section{model_profile_sec, Model profile specification}
!>     @item{pp_ne_ptype,      Model electron density profile\, parameterized profile type.,      v3fit_input::pp_ne_ptype}
!>     @item{pp_ne_b,          Array of b_coefficients for the electron density profile.,         v3fit_input::pp_ne_b}
!>     @item{pp_ne_as,         Array of as_coefficients electron density splines.,                v3fit_input::pp_ne_as}
!>     @item{pp_ne_af,         Array of af_coefficients electron density splines.,                v3fit_input::pp_ne_af}
!>     @item{pp_sxrem_ptype,   Model soft x-ray emissivity profile\, parameterized profile type., v3fit_input::pp_sxrem_ptype}
!>     @item{pp_sxrem_b,       Array of b_coefficients for the soft x-ray emissivity profile.
!>                             @b DEPRICATED @b,                                                  v3fit_input::pp_sxrem_b}
!>     @item{pp_sxrem_as,      Array of as_coefficients soft x-ray emissivity splines.
!>                             @b DEPRICATED @b,                                                  v3fit_input::pp_sxrem_as}
!>     @item{pp_sxrem_af,      Array of af_coefficients soft x-ray emissivity splines.
!>                             @b DEPRICATED @b,                                                  v3fit_input::pp_sxrem_af}
!>     @item{num_sxrem_p,      Number of soft x-ray emissivity profiles.,                         v3fit_input::num_sxrem_p}
!>     @item{pp_sxrem_ptype_a, Model soft x-ray emissivity profile\, parameterized profile type., v3fit_input::pp_sxrem_ptype_a}
!>     @item{pp_sxrem_b_a,     Array of b_coefficients for the soft x-ray emissivity profile.,    v3fit_input::pp_sxrem_b_a}
!>     @item{pp_sxrem_as_a,    Array of as_coefficients soft x-ray emissivity splines.,           v3fit_input::pp_sxrem_as_a}
!>     @item{pp_sxrem_af_a,    Array of af_coefficients soft x-ray emissivity splines.,           v3fit_input::pp_sxrem_af_a}
!>     @item{pp_te_ptype,      Model electron temperature profile\, parameterized profile type.,  v3fit_input::pp_te_ptype}
!>     @item{pp_te_b,          Array of b_coefficients for the e-temperature profile.,            v3fit_input::pp_te_b}
!>     @item{pp_te_as,         Array of as_coefficients e-temperature splines.,                   v3fit_input::pp_te_as}
!>     @item{pp_te_af,         Array of af_coefficients e-temperature splines.,                   v3fit_input::pp_te_af}
!>     @item{sxrem_te,         Array of temperature points for the sxrem ratio function.,         v3fit_input::sxrem_te_a}
!>     @item{sxrem_ratio,      Array of ratio points for the sxrem ratio function.,               v3fit_input::sxrem_ratio_a}
!>  @end_table
!>
!>  @table_section{model_spec_sec, Model Specification Variables}
!>     @item{model_ne_type, Character variable to specify how the electron density will be calculated by the model.
!>        -# @fixed_width{'none'}    Not part of the model.
!>        -# @fixed_width{'pp_ne'}   Calculated from the parameterized profile pp_ne.
!>        -# @fixed_width{'pp_te_p'} Calculated from the pp_te and the VMEC pressure,
!>        v3fit_input::model_ne_type}
!>     @item{model_sxrem_type, Character variable to specify how the Soft x-ray emissivity will be calculated by the model.
!>        -# @fixed_width{'none'}     Not part of the model.
!>        -# @fixed_width{'pp_sxrem'} Calculated from the parameterized profile pp_sxrem.
!>
!>        @b DEPRICATED @b,
!>        v3fit_input::model_sxrem_type}
!>     @item{model_sxrem_type_a, Character array to specify how the Soft x-ray emissivity will be calculated by the model.
!>        -# @fixed_width{'none'}           Not part of the model.
!>        -# @fixed_width{'pp_sxrem'}       Calculated from the parameterized profile pp_sxrem.
!>        -# @fixed_width{'pp_sxrem_te_ne'} Calculated from pp_te and pp_ne.,
!>        v3fit_input::model_sxrem_type_a}
!>     @item{model_te_type, Character variable to specify how the e-temperature will be calculated by the model.
!>        -# @fixed_width{'none'}                    Not part of the model.
!>        -# @fixed_width{'pp_te'}                   Calculated from the parameterized profile pp_te.
!>        -# @fixed_width{'pp_ne_vmec_p'\,'pp_ne_p'} Calculated from the pp_ne and the VMEC pressure.,
!>        v3fit_input::model_te_type}
!>     @item{ne_pp_unit,          Scaling/conversion factor for the density. ne * ne_pp_unit is assumed to have units of m^-3.,                                 v3fit_input::ne_pp_unit}
!>     @item{ne_min,              Minimum electron density\, m^-3.,                                                                                             v3fit_input::ne_min}
!>     @item{te_min,              Minimum electron temperature\, eV.,                                                                                           v3fit_input::te_min}
!>     @item{e_pressure_fraction, Electron pressure fraction of the total pressure. Used when model_te_type = 'pp_ne_vmec_p' or model_ne_type = 'pp_te_vmec_p', v3fit_input::e_pressure_fraction}
!>     @item{emission_file,       Path to file defining the soft x-ray emission. @see @ref emission_file_sec ,                                                  v3fit_input::emission_file}
!>  @end_table
!>
!>  @table_section{recon_con_sec, Reconstruction constraints (DEPRICATED)}
!>     @item{n_rc,     Number of reconstruction constraints,       v3fit_input::n_rc}
!>     @item{rc_type,  Array of reconstruction constraint types,   v3fit_input::rc_type}
!>     @item{rc_index, Array of reconstruction constraint indices, v3fit_input::rc_index}
!>     @item{rc_value, Array of reconstruction constraint values,  v3fit_input::rc_value}
!>  @end_table
!>
!>  @table_section{derived_param_sec, Derived parameter specification}
!>     @item{n_dp,     Number of derived parameters,                            v3fit_input::n_dp}
!>     @item{dp_type,  Array of derived parameter types. Any parameter can be a
!>                     * @ref vmec_equilibrium_derive_param_sec
!>                     * @ref model_derived_param_sec,                          v3fit_input::dp_type}
!>     @item{dp_index, Array of derived parameter indices,                      v3fit_input::dp_index}
!>  @end_table
!>
!>  @table_section{recon_param_sec, Reconstruction parameter specification}
!>     @item{n_rp,           Number of reconstruction parameters,                                v3fit_input::n_rp}
!>     @item{rp_type,        Array of reconstruction parameter types. Valid reconstruction
!>                           parameters are defined by the @ref model or the @ref equilibirum.
!>                           * @ref model_recon_param_sec
!>                           * @ref vmec_equilibrium_recon_param_sec
!>                           * @ref vacuum_equilibrium_recon_param_sec,                          v3fit_input::rp_type}
!>     @item{rp_index,       Array of reconstruction parameter indices,                          v3fit_input::rp_index}
!>     @item{rp_index2,      Array of reconstruction parameter 2nd indices,                      v3fit_input::rp_index2}
!>     @item{rp_vrnc,        Array of reconstruction parameter variances,                        v3fit_input::rp_vrnc}
!>     @item{rp_range_type,  2D Array of reconstruction parameter range types
!>                           -# @fixed_width{'infinity'} Unbounded
!>                           -# @fixed_width{'value'}    Bounded by the value of rp_range_value
!>                           -# @fixed_width{\<other\>}  Bounded by an equilibrium parameter.
!>
!>                           Any parameter maybe used as a constraint.
!>                           * @ref model_sec
!>                           * @ref vmec_equilibrium_sec
!>                           * @ref vacuum_equilibrium_sec
!>
!>                           The second index contains the low and high range types
!>                           respectively.,                                                      v3fit_input::rp_range_type}
!>     @item{rp_range_value, 2D Array of values used if rp_range_type = @fixed_width{'value'}.
!>                           The second index contains the low and high values repectively.,     v3fit_input::rp_range_value}
!>     @item{rp_range_index, 3D Array of values of indicies used if\,
!>                           rp_range_type = @fixed_width{\<other\>}.
!>                           The second index specifies the low and high ranges repectively.
!>                           The third index contains the first and second index of the
!>                           inequality constraint parameter.,                                   v3fit_input::rp_range_index}
!>  @end_table
!>
!>  @table_section{signal_spec_sec, Signal Data specification}
!>     @item{n_sdata_o,               Number of signal_data observations @b DEPRICATED @b,                                                 v3fit_input::n_sdata_o}
!>     @item{iw_sdo_verbose,          Integer to control write out of signal_data observations. (-1\, no write\, otherwise\, verbosity),   v3fit_input::iw_sdo_verbose}
!>     @item{sdo_data_a,              Array of data for observed signals,                                                                  v3fit_input::sdo_data_a}
!>     @item{sdo_sigma_a,             Array of sigmas for observed signals,                                                                v3fit_input::sdo_sigma_a}
!>     @item{sdo_weight_a,            Array of weight for observed signals,                                                                v3fit_input::sdo_weight_a}
!>     @item{mag_a,                   Array of flags for magnetc signals. This flags controls if the induced signal is included in the
!>                                    total signal.,                                                                                       v3fit_input::mag_a}
!>     @item{mag_force_coil,          Flag to force the computation of induced signals. This has no affect on the total signal if
!>                                    @ref mag_a is false.,                                                                                v3fit_input::mag_force_coil}
!>     @item{sdo_spec_can_overwrite, Logical controlling if the sdo_*_spec_ overwrites values in the sdo_sigma_a array. Defaults to true., v3fit_input::sdo_spec_can_overwrite}
!>     @table_subsection{sigma_spec_sec, Signal Data sigma specification}
!>        @item{sdo_s_spec_imin,      Array of minimum index values for specifying sigmas,                                                 v3fit_input::sdo_s_spec_imin}
!>        @item{sdo_s_spec_imax,      Array of maximum index values for specifying sigmas,                                                 v3fit_input::sdo_s_spec_imax}
!>        @item{sdo_s_spec_floor,     Array of floor values for specifying sigmas,                                                         v3fit_input::sdo_s_spec_floor}
!>        @item{sdo_s_spec_fraction,  Array of fractional values for specifying sigmas
!>            - floor \> 0\, sigma = max(floor\,fraction * sdo)
!>            - floor \< 0\, sigma = sqrt(floor ** 2 + (fraction * sdo) ** 2),                                                             v3fit_input::sdo_s_spec_fraction}
!>     @table_subsection{weight_spec_sec, Signal Data weight specification}
!>        @item{sdo_w_spec_imin,      Array of minimum index values for specifying weights,                                                v3fit_input::sdo_w_spec_imin}
!>        @item{sdo_w_spec_imax,      Array of maximum index values for specifying weights,                                                v3fit_input::sdo_w_spec_imax}
!>        @item{sdo_w_spec_weight,    Array of weight values for specifying weights,                                                       v3fit_input::sdo_w_spec_weight}
!>     @table_subsection{magnetic_spec_sec, Magnetic signal specification}
!>        @item{mag_spec_imin,        Array of minimum index values for specifying plasma only signals,                                    v3fit_input::mag_spec_imin}
!>        @item{mag_spec_imax,        Array of maximum index values for specifying plasma only signals,                                    v3fit_input::mag_spec_imax}
!>        @item{mag_spec_use_induced, Array of logicals for specifying plasma only signals,                                                v3fit_input::mag_spec_use_induced}
!>  @end_table
!>
!>  @table_section{magnetic_sec, Magnetic signal parameters}
!>     @item{pol_rad_ratio, Ratio of the radial grid points to the poloidal grid points, v3fit_input::pol_rad_ratio}
!>  @end_table
!>
!>  @todo FIXME: S and U are VMEC specific. Make these a radial position and poloidal angle.
!>  @table_section{extcurz_sec, These two are used by the extcurz diagnostic signal.}
!>     @item{extcurz_s0, S position to integrate about., v3fit_input::extcurz_s0}
!>     @item{extcurz_u0, U position to integrate about., v3fit_input::extcurz_u0}
!>  @end_table
!>
!>  @table_section{geometric_sec, Geometric information}
!>     @item{r_major_radius,    For circular torus geometric,                       v3fit_input::r_major_radius}
!>     @item{a_minor_radius,    For circular torus geometric,                       v3fit_input::a_minor_radius}
!>     @table_subsection{lif_sec, Limiter_Iso Function.
!>                                First dimension of arrays is dimensioned na_lif}
!>        @item{n_lif,          number of limiter_iso functions,                    v3fit_input::n_lif}
!>        @item{n_phi_lif,      array of number of phi values (:),                  v3fit_input::n_phi_lif}
!>        @item{lif_arz,        array of r-z polynomial coefficients (:\,0:4\,0:4), v3fit_input::lif_arz}
!>        @item{lif_rc,         array of r offset values (:),                       v3fit_input::lif_rc}
!>        @item{lif_zc,         array of z offset values (:),                       v3fit_input::lif_zc}
!>        @item{lif_sigma,      array of sigma values (:),                          v3fit_input::lif_sigma}
!>        @item{lif_phi_degree, array of phi values (:\,na_phis_lif),               v3fit_input::lif_phi_degree}
!>        @item{lif_on_edge,    array of logical values (:)
!>           - true  - edge of plasma is on limiter
!>           - false - edge of plasma is within limiter,
!>           v3fit_input::lif_on_edge}
!>  @end_table
!>
!>  @table_section{prior_sec, Prior information}
!>     @item{n_prior,          Number of specified priors.,             v3fit_input::n_prior}
!>     @item{prior_name,       Name of the specified priors.,           v3fit_input::prior_name}
!>     @item{prior_param_name, Name of the prior parameter. Any model
!>                             or equilirbium parameter maybe used as
!>                             a prior.
!>                             * @ref model_sec
!>                             * @ref vmec_equilibrium_sec
!>                             * @ref vacuum_equilibrium_sec,           v3fit_input::prior_param_name}
!>     @item{prior_indices,    2D array of the ith and jth index of
!>                             the parameter. The second index
!>                             contains the first and second indicies
!>                             repectively.,                            v3fit_input::prior_indices}
!>     @item{prior_units,      The units of the prior signal.,          v3fit_input::prior_units}
!>  @end_table
!>
!>  @table_section{coosig_sec, COOSIG - Combination Of Other SIGnals NLI Variables}
!>     @item{n_coosig,       Number of new signals of signal_type coosig,                                                                                                           v3fit_input::n_coosig}
!>     @item{n_sig_coosig,   Array of number of signals in a combinaton.,                                                                                                           v3fit_input::n_sig_coosig}
!>     @item{coosig_indices, 2d Integer Array of indices (to list of signals) for combinations. First index is for which new signal\, second index is for terms in the combination, v3fit_input::coosig_indices}
!>     @item{coosig_coeff,   2d real array of coefficients for the combination (a's in the type),                                                                                   v3fit_input::coosig_coeff}
!>     @item{coosig_type,    Character array of combination types
!>        -# @fixed_width{'sum'}  - linear combination
!>        -# @fixed_width{'max'}  - maximum of weighted signals
!>        -# @fixed_width{'min'}  - minimum of weighted signals,
!>        v3fit_input::coosig_type}
!>     @item{coosig_name,    Character array of coosig names,                                                                                                                       v3fit_input::coosig_name}
!>     @item{coosig_units,   Character array of coosig units,                                                                                                                       v3fit_input::coosig_units}
!>  @end_table
!>
!>  @table_section{plasma_off_sec, Plasma Offset Section}
!>     @item{phi_offset, Initial phi offset value to rotate the plasma relative
!>                       to the diagnostics. This value is in radians.,           v3fit_input::phi_offset}
!>     @item{z_offset,   Initial z offset value to vertical shift plasma relative
!>                       to the center. This value is in meters.,                 v3fit_input::z_offset}
!>  @end_table
!>
!>  @table_section{recon_sec, Reconstruction step control NLI variables}
!>     @item{nrstep,       max number of reconstruction steps to perform,      v3fit_input::nrstep}
!>     @item{dg2_stop,     Stopping criterion on change in g^2,                v3fit_input::dg2_stop}
!>     @item{cut_svd,      cutoff value for relative Singular Values,          v3fit_input::cut_svd}
!>     @item{cut_eff,      cutoff value for expected step efficiency,          v3fit_input::cut_eff}
!>     @item{cut_marg_eff, cutoff value for expected marginal step efficiency, v3fit_input::coosig_units}
!>     @item{cut_delta_a,  cutoff value for expected step size,                v3fit_input::cut_delta_a}
!>     @item{cut_dg2,      cutoff value for expected change in g^2,            v3fit_input::cut_dg2}
!>     @item{astep_max,    maximum allowable normalized step size,             v3fit_input::astep_max}
!>     @item{step_type,    character specification of reconstruction step type
!>        -# @fixed_width{'sl'}  - straight line
!>        -# @fixed_width{'seg'} - segmented
!>        -# @fixed_width{'lm'}  - Levenberg-Marquardt,
!>        v3fit_input::step_type}
!>     @item{cut_inv_svd,  Cut off value for the singular values used for the
!>                         pseudo matrix inverse,                              v3fit_input::cut_inv_svd}
!>  @end_table
!>
!>  @table_section{xcdot_sec, VMEC zero out xcdot control}
!>     @item{l_zero_xcdot, reset the xcdot array to zero when converging vmec., v3fit_input::l_zero_xcdot}
!>
!>  @section namelist_exam_sec Example Files
!>  @subsection namelist_exam_v3post_sec V3POST Example
!>  @code
!>  ! Example V3POST input file.
!>  &v3fit_main_nli
!>
!>  ! Equilibrium to use.
!>  vmec_nli_filename = 'input.example.vmec'
!>
!>  ! Diagnostic inout files.
!>  mdsig_list_file = 'diagnostic.example_mdsig.LIST'
!>
!>  my_task = 'v3post'
!>  /
!>  @endcode
!>  @subsection namelist_exam_reconstruction_sec Reconstruction Example
!>  @code
!>  ! Example Reconstruction input file.
!>  &v3fit_main_nli
!>
!>  ! Equilibrium to use.
!>  vmec_nli_filename = 'input.example.vmec'
!>
!>  ! Diagnostic inout files.
!>  mdsig_list_file = 'diagnostic.example_mdsig.LIST'
!>  sxrch_dot_filename = 'sxrch.example'
!>  thscte_dot_filename = 'thscte.example'
!>
!>  my_task = 'reconstruct'
!>
!>  ! Model Specification
!>
!>  ! Electron Density
!>  pp_ne_ptype = 'two_power'
!>
!>  ! b(0) + b(1)*[1 - s^b(2)]^b(3)
!>  pp_ne_b(0:3) = 0.0, 1.0, 5.0, 10.0
!>  ne_pp_unit = 1.0E18
!>  e_pressure_fraction=0.75
!>  
!>  ! Soft X-ray emissivity
!>  num_sxrem_p = 1
!>  model_sxrem_type_a(1) = 'pp_sxrem'
!>  pp_sxrem_ptype_a(1) = 'two_power'
!>
!>  !>  ! b(0) + b(1)*[1 - s^b(2)]^b(3)
!>  pp_sxrem_b_a(1,0:3) = 0.0, 1.0, 3.0, 10.0
!>
!>  ! Electron Temperature
!>  model_te_type = 'pp_ne_vmec_p'
!>
!>  ! Reconstruction control parameters.
!>  nrstep = 20
!>  dg2_stop = 0.1
!>  cut_svd = 0.002
!>  cut_eff = 0.1
!>  cut_marg_eff = 0.0
!>  cut_delta_a = 0.0
!>  cut_dg2 = 0.0
!>  astep_max = 100.0
!>  step_type = 'lm'
!>
!>  ! Reconstruction parameters
!>  n_rp = 5
!>
!>  rp_type(1) = 'phiedge'
!>  rp_vrnc(1) = 1.E-04
!>
!>  rp_type(2) = 'ac'
!>  rp_index(2) = 1
!>  rp_vrnc(2) = 0.01
!>
!>  rp_type(3) = 'pres_scale'
!>  rp_vrnc(3) = 50.
!>
!>  rp_type(4)='pp_ne_b'
!>  rp_index(4)=1
!>  rp_vrnc(4)=0.2
!>
!>  rp_type(5)='pp_ne_b'
!>  rp_index(5)=2
!>  rp_vrnc(5)=0.04
!>
!>  rp_type(5)='pp_sxrem_b'
!>  rp_index(5)=1
!>  rp_vrnc(5)=0.1
!>
!>  rp_type(6)='pp_sxrem_b'
!>  rp_index(6)=2
!>  rp_vrnc(6)=0.1
!>
!>  rp_type(7)='pp_sxrem_b'
!>  rp_index(7)=3
!>  rp_vrnc(7)=0.1
!>  
!>  ! Limiter_iso Functions
!>  n_lif = 2
!>  
!>  n_phi_lif(1) = 1
!>  lif_rc(1) = 0.7493
!>  lif_on_edge(1) = T
!>  lif_sigma(1) = 0.001
!>  lif_arz(1,0,0) = -0.06002,
!>  lif_arz(1,2,0) = 1.0,
!>  lif_arz(1,0,2) = 1.0,
!>  lif_phi_degree(1,1) = 12.0,
!>  
!>  n_phi_lif(2) = 1
!>  lif_rc(2) = 0.7493
!>  lif_on_edge(2) = T
!>  lif_arz(2,0,0) = -0.0618,
!>  lif_arz(2,2,0) = 1.0,
!>  lif_arz(2,0,2) = 1.0,
!>  lif_phi_degree(2,1) = 60.0,
!>  
!>  ! Signal Weight specification - set weights for coosig parts to zero
!>  sdo_w_spec_imin = 12 16 173 174
!>  sdo_w_spec_imax = 12 16 173 174
!>  sdo_w_spec_weight = 0. 0. 0. 0.
!>  
!>  ! Combination Of Other SIGnals (COOSIG)
!>  n_coosig = 2
!>  
!>  coosig_type(1) = 'sum'     ! Linear combination, to simulate a compensated diamagnetic loop
!>  coosig_indices(1,1)=12     !  Here I use the two bad signals, that already have their
!>  coosig_indices(1,2)=16     !  weights set to zero
!>  coosig_coeff(1,1)=.389976  !  Set to make zero signal, using *e model values
!>  coosig_coeff(1,2)=-.850293
!>  coosig_name(1) = 'simCompDiaLoop'
!>  coosig_units(1) = '?'
!>  
!>  coosig_type(2) = 'max'     ! Maximum, for use with edge limits
!>  coosig_indices(2,1)=173     ! Signal numbers for the two edge limits
!>  coosig_indices(2,2)=174
!>  coosig_coeff(2,1)=1        ! No scaling of the relative signals before the max
!>  coosig_coeff(2,2)=1
!>  coosig_name(2) = 'OneLimiter'
!>  coosig_units(2) = 'm'
!>  
!>  ! Observed Data - magnetic diagnostics
!>  n_sdata_o = 25
!>  sdo_data_a(1)= 0.099877    sdo_sigma_a(1)=0.00177   !  First Magnetic Diagnostics
!>  sdo_data_a(2)= 0.170352    sdo_sigma_a(2)=0.00177
!>  snip - boring observations
!>  sdo_data_a(25)= 0.031115   sdo_sigma_a(25)=0.00007   ! Last Magnetic Diagnostic
!>  
!>  ! Observed data for the soft x-ray chords
!>  ! Simulated data, see Excel "Simulated Observations 2011-10-26"
!>  n_sdata_o = 125
!>  sdo_s_spec_imin(1) = 26
!>  sdo_s_spec_imax(1) = 125
!>  sdo_s_spec_floor(1) = 1.
!>  sdo_s_spec_fraction(1) = 0.05
!>  
!>  sdo_data_a(26) = 6.09338E+01      ! sxr  1
!>  snip - boring simulated obserations
!>  sdo_data_a(125) =7.88700E+00       ! sxr100
!>  
!>  ! Observed data for the Thomson Scattering points
!>  ! Simulated data, see Excel "Simulated Observations 2011-10-26"
!>  n_sdata_o = 174
!>  sdo_s_spec_imin(2) = 126
!>  sdo_s_spec_imax(2) = 172
!>  sdo_s_spec_floor(2) = 1.
!>  sdo_s_spec_fraction(2) = 0.05
!>  sdo_w_spec_imin(5) = 126
!>  sdo_w_spec_imax(5) = 172
!>  sdo_w_spec_weight(5) = 0.
!>  
!>  sdo_data_a(126) = 6.82686E+00    ! thscte
!>  snip - boring simulated observations
!>  sdo_data_a(170) = 1.10650E+00    ! thscte
!>  sdo_data_a(171) = 2.26904E+00    ! thscte
!>  sdo_data_a(172) = 3.86471E+00    ! thscte
!>  sdo_data_a(173)= 0.000000   sdo_sigma_a(173)=0.00100  !  Limiter
!>  sdo_data_a(174)= 0.000000   sdo_sigma_a(174)=0.00100
!>  
!>  !  Add on observations for the COOSIGs
!>  n_sdata_o = 176
!>  sdo_data_a(175) = 0.        !  'simCompDiaLoop', coosig_coeff set to make signal zero
!>  sdo_sigma_a(175) = 1.5E-03
!>  sdo_data_a(176) = 0.        !  'OneLimiter', signal should be zero
!>  sdo_sigma_a(176) = 0.001    !   1 mm
!>  /
!>  @endcode
!>
!>  @section namelist_prog_ref_sec Programmers Reference
!>  Reference material for the coding to implement this namelist is found in the
!>  @ref v3fit_input module.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file v3fit_input.f
!>  @brief Contains module @ref v3fit_input.
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

      IMPLICIT NONE

!*******************************************************************************
!  v3fit input module parameters
!*******************************************************************************
!>  Maximum number of diagnostic signals.
      INTEGER, PARAMETER :: v3fit_max_diagnostics = 1000
!>  Maximum number of geometric signals.
      INTEGER, PARAMETER :: v3fit_max_geometrics = 1000
!>  Maximum number of prior signals.
      INTEGER, PARAMETER :: v3fit_max_priors = 1000
!>  Maximum number of combination signals.
      INTEGER, PARAMETER :: v3fit_max_combinations = 1000
!>  Maximum number of total signals.
      INTEGER, PARAMETER :: v3fit_max_signals = v3fit_max_diagnostics          &
     &                                        + v3fit_max_geometrics           &
     &                                        + v3fit_max_priors               &
     &                                        + v3fit_max_combinations

!>  Maximum number of reconstruction parameters.
      INTEGER, PARAMETER :: v3fit_max_parameters = 100

!>  Maximum size of the easy specification arrays.
!>  @see sdo_s_spec_imin
!>  @see sdo_s_spec_imax
!>  @see sdo_s_spec_floor
!>  @see sdo_s_spec_fraction
!>  @see sdo_w_spec_imin
!>  @see sdo_w_spec_imax
!>  @see sdo_w_spec_weight
      INTEGER, PARAMETER :: v3fit_max_spec_size = 150

!>  Maximum number of phi positions for limiters.
      INTEGER, PARAMETER :: v3fit_max_lif_size = 100

!>  Maximum number of signals a combination signal can refer to.
      INTEGER, PARAMETER :: v3fit_max_cos_size = 100

!>  Reconstruction step name length.
      INTEGER, PARAMETER :: v3fit_step_name_length = 4

!>  Maximum number of sxrem profiles.
      INTEGER, PARAMETER :: v3fit_max_sxrem_profiles = 10

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) v3fit_namelist_class
!
!*******************************************************************************
!>  File name for main namelist input.
!>  @ref namelist_sec
      CHARACTER (len=data_filename_length) :: main_nli_filename = ''

!  Equilibrium input files
!>  File name for VMEC namelist input
      CHARACTER (len=data_filename_length) :: vmec_nli_filename = ''
!>  File name for the VMEC wout input.
      CHARACTER (len=data_filename_length) :: vmec_wout_input = ''
!>  File name for a vacuum namelist input
      CHARACTER (len=data_filename_length) :: vacuum_nli_filename = ''

!  Signal input files
!>  File name for list of MDSIG files
      CHARACTER (len=data_filename_length) :: mdsig_list_filename = ''
!>  File holding sxr chord information
      CHARACTER (len=data_filename_length) :: sxrch_dot_filename = ''
!>  File holding interferometry-polarimetry information.
!>  @ref intpol_dot_sec
      CHARACTER (len=data_filename_length) :: ipch_dot_filename = ''
!>  File holding Thomson Scattering information
      CHARACTER (len=data_filename_length) :: thscte_dot_filename = ''
!>  File holding motional stark effect information
      CHARACTER (len=data_filename_length) :: mse_dot_filename = ''
!>  File holding soft x-ray ratio information
      CHARACTER (len=data_filename_length) ::                                  &
     &   sxrem_ratio_dot_filename = ''

!  Array allocation sizes. These should be set to the maxmium size.
!>  Maximum number of signals to create.
!>  @note this needs to be greater than or equal to total number of signals
!>  used. The signal arrays are resized to free up memory if this value is
!>  larger than the number of signals constructed. If this value is equal to the
!>  number of signals, the resize step can be avoided.
      INTEGER :: na_s_desc = v3fit_max_signals

!  Task specification and work variables
!>  The v3fit task. Possible values are.
!>  -# @fixed_width{equilibrium} Solve the equilibrium.
!>  -# @fixed_width{v3post}      Solve the equilibrium and compute modeled signals.
!>  -# @fixed_width{reconstruct} Reconstruct the equilibirum.
!>  -# @fixed_width{units_tests} Run the internal unit tests.
      CHARACTER (len=data_name_length)         :: my_task = ''

!  Model profile specification
!  Electron Denisty
!>  Model electron density profile, parameterized profile type.
!>  -# @fixed_width{'none'}    Not part of the model.
!>  -# @fixed_width{'pp_ne'}   Calculated from the parameterized profile pp_ne.
!>  -# @fixed_width{'pp_te_p'} Calculated from the pp_te and the VMEC pressure,
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
!>  @deprecated
!>  Model soft x-ray emissvity profile, parameterized profile type.
!>  -# @fixed_width{'none'}     Not part of the model.
!>  -# @fixed_width{'pp_sxrem'} Calculated from the parameterized profile pp_sxrem.
!>  @see pprofile_T
      CHARACTER (len=p_type_len)           :: pp_sxrem_ptype = 'none'
!>  @deprecated
!>  Array of b_coefficients for the soft x-ray emissvity profile
!>  @see pprofile_T
      REAL (rprec), DIMENSION(ilb_b:iub_b) :: pp_sxrem_b = 0.0
!>  @deprecated
!>  Array of as_coefficients soft x-ray emissvity profile splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(iub_asf)     :: pp_sxrem_as = 0.0
!>  @deprecated
!>  Array of af_coefficients soft x-ray emissvity profile splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(iub_asf)     :: pp_sxrem_af = 0.0

!  Multi soft x-ray emissvity profile specification.
!>  Number of sxrem profiles.
      INTEGER                              :: num_sxrem_p = 1
!>  Model soft x-ray emissvity profile, parameterized profile type.
!>  -# @fixed_width{'none'}     Not part of the model.
!>  -# @fixed_width{'pp_sxrem'} Calculated from the parameterized profile pp_sxrem.
!>  @see pprofile_T
      CHARACTER (len=p_type_len), DIMENSION(v3fit_max_sxrem_profiles)          &
     &   :: pp_sxrem_ptype_a = 'none'
!>  Array of b_coefficients for the soft x-ray emissvity profile
!>  @see pprofile_T
      REAL (rprec), DIMENSION(v3fit_max_sxrem_profiles,ilb_b:iub_b)            &
     &   :: pp_sxrem_b_a = 0.0
!>  Array of as_coefficients soft x-ray emissvity profile splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(v3fit_max_sxrem_profiles,iub_asf)                &
     &   :: pp_sxrem_as_a = 0.0
!>  Array of af_coefficients soft x-ray emissvity profile splines
!>  @see pprofile_T
      REAL (rprec), DIMENSION(v3fit_max_sxrem_profiles,iub_asf)                &
     &   :: pp_sxrem_af_a = 0.0

!  Electron temperature
!>  Model electron temperature profile, parameterized profile type.
!>  -# @fixed_width{'none'}                    Not part of the model.
!>  -# @fixed_width{'pp_te'}                   Calculated from the parameterized profile pp_te.
!>  -# @fixed_width{'pp_ne_vmec_p'\,'pp_ne_p'} Calculated from the pp_ne and the VMEC pressure.
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

!  Soft X-ray emission ratio
!>  Array of temperature points for the sxrem ratio function.
      REAL (rprec), DIMENSION(iub_asf) :: sxrem_te_a = 0.0
!>  Array of ratio points for the sxrem ratio function.
      REAL (rprec), DIMENSION(iub_asf) :: sxrem_ratio_a = 0.0

!  Model Specification Variables
!>  Specify how electron density is computed by the model.
!>  @see model
!>  @todo FIXME: New coding contains a model_ne_type that old code did not have.
!>  Make the default value pp_ne for compatibility. Eventially this shoud be set
!>  back to none and users should set model_ne_type in the namelist input.
      CHARACTER (len=data_name_length) :: model_ne_type = 'pp_ne'
!>  @deprecated
!>  Specify how soft x-ray emissivity is computed by the model.
!>  @see model
      CHARACTER (len=data_name_length) :: model_sxrem_type = 'none'
!>  Specify how soft x-ray emissivity is computed by the model.
!>  @see model
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_sxrem_profiles) ::                                &
     &      model_sxrem_type_a = 'none'
!>  Specify how electron temperature is computed by the model.
!>  @see model
      CHARACTER (len=data_name_length) :: model_te_type = 'none'
!>  Normalization factor for the electron density.
      REAL (rprec)                     :: ne_pp_unit = 1.0E18
!>  Minimum electron density.
      REAL (rprec)                     :: ne_min = 0.0
!>  Minimum electron temperature.
      REAL (rprec)                     :: te_min = 0.0
!>  Specifies the fraction of the pressure constributed by the electrons.
      REAL (rprec)                     :: e_pressure_fraction = 0.5
!>  Specifies the path to the file contain the soft x-ray emission function.
      CHARACTER (len=data_name_length) :: emission_file = ''

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

!  Derived parameters
!>  Number of derived parameters.
      INTEGER                              :: n_dp = 0
!>  Names of derived parameters.
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_parameters)   :: dp_type = ''
!>  Indices of derived parameters.
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters,2) :: dp_index = 0

!  Reconstruction parameters
!>  Number of reconstruction parameters.
      INTEGER                              :: n_rp = 0
!>  Names of reconstruction parameters.
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_parameters)   :: rp_type = ''
!>  Ith index of reconstruction parameters.
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters)   :: rp_index = 0
!>  Jth index of reconstruction parameters.
      INTEGER,                                                                 &
     &   DIMENSION(v3fit_max_parameters)   :: rp_index2 = 0
!>  Maximum reconstruction increment.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_parameters)   :: rp_vrnc = 0.0
!>  Parameter range types. The first index is the lower range and second index
!>  the upper range.
!>  -# @fixed_width{'infinity'} Unbounded
!>  -# @fixed_width{'value'}    Bounded by the value of rp_range_value
!>  -# @fixed_width{\<other\>}  Bounded by an equilibrium parameter.
      CHARACTER (len=data_name_length),                                        &
     &   DIMENSION(v3fit_max_parameters,2) :: rp_range_type = 'infinity'
!>  Parameter range value. The first index is the lower range and second index
!>  the upper range.
      REAL (rprec),                                                            &
     &   DIMENSION(v3fit_max_parameters,2) :: rp_range_value = 0.0
!>  Parameter range indices. The first index is the lower range and second index
!>  the upper range.
      INTEGER, DIMENSION(v3fit_max_parameters,2,data_max_indicies) ::          &
     &   rp_range_index = 0

!  Signal data
!>  @deprecated
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
!>  Allow the sigma spec to overwrite the sigma array.
      LOGICAL                         :: sdo_spec_can_overwrite = .true.
!>  Magnetic signal flags. Controls if the induced signal is included in the
!>  total signal.
      LOGICAL, DIMENSION(v3fit_max_diagnostics) :: mag_a = .true.
!>  Forces the magnetic signal to compute the induced signal.
      LOGICAL                         :: mag_force_coil = .false.

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

!  Magnetic specification
!>  Magnetic specification lower index.
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   mag_spec_imin = 0
!>  Magnetic specification upper index.
      INTEGER, DIMENSION(v3fit_max_spec_size)      ::                          &
     &   mag_spec_imax = 0
!>  Magnetic specification toggle induced coil signal.
      LOGICAL, DIMENSION(v3fit_max_spec_size) ::                               &
     &   mag_spec_use_induced = .true.

!  Magnetic integration parameters
!>  Ratio of the number of poloidal grid points to radial grid points.
      REAL (rprec) :: pol_rad_ratio = 1.0

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

!  Prior signal specification.
!>  Number of specified priors.
      INTEGER                                       :: n_prior = 0
!>  Prior names.
      CHARACTER (len=data_name_length), DIMENSION(v3fit_max_priors)            &
     &   :: prior_name = ''
!>  Prior parameter names.
      CHARACTER (len=data_name_length), DIMENSION(v3fit_max_priors)            &
     &   :: prior_param_name = ''
!>  Prior parameter indicies.
      INTEGER, DIMENSION(v3fit_max_priors,2) :: prior_indices = 0
!>  Prior Units
      CHARACTER (len=data_short_name_length),                                  &
     &   DIMENSION(v3fit_max_priors) :: prior_units = ''

!  Combinations of other signals specification.
!>  Number of combination signals.
      INTEGER :: n_coosig = 0
!>  Number of signals in a combination signal.
      INTEGER, DIMENSION(v3fit_max_combinations) :: n_sig_coosig = 0
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
      CHARACTER (len=data_short_name_length),                                  &
     &   DIMENSION(v3fit_max_combinations) :: coosig_units = ''

!  Initial Plasma Offsets
!>  Offset of plasma relative to the diagnostics.
      REAL (rprec) :: phi_offset = 0
!>  Vertical shift of plasma relative to the center.
      REAL (rprec) :: z_offset = 0

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
!>  -# @fixed_width{'sl'}  - straight line
!>  -# @fixed_width{'seg'} - segmented
!>  -# @fixed_width{'lm'}  - Levenberg-Marquardt,
!>  @see reconstruction
      CHARACTER (len=v3fit_step_name_length) :: step_type = 'sl'
!>  Cutoff value for singular values in matrix inverse.
      REAL (rprec)                           :: cut_inv_svd = 1.0E-10

!>  Zero out vmec xcdot array
!>  @todo FIXME: VMEC specific. Not implemented.
      LOGICAL :: l_zero_xcdot = .true.

!  Declare Namelist
      NAMELIST/v3fit_main_nli/                                                 &
     &   main_nli_filename,                                                    &
!  Equilibrium input files
     &   vmec_nli_filename, vmec_wout_input, vacuum_nli_filename,              &
!  Signal input files
     &   mdsig_list_filename, sxrch_dot_filename, thscte_dot_filename,         &
     &   ipch_dot_filename, mse_dot_filename, sxrem_ratio_dot_filename,        &
!  Array allocation sizes
     &   na_s_desc,                                                            &
!  Task specification and work variables
     &   my_task,                                                              &
!  Model profile specification
!  Electron Denisty
     &   pp_ne_ptype, pp_ne_b, pp_ne_as, pp_ne_af,                             &
!  Soft x-ray emissvity
     &   pp_sxrem_ptype, pp_sxrem_b, pp_sxrem_as, pp_sxrem_af,                 &
     &   num_sxrem_p, pp_sxrem_ptype_a, pp_sxrem_b_a, pp_sxrem_as_a,           &
     &   pp_sxrem_af_a,                                                        &
!  Electron temperature
     &   pp_te_ptype, pp_te_b, pp_te_as, pp_te_af,                             &
!  Soft x-ray ratio
     &   sxrem_te_a, sxrem_ratio_a,                                            &
!  Model Specification Variables
     &   model_ne_type, model_sxrem_type, model_sxrem_type_a,                  &
     &   model_te_type, ne_pp_unit, ne_min, te_min, e_pressure_fraction,       &
     &   emission_file,                                                        &
!  Reconstruction constraints
     &   n_rc, rc_type, rc_index, rc_value,                                    &
!  Derived parameters
     &   n_dp, dp_type, dp_index,                                              &
!  Reconstruction parameters
     &   n_rp, rp_type, rp_index, rp_index2, rp_vrnc, rp_range_type,           &
     &   rp_range_value, rp_range_index,                                       &
!  Signal Data specification
     &   n_sdata_o, iw_sdo_verbose, sdo_data_a, sdo_sigma_a,                   &
     &   sdo_weight_a, mag_a, mag_force_coil,                                  &
!  Magnetic integration parameters
     &   pol_rad_ratio,                                                        &
!  External current signal specification
!>  @todo FIXME: The label of s0 and u0 are vmec specific. These should be
!>  changed to a radial coordinate and poloidal angle.
     &   extcurz_s0, extcurz_u0,                                               &
!  Signal data sigma specification
     &   sdo_s_spec_imin, sdo_s_spec_imax, sdo_s_spec_floor,                   &
     &   sdo_s_spec_fraction,                                                  &
!  Signal data weight specification
     &   sdo_w_spec_imin, sdo_w_spec_imax, sdo_w_spec_weight,                  &
!  Magnetic signal specification
     &   mag_spec_imin, mag_spec_imax, mag_spec_use_induced,                   &
!  Geometric information
     &   r_major_radius, a_minor_radius,                                       &
!  Limiter iso function specification.
     &   n_lif, n_phi_lif, lif_arz, lif_rc, lif_zc, lif_sigma,                 &
     &   lif_phi_degree, lif_on_edge,                                          &
!  Prior signal specification.
     &   n_prior, prior_name, prior_param_name, prior_indices,                 &
     &   prior_units,                                                          &
!  Combinations of other signals specification.
     &   n_coosig, n_sig_coosig, coosig_indices, coosig_coeff,                 &
     &   coosig_type, coosig_name, coosig_units,                               &
!  Initial PLasma Offsets
     &   phi_offset, z_offset,                                                 &
!  Reconstruction step control specification.
     &   nrstep, dg2_stop, cut_svd, cut_eff, cut_marg_eff, cut_delta_a,        &
     &   cut_dg2, astep_max, step_type, cut_inv_svd,                           &
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
      CHARACTER (len=*), INTENT(in) :: namelist_file

!  local variables
      INTEGER                       :: iou_mnli, status
      INTEGER                       :: j, i
      INTEGER                       :: imin, imax
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Initalize a default value of the I\O unit. V3FIT increments from there.
      iou_mnli = 0
      CALL safe_open(iou_mnli, status, TRIM(namelist_file),                    &
     &               'old', 'formatted')
      CALL assert_eq(0, status, 'v3fit_input_read_namelist' //                 &
     &   ': Safe_open of ' // TRIM(namelist_file) // ' failed')

!  Read the namelist input file.
      READ (iou_mnli, nml=v3fit_main_nli)
      CLOSE (iou_mnli, iostat=status)
      CALL assert_eq(0, status, 'v3fit_input_read_namelist' //                 &
     &   ': Error losing ' // TRIM(namelist_file) // ' failed')

!  Toggle use coil signals for the magnetics from the specification
      DO j = 1, SIZE(mag_spec_imin)
         IF (mag_spec_imin(j) .le. 0) EXIT
         imin = MAX(MIN(mag_spec_imin(j), v3fit_max_diagnostics), 1)
         imax = MIN(MAX(mag_spec_imax(j), imin), v3fit_max_diagnostics)
         mag_a(imin:imax) = mag_spec_use_induced(j)
      END DO

!  Set weights from weight specification.
      DO j = 1, SIZE(sdo_w_spec_imin)
         IF (sdo_w_spec_imin(j) .le. 0) EXIT
         imin = MAX(MIN(sdo_w_spec_imin(j), v3fit_max_signals), 1)
         imax = MIN(MAX(sdo_w_spec_imax(j), imin), v3fit_max_signals)
         sdo_weight_a(imin:imax) = sdo_w_spec_weight(j)
      END DO

!  Set sigmas from the sigma specification.
      DO j = 1,SIZE(sdo_s_spec_imin)
         IF (sdo_s_spec_imin(j) .le. 0) EXIT
         imin = MAX(MIN(sdo_s_spec_imin(j), v3fit_max_signals), 1)
         imax = MIN(MAX(sdo_s_spec_imax(j), imin), v3fit_max_signals)
         DO i = imin, imax
            IF ((sdo_sigma_a(i) .eq. 0) .or.                                   &
     &          sdo_spec_can_overwrite) THEN
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


!-------------------------------------------------------------------------------
!  All the coding between these boxes can be removed once the non _a sxrem
!  values are removed. ***
!-------------------------------------------------------------------------------
!  Check the non _a sxrem profile variables. If pp_sxrem_ptype does not equal
!  'none' then over write the the first elements of the _a versions.
      SELECT CASE (TRIM(model_sxrem_type))

         CASE ('none') ! Do nothing

         CASE DEFAULT
            WRITE (*,*) ' *** Overwiting the sxrem profile first ' //          &
     &                  'element.'
            WRITE (*,*) '     Non array soft xray profiles have ' //           &
     &                  'been deprecated'
            model_sxrem_type_a(1) = model_sxrem_type
            pp_sxrem_ptype_a(1) = pp_sxrem_ptype
            pp_sxrem_b_a(1,:) = pp_sxrem_b
            pp_sxrem_as_a(1,:) = pp_sxrem_as
            pp_sxrem_af_a(1,:) = pp_sxrem_af

!  We'll assume that any reconstuction parameters that refer to the non _a also
!  are the first element.
            DO i = 1, n_rp
               SELECT CASE(TRIM(rp_type(i)))

                  CASE ('pp_sxrem_b')
                     WRITE (*,1000) 'pp_sxrem_b', 'pp_sxrem_b_a'
                     rp_type(i) = 'pp_sxrem_b_a'
                     rp_index2(i) = rp_index(i)
                     rp_index(i) = 1

                  CASE ('pp_sxrem_as')
                     WRITE (*,1000) 'pp_sxrem_as', 'pp_sxrem_as_a'
                     rp_type(i) = 'pp_sxrem_as_a'
                     rp_index2(i) = rp_index(i)
                     rp_index(i) = 1

                  CASE ('pp_sxrem_af')
                     WRITE (*,1000) 'pp_sxrem_af', 'pp_sxrem_af_a'
                     rp_type(i) = 'pp_sxrem_af_a'
                     rp_index2(i) = rp_index(i)
                     rp_index(i) = 1

               END SELECT

!  Need to update ranges as well starting with the lower range.
               SELECT CASE(TRIM(rp_range_type(i,1)))

                  CASE ('pp_sxrem_b')
                     WRITE (*,1001) 'pp_sxrem_b', 'pp_sxrem_b_a'
                     rp_range_type(i,1) = 'pp_sxrem_b_a'
                     rp_range_index(i,1,2) = rp_range_index(i,1,1)
                     rp_range_index(i,1,1) = 1

                  CASE ('pp_sxrem_as')
                     WRITE (*,1001) 'pp_sxrem_as', 'pp_sxrem_as_a'
                     rp_range_type(i,1) = 'pp_sxrem_as_a'
                     rp_range_index(i,1,2) = rp_range_index(i,1,1)
                     rp_range_index(i,1,1) = 1

                  CASE ('pp_sxrem_af')
                     WRITE (*,1001) 'pp_sxrem_af', 'pp_sxrem_af_a'
                     rp_range_type(i,1) = 'pp_sxrem_af_a'
                     rp_range_index(i,1,2) = rp_range_index(i,1,1)
                     rp_range_index(i,1,1) = 1

               END SELECT

!  Update the lower ranges.
               SELECT CASE(TRIM(rp_range_type(i,2)))

                  CASE ('pp_sxrem_b')
                     WRITE (*,1002) 'pp_sxrem_b', 'pp_sxrem_b_a'
                     rp_range_type(i,2) = 'pp_sxrem_b_a'
                     rp_range_index(i,2,2) = rp_range_index(i,2,1)
                     rp_range_index(i,2,1) = 1

                  CASE ('pp_sxrem_as')
                     WRITE (*,1002) 'pp_sxrem_as', 'pp_sxrem_as_a'
                     rp_range_type(i,2) = 'pp_sxrem_as_a'
                     rp_range_index(i,2,2) = rp_range_index(i,2,1)
                     rp_range_index(i,2,1) = 1

                  CASE ('pp_sxrem_af')
                     WRITE (*,1002) 'pp_sxrem_af', 'pp_sxrem_af_a'
                     rp_range_type(i,2) = 'pp_sxrem_af_a'
                     rp_range_index(i,2,2) = rp_range_index(i,2,1)
                     rp_range_index(i,2,1) = 1

               END SELECT
            END DO

!  Update priors to use _a arrays.
            DO i = 1, n_prior
               SELECT CASE(TRIM(prior_name(i)))
                  CASE ('pp_sxrem_b')
                     WRITE (*,1003) 'pp_sxrem_b', 'pp_sxrem_b_a'
                     prior_name(i) = 'pp_sxrem_b_a'
                     prior_indices(i,2) = prior_indices(i,1)

                  CASE ('pp_sxrem_as')
                     WRITE (*,1003) 'pp_sxrem_as', 'pp_sxrem_as_a'
                     prior_name(i) = 'pp_sxrem_as_a'
                     prior_indices(i,2) = prior_indices(i,1)

                  CASE ('pp_sxrem_af')
                     WRITE (*,1003) 'pp_sxrem_af', 'pp_sxrem_af_a'
                     prior_name(i) = 'pp_sxrem_af_a'
                     prior_indices(i,2) = prior_indices(i,1)
               END SELECT
            END DO
      END SELECT

1000  FORMAT('Changing reconstruction parameter ',a,' to ', a)
1001  FORMAT('Changing lower parameter constraint ',a,' to ', a)
1002  FORMAT('Changing upper parameter constraint ',a,' to ', a)
1003  FORMAT('Changing prior signal ',a,' to ', a)
!-------------------------------------------------------------------------------
!  All the coding between these boxes can be removed once the non _a sxrem
!  values are removed. ***
!-------------------------------------------------------------------------------

      WRITE (*,*) ' *** V3FIT namelist input read from ' //                    &
     &            TRIM(namelist_file)

      CALL profiler_set_stop_time('v3fit_input_read_namelist',                 &
     &                            start_time)

      END SUBROUTINE

      END MODULE
