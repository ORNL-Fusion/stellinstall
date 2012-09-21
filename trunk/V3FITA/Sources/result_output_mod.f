!*******************************************************************************
!  File result_output_mod.f
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  This modules containes code to output the results of a v3fit reconstruction 
!  into a computer readable netcdf file.
!
!*******************************************************************************

      MODULE result_output
      USE gsq_mod
      USE signal_T
      USE recon_param_T
      USE ezcdf
      USE stel_kinds, only : rprec
      USE v3f_vmec_comm, only : vvc_smaleli

!*******************************************************************************
!  Object representing the results file. This is mean to be an opaque. All 
!  fields private and can only be accessed from inside this module.
!
!  NetCDF variables:
!    ncid    Reference to the netcdf file.
!    status  Status value of the netcdf file. If this evaluates to anything
!            that nf_noerror, then an assert is triggered.
!
!  Dimensions: Reference IDs for the created dimensions.
!    maxnsteps_dim_id  The max number of reconstruction steps.
!    nparam_dim_id     Number of reconstruction parameters.
!    nsignal_dim_id    Number of signals.
!
!  Variables: Reference IDs for the created variables.
!    nsteps_id (scalar) Actual number of steps taken.
!    
!    g2_id (maxnsteps_dim_id) g^2 value for each reconstruction step.
!
!    vmec_rbc00_id (maxnsteps_dim_id) VMEC value of inital major radus
!                                      position.
!    vmec_rbc01_id (maxnsteps_dim_id) VMEC value of the inital oval radius in
!                                      the r direction.
!    vmec_zbs01_id (maxnsteps_dim_id) VMEC value of the inital oval radius in
!                                      the z direction.
!    vmec_i_sub_i_id (maxnsteps_dim_id) VMEC value of the internal inductance.
!
!    param_value_id (maxnsteps_dim_id, nparam_dim_id) Parameter values for each
!                   reconstruction step.
!    param_sigma_id (maxnsteps_dim_id, nparam_dim_id) Parameter sigmas for each
!                   reconstruction step.
!
!    signal_obser_value_id (maxnsteps_dim_id, nsignal_dim_id) Observed signals 
!                          for each reconstruction step.
!    signal_model_value_id (maxnsteps_dim_id, nsignal_dim_id) Modeled signals 
!                          for each reconstruction step.
!    signal_sigma_id (maxnsteps_dim_id, nsignal_dim_id) Signal sigmas 
!                          for each reconstruction step.
!    signal_eff_matrix_id (maxnsteps_dim_id, nparam_dim_id, nsignal_dim_id)
!                          Signal effectiveness matrix.
!
!*******************************************************************************
      TYPE result_file
         PRIVATE
         INTEGER :: ncid, status
! Dimension IDs
         INTEGER :: maxnsteps_dim_id, nparam_dim_id, nsignal_dim_id
! Variable IDs
         INTEGER :: vmec_rbc00_id, vmec_rbc01_id, vmec_zbs01_id,               &
     &              vmec_i_sub_i_id
         INTEGER :: nsteps_id
         INTEGER :: g2_id
         INTEGER :: param_value_id, param_sigma_id
         INTEGER :: signal_obser_value_id, signal_model_value_id,              &
     &              signal_sigma_id, signal_eff_matrix_id
      END TYPE result_file

      PUBLIC  :: result_output_file_create, result_output_file_write_step,     &
     &           result_output_file_close
      PRIVATE :: result_output_check_error
      
      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Create and setup a new netcdf file.
!
!  A discription of the resulting NetCDF file is provided.
!    
!  Dimensions:
!    nextcur    Number of external currents.
!    maxnsteps  Maximum number of steps.
!    nparam     Number of parameters.
!    nsignal    Number of signals.
!    string_len Max length of strings. Set to 80.
!
!  Variables:
!    nsteps (INT) Actual number of reconstruction steps. This starts from zero.
!
!    g2 (DOUBLE) g^2 of each reconstruction step.
!
!    vmec_extcur External currents.
!    vmec_curtor Non-reconstructed plasma current.
!    vmec_pres_scale Non-reconstructed plasma pressure.
!
!    vmec_rbc00   Inital major radus position.
!    vmec_rbc01   Inital oval radius r direction.
!    vmec_zbs01   Inital oval radius z direction.
!    vmec_i_sub_i Internal inducance.
!
!    param_name (STRING) Name of the reconstruction parameter.
!    param_index (INT) Index of the reconstruction parameter.
!    param_value (DOUBLE) Value of the reconstruction parameter at each 
!                         reconstruction step.
!    param_sigma (DOUBLE) Sigma of the reconstruction parameter at each 
!                         reconstruction step.
!
!    signal_name (STRING) Long name of the signal.
!    signal_type (STRING) Description of the signal.
!    signal_observed_value (DOUBLE) Observed value of the signal at each 
!                                   reconstruction step.
!    signal_model_value (DOUBLE) Modeled value of the signal at each 
!                                reconstruction step.
!    signal_sigma (DOUBLE) Sigma value of the reconstruction step.
!                          sig = Sqrt( Model_sig^2 + Observe_sig^2 )
!    signal_eff_matrix (DOUBLE) Matrix of signal effectiveness values.
!
!*******************************************************************************
      FUNCTION result_output_file_create(path,                                 &
     &                                   params,                               &
     &                                   sd_model,                             &
     &                                   sd_observed)
      USE v3f_global, ONLY: nrstep, n_rparam, n_sdata_o, s_desc, model_a
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER (len = *), INTENT(in) :: path
      TYPE (recon_param), DIMENSION(n_rparam), INTENT(in) :: params
      TYPE (signal_data), DIMENSION(:), INTENT (in) ::                         &
     &   sd_model, sd_observed
      TYPE (result_file) :: result_output_file_create
      INTEGER :: i, beginIndex, endIndex

! Netcdf ids that are not needed after this subroutine is finished.
! Dimension IDs
      INTEGER :: string_len_id, nextcur_dim_id
! Variable IDs
      INTEGER :: vmec_extcur_id, vmec_curtor_id, vmec_pres_scale_id
      INTEGER :: param_name_id, param_index_id
      INTEGER :: signal_name_id, signal_type_id

! Working variables
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: sem_a
      REAL(rprec), ALLOCATABLE, DIMENSION(:)   :: psigma_v
      INTEGER                                  :: err
      CHARACTER(len=*), PARAMETER              :: sub_name =                   &
     &  'result_output_file_create: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Find the start and end index of the basename. The basename is defined as
! input.BASENAME.v3fit if the prefix "input" and postfix "v3fit" are present.
! Other wise the entire filename is used.
      beginIndex = INDEX(path, "input.") + 6
      endIndex = INDEX(path, ".v3fit", .TRUE.) - 1

      result_output_file_create%status =                                       &
     &   nf_create("result." //                                                &
     &             TRIM(path(beginIndex:endIndex)) //                          &
     &             ".nc",                                                      &
     &             nf_clobber,                                                 &
     &             result_output_file_create%ncid)
      CALL result_output_check_error(result_output_file_create%status)

! Define dimensions for data.
      result_output_file_create%status =                                       &
     &   nf_def_dim(result_output_file_create%ncid,                            &
     &              "nextcur",                                                 &
     &              model_a%eqstate%fixp%nextcur,                              &
     &              nextcur_dim_id)
      CALL result_output_check_error(result_output_file_create%status)

! Add space dimension for all reconstruction steps including a step before
! reconstruction.
      result_output_file_create%status =                                       &
     &   nf_def_dim(result_output_file_create%ncid,                            &
     &              "maxnsteps",                                               &
     &              nrstep + 1,                                                &
     &              result_output_file_create%maxnsteps_dim_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_dim(result_output_file_create%ncid,                            &
     &              "nparam",                                                  &
     &              n_rparam,                                                  &
     &              result_output_file_create%nparam_dim_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_dim(result_output_file_create%ncid,                            &
     &              "nsignal",                                                 &
     &              n_sdata_o,                                                 &
     &              result_output_file_create%nsignal_dim_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_dim(result_output_file_create%ncid,                            &
     &              "string_len",                                              &
     &              80,                                                        &
     &              string_len_id)
      CALL result_output_check_error(result_output_file_create%status)

! Define variables for data.
      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "nsteps",                                                  &
     &              nf_int,                                                    &
     &              0, 0,                                                      &
     &              result_output_file_create%nsteps_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "g2",                                                      &
     &              nf_double,                                                 &
     &              1,                                                         &
     &              (/ result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%g2_id)
      CALL result_output_check_error(result_output_file_create%status)

! Define vmec variables.
      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "vmec_extcur",                                             &
     &              nf_double,                                                 &
     &              1,                                                         &
     &              (/ nextcur_dim_id /),                                      &
     &              vmec_extcur_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "vmec_curtor",                                             &
     &              nf_double,                                                 &
     &              0, 0,                                                      &
     &              vmec_curtor_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "vmec_pres_scale",                                         &
     &              nf_double,                                                 &
     &              0, 0,                                                      &
     &              vmec_pres_scale_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "vmec_rbc00",                                              &
     &              nf_double,                                                 &
     &              1,                                                         &
     &              (/ result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%vmec_rbc00_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "vmec_rbc01",                                              &
     &              nf_double,                                                 &
     &              1,                                                         &
     &              (/ result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%vmec_rbc01_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "vmec_zbs01",                                              &
     &              nf_double,                                                 &
     &              1,                                                         &
     &              (/ result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%vmec_zbs01_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "vmec_i_sub_i",                                            &
     &              nf_double,                                                 &
     &              1,                                                         &
     &              (/ result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%vmec_i_sub_i_id)
      CALL result_output_check_error(result_output_file_create%status)

! Define reconstruction parameter variables.
      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "param_name",                                              &
     &              nf_char,                                                   &
     &              2,                                                         &
     &              (/ string_len_id,                                          &
     &                 result_output_file_create%nparam_dim_id /),             &
     &              param_name_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "param_index",                                             &
     &              nf_int,                                                    &
     &              1,                                                         &
     &              (/ result_output_file_create%nparam_dim_id /),             &
     &              param_index_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "param_value",                                             &
     &              nf_double,                                                 &
     &              2,                                                         &
     &              (/ result_output_file_create%nparam_dim_id,                &
     &                 result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%param_value_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "param_sigma",                                             &
     &              nf_double,                                                 &
     &              2,                                                         &
     &              (/ result_output_file_create%nparam_dim_id,                &
     &                 result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%param_sigma_id)
      CALL result_output_check_error(result_output_file_create%status)

! Define signal variables.
      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "signal_name",                                             &
     &              nf_char,                                                   &
     &              2,                                                         &
     &              (/ string_len_id,                                          &
     &                 result_output_file_create%nsignal_dim_id /),            &
     &              signal_name_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "signal_type",                                             &
     &              nf_char,                                                   &
     &              2,                                                         &
     &              (/ string_len_id,                                          &
     &                 result_output_file_create%nsignal_dim_id /),            &
     &              signal_type_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "signal_observed_value",                                   &
     &              nf_double,                                                 &
     &              2,                                                         &
     &              (/ result_output_file_create%nsignal_dim_id,               &
     &                 result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%signal_obser_value_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "signal_model_value",                                      &
     &              nf_double,                                                 &
     &              2,                                                         &
     &              (/ result_output_file_create%nsignal_dim_id,               &
     &                 result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%signal_model_value_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "signal_sigma",                                            &
     &              nf_double,                                                 &
     &              2,                                                         &
     &              (/ result_output_file_create%nsignal_dim_id,               &
     &                 result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%signal_sigma_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_def_var(result_output_file_create%ncid,                            &
     &              "signal_eff_matrix",                                       &
     &              nf_double,                                                 &
     &              3,                                                         &
     &              (/ result_output_file_create%nsignal_dim_id,               &
     &                 result_output_file_create%nparam_dim_id,                &
     &                 result_output_file_create%maxnsteps_dim_id /),          &
     &              result_output_file_create%signal_eff_matrix_id)
      CALL result_output_check_error(result_output_file_create%status)

! Finished defining netcdf file. Exit out of define mode.
      result_output_file_create%status =                                       &
     &   nf_enddef(result_output_file_create%ncid)
      CALL result_output_check_error(result_output_file_create%status)

! Fill in initial variable values.
      result_output_file_create%status =                                       &
     &   nf_put_var1_int(result_output_file_create%ncid,                       &
     &                   result_output_file_create%nsteps_id,                  &
     &                   1, 0)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_put_var1_double(result_output_file_create%ncid,                    &
     &                      result_output_file_create%g2_id,                   & 
     &                      1, gsq_g2)
      CALL result_output_check_error(result_output_file_create%status)

! Fill in vmec variables.
! Only the number of currents specified in the mgrid file should be written out.
      result_output_file_create%status =                                       &
     &   nf_put_var_double(result_output_file_create%ncid,                     &
     &                     vmec_extcur_id,                                     &
     &                     model_a%eqstate%varp%extcur(                        &
     &                       1:model_a%eqstate%fixp%nextcur                    &
     &                     ))
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_put_var1_double(result_output_file_create%ncid,                    &
     &                      vmec_curtor_id,                                    &
     &                      1, model_a%eqstate%varp%curtor)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_put_var1_double(result_output_file_create%ncid,                    &
     &                      vmec_pres_scale_id,                                &
     &                      1, model_a%eqstate%varp%pres_scale)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_put_var1_double(result_output_file_create%ncid,                    &
     &                      result_output_file_create%vmec_rbc00_id,           &
     &                      1, model_a%eqstate%varp%rbc(0,0))
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_put_var1_double(result_output_file_create%ncid,                    &
     &                      result_output_file_create%vmec_rbc01_id,           &
     &                      1, model_a%eqstate%varp%rbc(0,1))
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_put_var1_double(result_output_file_create%ncid,                    &
     &                      result_output_file_create%vmec_zbs01_id,           &
     &                      1, model_a%eqstate%varp%zbs(0,1))
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf_put_var1_double(result_output_file_create%ncid,                    &
     &                      result_output_file_create%vmec_i_sub_i_id,         &
     &                      1, vvc_smaleli)
      CALL result_output_check_error(result_output_file_create%status)

! Fill in reconstuction parameters.
      DO i = 1, n_rparam
         result_output_file_create%status =                                    &
     &      nf_put_vara_text(result_output_file_create%ncid,                   &
     &                       param_name_id,                                    &
     &                       (/ 1, i /), (/ 20, 1 /),                          &
     &                       params(i)%p_type)
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf_put_var1_int(result_output_file_create%ncid,                    &
     &                      param_index_id,                                    &
     &                      i,                                                 &
     &                      params(i)%index)
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf_put_var1_double(result_output_file_create%ncid,                 &
     &                         result_output_file_create%param_value_id,       &
     &                         (/ i, 1 /),                                     &
     &                         params(i)%value)
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

! There is no sigma to be calculated for the first reconstruction step. Fill
! this step in with zeros.
         result_output_file_create%status =                                    &
     &      nf_put_var1_double(result_output_file_create%ncid,                 &
     &                         result_output_file_create%param_sigma_id,       &
     &                         (/ i, 1 /),                                     &
     &                         0.0)
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)
      END DO

! Fill in signals
      DO i = 1, n_sdata_o
         result_output_file_create%status =                                    &
     &      nf_put_vara_text(result_output_file_create%ncid,                   &
     &                       signal_name_id,                                   &
     &                       (/ 1, i /), (/ 30, 1 /),                          &
     &                       s_desc(i)%s_name)
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf_put_vara_text(result_output_file_create%ncid,                   &
     &                       signal_type_id,                                   &
     &                       (/ 1, i /), (/ 10, 1 /),                          &
     &                       s_desc(i)%s_type)
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf_put_var1_double(                                                &
     &         result_output_file_create%ncid,                                 &
     &         result_output_file_create%signal_obser_value_id,                &
     &         (/ i, 1 /),                                                     &
     &         sd_observed(i)%data(1))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf_put_var1_double(                                                &
     &         result_output_file_create%ncid,                                 &
     &         result_output_file_create%signal_model_value_id,                &
     &         (/ i, 1 /),                                                     &
     &         sd_model(i)%data(1))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf_put_var1_double(                                                &
     &         result_output_file_create%ncid,                                 &
     &         result_output_file_create%signal_sigma_id,                      &
     &         (/ i, 1 /),                                                     &
     &         SQRT(sd_model(i)%sigma(1)**2 +                                  &
     &              sd_observed(i)%sigma(1)**2))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)
      END DO

! Calculating the signal effective causes a crash on the D3D test case. For this
! initial reconstruction step, just fill in zeros.
      ALLOCATE(sem_a(gsq_nrp,gsq_nsig),STAT=err)
      CALL assert_eq(0,err,sub_name // 'allocate sem_a')
      result_output_file_create%status =                                       &
     &   nf_put_vara_double(                                                   &
     &      result_output_file_create%ncid,                                    &
     &      result_output_file_create%signal_eff_matrix_id,                    &
     &      (/ 1, 1, 1 /), (/ n_sdata_o, n_rparam, 1 /),                       &
     &      TRANSPOSE(sem_a))
      DEALLOCATE(sem_a)

      END FUNCTION result_output_file_create

!*******************************************************************************
!  Write reconstruction step to the netcdf file.
!*******************************************************************************
      SUBROUTINE result_output_file_write_step(result_file_ref,                &
     &                                         params,                         &
     &                                         sd_model,                       &
     &                                         sd_observed)
      USE v3f_global, ONLY: n_rparam, n_sdata_o, model_a
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE (result_file) :: result_file_ref
      TYPE (recon_param), DIMENSION(n_rparam), INTENT(in) :: params
      TYPE (signal_data), DIMENSION(:), INTENT (in) ::                         &
     &   sd_model, sd_observed
      INTEGER :: nsteps, i

! Working variables
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: gsq_1pcl, sem_a
      REAL(rprec), ALLOCATABLE, DIMENSION(:)   :: psigma_v
      INTEGER                                  :: err
      CHARACTER(len=*), PARAMETER              :: sub_name =                   &
     &  'result_output_file_write_step: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Read and update current step.
      result_file_ref%status =                                                 &
     &   nf_get_var1_int(result_file_ref%ncid,                                 &
     &                   result_file_ref%nsteps_id,                            &
     &                   1, nsteps)
      CALL result_output_check_error(result_file_ref%status)
      nsteps = nsteps + 1
      result_file_ref%status =                                                 &
     &   nf_put_var1_int(result_file_ref%ncid,                                 &
     &                   result_file_ref%nsteps_id,                            &
     &                   1, nsteps)
      CALL result_output_check_error(result_file_ref%status)

! Steps start at zero however array offsets in the netcdf routines start at one. 
! Increment nsteps to reflect the correct array position.
      nsteps = nsteps + 1

      result_file_ref%status =                                                 &
     &   nf_put_var1_double(result_file_ref%ncid,                              &
     &                      result_file_ref%g2_id,                             &
     &                      nsteps, gsq_g2)
      CALL result_output_check_error(result_file_ref%status)

! Fill in vmec parameters.
      result_file_ref%status =                                                 &
     &   nf_put_var1_double(result_file_ref%ncid,                              &
     &                      result_file_ref%vmec_rbc00_id,                     &
     &                      nsteps, model_a%eqstate%varp%rbc(0,0))
      CALL result_output_check_error(result_file_ref%status)

      result_file_ref%status =                                                 &
     &   nf_put_var1_double(result_file_ref%ncid,                              &
     &                      result_file_ref%vmec_rbc01_id,                     &
     &                      nsteps, model_a%eqstate%varp%rbc(0,1))
      CALL result_output_check_error(result_file_ref%status)

      result_file_ref%status =                                                 &
     &   nf_put_var1_double(result_file_ref%ncid,                              &
     &                      result_file_ref%vmec_zbs01_id,                     &
     &                      nsteps, model_a%eqstate%varp%zbs(0,1))
      CALL result_output_check_error(result_file_ref%status)

      result_file_ref%status =                                                 &
     &   nf_put_var1_double(result_file_ref%ncid,                              &
     &                      result_file_ref%vmec_i_sub_i_id,                   &
     &                      nsteps, vvc_smaleli)
      CALL result_output_check_error(result_file_ref%status)

! Calculate parameter sigma.
      ALLOCATE(gsq_1pcl(gsq_nrp,gsq_nsv),STAT=err)
      CALL assert_eq(0,err,sub_name // 'allocate gsq_1pcl')
      CALL gsq_ps1(gsq_1pcl)
! Fill in reconstuction parameters.
      DO i = 1, n_rparam
         result_file_ref%status =                                              &
     &      nf_put_var1_double(result_file_ref%ncid,                           &
     &                         result_file_ref%param_value_id,                 &
     &                         (/ i, nsteps /), params(i)%value)
         CALL result_output_check_error(result_file_ref%status)

         result_file_ref%status =                                              &
     &      nf_put_var1_double(result_file_ref%ncid,                           &
     &                         result_file_ref%param_sigma_id,                 &
     &                         (/ i, nsteps /), gsq_1pcl(i, gsq_nrp))
         CALL result_output_check_error(result_file_ref%status)
      END DO
      DEALLOCATE(gsq_1pcl)

! Fill in signals
      DO i = 1, n_sdata_o
         result_file_ref%status =                                              &
     &      nf_put_var1_double(                                                &
     &         result_file_ref%ncid,                                           &
     &         result_file_ref%signal_obser_value_id,                          &
     &         (/ i, nsteps /), sd_observed(i)%data(1))
         CALL result_output_check_error(result_file_ref%status)

         result_file_ref%status =                                              &
     &      nf_put_var1_double(                                                &
     &         result_file_ref%ncid,                                           &
     &         result_file_ref%signal_model_value_id,                          &
     &         (/ i, nsteps /), sd_model(i)%data(1))
         CALL result_output_check_error(result_file_ref%status)

         result_file_ref%status =                                              &
     &      nf_put_var1_double(
     &         result_file_ref%ncid,                                           &
     &         result_file_ref%signal_sigma_id,                                &
     &         (/ i, nsteps /),                                                &
     &         SQRT(sd_model(i)%sigma(1)**2 +                                  &
     &              sd_observed(i)%sigma(1)**2))
         CALL result_output_check_error(result_file_ref%status)
      END DO

! Calculate signal effectiveness matrix
      ALLOCATE(psigma_v(gsq_nrp),STAT=err)
      CALL assert_eq(0,err,sub_name // 'allocate psigma_v')
      ALLOCATE(sem_a(gsq_nrp,gsq_nsig),STAT=err)
      CALL assert_eq(0,err,sub_name // 'allocate sem_a')
      CALL gsq_sem(gsq_ssigma,gsq_jacobian,psigma_v,sem_a)
      result_file_ref%status = nf_put_vara_double(                             &
     &   result_file_ref%ncid,                                                 &
     &   result_file_ref%signal_eff_matrix_id,                                 &
     &   (/ 1, 1, nsteps /), (/ n_sdata_o, n_rparam, 1 /),                     &
     &   TRANSPOSE(sem_a))
      DEALLOCATE(sem_a, psigma_v)

      END SUBROUTINE

!*******************************************************************************
!  Close the netcdf file.
!*******************************************************************************
      SUBROUTINE result_output_file_close(result_file_ref)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE (result_file) :: result_file_ref

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      result_file_ref%status = nf_close(result_file_ref%ncid)
      CALL result_output_check_error(result_file_ref%status)

      END SUBROUTINE

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
! Check and handle errors.
!*******************************************************************************
      SUBROUTINE result_output_check_error(status)
      USE v3_utilities                                  ! subroutine assert_eq
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: status

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      END SUBROUTINE

      END MODULE result_output
