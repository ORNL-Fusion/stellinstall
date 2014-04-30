!-------------------------------------------------------------------------------
!  The @header, @table_section, @table_subsection, @item and @end_table commands
!  are custom defined commands in Doxygen.in. They are defined under ALIASES.
!  For the page created here, the 80 column limit is exceeded. Arguments of
!  aliases are separated by ','. If you intended ',' to be a string you must use
!  an escaped comma '\,'.
!
!>  @page vmec_equilibrium_sec VMEC Equilibrium Manual
!>
!>  @tableofcontents
!>  @section vmec_equilibrium_intro_sec Introduction
!>  This page documents the V3FIT interface to the VMEC equilibrium solver. This
!>  lists all the parameters associated with a VMEC equilibrium. These
!>  parameters are divided input two types, reconstruction parameters and
!>  derived parameters. All parameters are documented in a table of the
!>  following form.
!>  @header{Input variable, Description, Code Reference}
!>
!>  @section vmec_equilibrium_recon_param_sec VMEC Reconstruction Parameters
!>  Reconstruction parameters are divided into two types, VMEC reconstruction
!>  parameters and auxillary parameters.
!>  @subsection vmec_equilibrium_vmec_recon_param_sec VMEC Equilibrium Reconstruction Parameters
!>  VMEC Equilibrium reconstruction parameters are parameters that maybe
!>  reconstructed that change the equilibirum. These parameters require
!>  the equilibrium to reconverge.
!>  @begin_table
!>     @item{ac,          1D Array of function current profile parameters.,     @ref vmec_input::ac}
!>     @item{ac_aux_s,    1D Array of segment current profile s poitions.,      @ref vmec_input::ac_aux_s}
!>     @item{ac_aux_f,    1D Array of segment current profile f values.,        @ref vmec_input::ac_aux_f}
!>     @item{ai,          1D Array of function current profile parameters.,     @ref vmec_input::ai}
!>     @item{ai_aux_s,    1D Array of segment current profile s poitions.,      @ref vmec_input::ai_aux_s}
!>     @item{ai_aux_f,    1D Array of segment current profile f values.,        @ref vmec_input::ai_aux_f}
!>     @item{am,          1D Array of function current profile parameters.,     @ref vmec_input::am}
!>     @item{am_aux_s,    1D Array of segment current profile s poitions.,      @ref vmec_input::am_aux_s}
!>     @item{am_aux_f,    1D Array of segment current profile f values.,        @ref vmec_input::am_aux_f}
!>     @item{bloat,       Extends the field into the vacuum region.,            @ref vmec_input::bloat}
!>     @item{rbc,         2D Array of outer boundary R cosine coeffiecents.,    @ref vmec_input::rbc}
!>     @item{zbs,         2D Array of outer boundary Z sine coeffiecents.,      @ref vmec_input::zbs}
!>     @item{rbs,         2D Array of outer boundary R sine coeffiecents.,      @ref vmec_input::rbs}
!>     @item{zbc,         2D Array of outer boundary Z cosine coeffiecents.,    @ref vmec_input::zbc}
!>     @item{extcur,      1D Array of external currents.,                       @ref vmec_input::extcur}
!>     @item{phiedge,     The edge magnetic flux.,                              @ref vmec_input::phiedge}
!>     @item{pres_scale,  The scaling parameter for the pressure.,              @ref vmec_input::pres_scale}
!>  @end_table
!>  @subsection vmec_equilibrium_aux_recon_param_sec VMEC Auxiliary Reconstruction Parameters
!>  VMEC Auxiliary reconstruction parameters are parameters that maybe
!>  reconstructed that do not change the equilibirum. These parameters do not
!>  require the equilibrium to reconverge.
!>  @begin_table
!>     @item{pp_ne_b,       1D Array of function density profile parameters.,     @ref vmec_class::ne::b}
!>     @item{pp_ne_as,      1D Array of segment density profile s poitions.,      @ref vmec_class::ne::as}
!>     @item{pp_ne_af,      1D Array of segment density profile f values.,        @ref vmec_class::ne::af}
!>     @item{pp_sxrem_b,    1D Array of function sxrem profile parameters.
!>                          DEPREICATED only use id when using the old profile
!>                          specification.,                                       @ref vmec_class::sxrem::b}
!>     @item{pp_sxrem_as,   1D Array of segment sxrem profile s poitions.
!>                          DEPREICATED only use id when using the old profile
!>                          specification.,                                       @ref vmec_class::sxrem::as}
!>     @item{pp_sxrem_af,   1D Array of segment sxrem profile f values.
!>                          DEPREICATED only use id when using the old profile
!>                          specification.,                                       @ref vmec_class::sxrem::af}
!>     @item{pp_sxrem_b_a,  2D Array of function sxrem profile parameters.
!>                          First index is the profile number.,                   @ref vmec_class::sxrem::b}
!>     @item{pp_sxrem_as_a, 2D Array of segment sxrem profile s poitions.
!>                          First index is the profile number.,                   @ref vmec_class::sxrem::as}
!>     @item{pp_sxrem_af_a, 2D Array of segment sxrem profile f values.
!>                          First index is the profile number.,                   @ref vmec_class::sxrem::af}
!>     @item{pp_te_b,       1D Array of function temperature profile parameters., @ref vmec_class::te::b}
!>     @item{pp_te_as,      1D Array of segment temperature profile s poitions.,  @ref vmec_class::te::as}
!>     @item{pp_te_af,      1D Array of segment temperature profile f values.,    @ref vmec_class::te::af}
!>     @item{phi_offset,    Phi angle offset in radians,                          @ref vmec_equilibrium::vmec_class::phi_offset}
!>  @end_table
!>
!>  @section vmec_equilibrium_derive_param_sec VMEC Derived Parameters
!>  Derived parameters are parameters that come from the equilibrium but are not
!>  reconstructable.
!>  @begin_table
!>     @item{rmnc,        2D Array of surface R cosine coeffiecents.,      @ref read_wout_mod::rmnc}
!>     @item{zmns,        2D Array of surface Z sine coeffiecents.,        @ref read_wout_mod::zmns}
!>     @item{lmns,        2D Array of surface lambda sine coeffiecents.,   @ref read_wout_mod::lmns}
!>     @item{gmnc,        2D Array of surface g cosine coeffiecents.,      @ref read_wout_mod::gmnc}
!>     @item{bsubsmns,    2D Array of surface B_s sine coeffiecents.,      @ref read_wout_mod::bsubsmns}
!>     @item{bsubumnc,    2D Array of surface B_u cosine coeffiecents.,    @ref read_wout_mod::bsubumnc}
!>     @item{bsubvmnc,    2D Array of surface B_v cosine coeffiecents.,    @ref read_wout_mod::bsubvmnc}
!>     @item{bsupumnc,    2D Array of surface B^u cosine coeffiecents.,    @ref read_wout_mod::bsupumnc}
!>     @item{bsupvmnc,    2D Array of surface B^v cosine coeffiecents.,    @ref read_wout_mod::bsupvmnc}
!>     @item{rmns,        2D Array of surface R sine coeffiecents.,        @ref read_wout_mod::rmns}
!>     @item{zmnc,        2D Array of surface Z cosine coeffiecents.,      @ref read_wout_mod::zmnc}
!>     @item{lmnc,        2D Array of surface lambda cosine coeffiecents., @ref read_wout_mod::lmnc}
!>     @item{gmns,        2D Array of surface g sine coeffiecents.,        @ref read_wout_mod::gmns}
!>     @item{bsubsmnc,    2D Array of surface B_s cosine coeffiecents.,    @ref read_wout_mod::bsubsmnc}
!>     @item{bsubumns,    2D Array of surface B_u sine coeffiecents.,      @ref read_wout_mod::bsubumns}
!>     @item{bsubvmns,    2D Array of surface B_v sine coeffiecents.,      @ref read_wout_mod::bsubvmns}
!>     @item{bsupumns,    2D Array of surface B^u sine coeffiecents.,      @ref read_wout_mod::bsupumns}
!>     @item{bsupvmns,    2D Array of surface B^v sine coeffiecents.,      @ref read_wout_mod::bsupvmns}
!>     @item{phi,         1D Array of Toroidal flux.,                      @ref read_wout_mod::phi}
!>     @item{iotaf,       1D Array of full grid iotabar.,                  @ref read_wout_mod::iotaf}
!>     @item{iotas,       1D Array of half grid iotabar.,                  @ref read_wout_mod::iotas}
!>     @item{vvc_smaleli, Internal inducance. Computed in @ref eqfor.,     @ref v3f_vmec_comm::vvc_smaleli}
!>     @item{vvc_kappa_p, Mean elongation. Computed in @ref eqfor.,        @ref v3f_vmec_comm::vvc_kappa_p}
!>     @item{betatot,     Total beta.,                                     @ref read_wout_mod::betatot}
!>     @item{betapol,     Poloidal beta.,                                  @ref read_wout_mod::betapol}
!>     @item{betator,     Toroidal beta.,                                  @ref read_wout_mod::betator}
!>     @item{betaxis,     Beta on axis.,                                   @ref read_wout_mod::betaxis}
!>  @end_table
!*******************************************************************************
!>  @file vmec_equilibrium.f
!>  @brief Contains module @ref vmec_equilibrium.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref vmec_class. This module contains
!>  all the code necessary to interface V3FIT with VMEC.
!>  @par Super Class:
!>  @ref equilibrium
!*******************************************************************************

      MODULE vmec_equilibrium
      USE stel_kinds
      USE data_parameters
      USE pprofile_T
      USE wout_file_cache
      USE threed1_file_cache

      IMPLICIT NONE

!*******************************************************************************
!  vmec equilibrium module parameters
!*******************************************************************************
!  Define id's only for values that can change. All others shouldn't be needed
!  outside of the vmec interface. Parameter id's start at 3 since 0-2 are
!  reserved for non equilibrium model parameters. These numbers will need to be
!  updated if any new model parameters are added.

!  Variable Parameters
!>  1D Array of function current profile parameters. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_id          =  3
!>  1D Array of segment current profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_aux_s_id    =  4
!>  1D Array of segment current profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_aux_f_id    =  5
!>  1D Array of function iotabar profile parameters. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_id          =  6
!>  1D Array of segment iotabar profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_aux_s_id    =  7
!>  1D Array of segment iotabar profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_aux_f_id    =  8
!>  1D Array of function pressure profile parameters. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_id          =  9
!>  1D Array of segment pressure profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_aux_s_id    = 10
!>  1D Array of segment pressure profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_aux_f_id    = 11
!>  Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_bloat_id       = 12
!>  2D Array of outer boundary R cosine coeffiecents. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_rbc_id         = 13
!>  2D Array of outer boundary Z sine coeffiecents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_zbs_id         = 14
!>  2D Array of outer boundary R sine coeffiecents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_rbs_id         = 15
!>  2D Array of outer boundary Z cosine coeffiecents. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_zbc_id         = 16
!>  1D Array of external currents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_extcur_id      = 17
!>  The toroidal current. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_curtor_id      = 18
!>  The edge magnetic flux. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_phiedge_id     = 19
!>  The scaling parameter for the pressure. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_pres_scale_id  = 20

!  Derived Parameters
!>  2D Array of surface R cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_rmnc_id        = 21
!>  2D Array of surface Z sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_zmns_id        = 22
!>  2D Array of surface lambda cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_lmns_id        = 23
!>  2D Array of surface g cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_gmnc_id        = 24
!>  2D Array of surface B_s sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubsmns_id    = 25
!>  2D Array of surface B_u cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubumnc_id    = 26
!>  2D Array of surface B_v cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubvmnc_id    = 27
!>  2D Array of surface B^u cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupumnc_id    = 28
!>  2D Array of surface B^v cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupvmnc_id    = 29
!>  2D Array of surface R sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_rmns_id        = 30
!>  2D Array of surface Z cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_zmnc_id        = 31
!>  2D Array of surface lambda sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_lmnc_id        = 32
!>  2D Array of surface g sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_gmns_id        = 33
!>  2D Array of surface B_s cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubsmnc_id    = 34
!>  2D Array of surface B_u sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubumns_id    = 35
!>  2D Array of surface B_v sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubvmns_id    = 36
!>  2D Array of surface B^u sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupumns_id    = 37
!>  2D Array of surface B^v sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupvmns_id    = 38
!>  Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_phi_id         = 39
!>  1D Array of full grid iotabar. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_iotaf_id       = 40
!>  1D Array of half grid iotabar. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_iotas_id       = 41

!  Auxiliary Model Parameters
!>  1D Array of function density profile parameters. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_ne_b_id     = 42
!>  1D Array of segment density profile s poitions. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_ne_as_id    = 43
!>  1D Array of segment density profile f values. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_ne_af_id    = 44
!>  1D Array of function sxrem profile parameters. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_sxrem_b_id  = 45
!>  1D Array of segment sxrem profile s poitions. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_sxrem_as_id = 46
!>  1D Array of segment sxrem profile f values. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_sxrem_af_id = 47
!>  1D Array of function temperature profile parameters. Defined in
!>  @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_te_b_id     = 48
!>  1D Array of segment temperature profile s poitions. Defined in
!>  @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_te_as_id    = 49
!>  1D Array of segment temperature profile f values. Defined in
!>  @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_te_af_id    = 50

!  Derived Parameters
!>  Internal inducance. Defined in @ref v3f_vmec_comm. Computed in @ref eqfor.
      INTEGER, PARAMETER :: vmec_vvc_smaleli_id = 51
!>  Mean elongation. Defined in @ref v3f_vmec_comm. Computed in @ref eqfor.
      INTEGER, PARAMETER :: vmec_vvc_kappa_p_id = 52
!>  Total Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betatot_id     = 53
!>  Poloidal Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betapol_id     = 54
!>  Toroidal Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betator_id     = 55
!>  Beta on axis. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betaxis_id     = 56

!>  Plasma Phi offset. This is a parameter to allow changing the phi angle of a
!>  quasi helical state in an RFP. Defined in @ref vmec_equilibrium.
      INTEGER, PARAMETER :: vmec_phi_offset_id  = 57

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) vmec magnetic cache
!  2) vmec base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Object to cache data needed to compute magnetic signals.
!-------------------------------------------------------------------------------
      TYPE vmec_magnetic_cache
!>  Volume differential elements
         REAL (rprec)                            :: ds
         REAL (rprec)                            :: du
         REAL (rprec)                            :: dv

!>  R and Z positions on the suv gird points.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: rsuv => null()
         REAL (rprec), DIMENSION(:,:,:), POINTER :: zsuv => null()

!>  The following have a 1/Sqrt(g) factor embeded.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jrsuv => null()
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jphisuv => null()
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jzsuv => null()
      END TYPE

!-------------------------------------------------------------------------------
!>  Base class representing a vmec_equilibrium.
!>  @par Super Class:
!>  @ref equilibrium
!-------------------------------------------------------------------------------
      TYPE vmec_class
!>  File name of the output of v3fit.
         CHARACTER (len=data_filename_length)   :: wout_file_name
!>  File name of the v3fit namelist inout file.
         CHARACTER (len=data_filename_length)   :: vmec_file_name
!>  @ref pprofile_T describing a flux surface constant electron denisty profile.
         TYPE (pprofile), POINTER               :: ne => null()
!>  @ref pprofile_T describing a flux surface constant electron temperature
!>  profile.
         TYPE (pprofile), POINTER               :: te => null()
!>  @ref pprofile_T describing a flux surface conatsnt soft x-ray emissivity
!>  profile.
         TYPE (pprofile_pointer), DIMENSION(:), POINTER ::                     &
     &      sxrem => null()
!>  Ratio of the poloidal grid points to the poloidal grid points.
         REAL (rprec)                           :: pol_rad_ratio
!>  Instance of a @ref vmec_magnetic_cache object.
         TYPE (vmec_magnetic_cache), POINTER    ::                             &
     &      magnetic_cache => null()
!>  Cache of the vmec internal state xc array.
         REAL (rprec), DIMENSION(:), POINTER    :: xc => null()
!>  Cache of the wout file state.
         TYPE (wout_file_cache_class), POINTER  :: wout_file => null()
!>  Cache of the wout file state.
         TYPE (threed1_file_cache_class), POINTER ::                           &
     &      threed1_file => null()

!  Result file NetCDF variables
!>  External current Netcdf variable id.
         INTEGER                               :: vmec_extcur_id
!>  Toroidal current Netcdf variable id.
         INTEGER                               :: vmec_curtor_id
!>  Pressure scale Netcdf variable id.
         INTEGER                               :: vmec_pres_scale_id
!>  Boundary R cosine coeffients Netcdf variable id.
         INTEGER                               :: vmec_rbc_id
!>  Boundary R sine coeffients Netcdf variable id.
         INTEGER                               :: vmec_rbs_id
!>  Boundary Z cosine coeffients Netcdf variable id.
         INTEGER                               :: vmec_zbc_id
!>  Boundary Z sine coeffients Netcdf variable id.
         INTEGER                               :: vmec_zbs_id
!>  Internal inducance Netcdf variable id.
         INTEGER                               :: vmec_l_sub_i_id
!>  Mean elongation Netcdf variable id.
         INTEGER                               :: vmec_kappa_p_id

!  Extra reconstruction parameters
!>  Plasma phi offset.
         REAL (rprec)                          :: phi_offset = 0

!>  Index of the largest ns value.
         INTEGER                               :: ns_index
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the setting the @ref vmec_magnetic_cache types using
!>  @ref vmec_set_magnetic_cache_init or @ref vmec_set_magnetic_cache_calc. 
!-------------------------------------------------------------------------------
      INTERFACE vmec_set_magnetic_cache
         MODULE PROCEDURE vmec_set_magnetic_cache_init,                        &
     &                    vmec_set_magnetic_cache_calc
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Forward declare the interface for @ref runvmec.
!>  @note If @ref runvmec.f was implimented as a module this would not be
!>  necessary.
!-------------------------------------------------------------------------------
      INTERFACE
         SUBROUTINE runvmec(control_array, vmec_input_file,                    &
     &                      use_screen, reset_file_name)
         INTEGER, INTENT(inout), TARGET :: control_array(5)
         LOGICAL, INTENT(in)            :: use_screen
         CHARACTER(LEN=*), INTENT(in)   :: vmec_input_file
         CHARACTER(LEN=*), OPTIONAL     :: reset_file_name
         END SUBROUTINE
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref vmec_class object.
!>
!>  Allocates memory and initializes a @ref vmec_class object. Runs vmec
!>  for one iteration to initialize all the input parameters. Defined in
!>  @ref vmec_input.
!>
!>  @param[in] file_name  Filename of the vmec namelist input file.
!>  @param[in] ne         @ref pprofile_T for the electron density.
!>  @param[in] te         @ref pprofile_T for the electron temperature.
!>  @param[in] sxrem      @ref pprofile_T for the soft x-ray emissivity.
!>  @param[in] phi_offset Initial phi offset of the plasma relative to the
!>                        diagnostics in radians.
!>  @param[in] iou        Input/output unit to log messages to.
!>  @returns A pointer to a constructed @ref vmec_class object.
!-------------------------------------------------------------------------------
      FUNCTION vmec_construct(file_name, wout_file_name,                       &
     &                        ne, te, sxrem, phi_offset,                       &
     &                        pol_rad_ratio, iou)
      USE vmec_params, only: restart_flag, readin_flag, timestep_flag,         &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag
      USE vmec_input, only: l_v3fit, ns_array
      USE xstuff, only: xc

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER           :: vmec_construct
      CHARACTER (len=*), INTENT(in)        :: file_name
      CHARACTER (len=*), INTENT(in)        :: wout_file_name
      TYPE (pprofile), POINTER             :: ne
      TYPE (pprofile), POINTER             :: te
      TYPE (pprofile_pointer), DIMENSION(:), POINTER :: sxrem
      REAL (rprec), INTENT(in)             :: phi_offset
      REAL (rprec), INTENT(in)             :: pol_rad_ratio
      INTEGER, INTENT(in)                  :: iou

!  local variables
      INTEGER, DIMENSION(5)                :: ictrl_array
      INTEGER                              :: index_dat, index_end

!  Start of executable code
      WRITE (*,*) ' *** Initializing VMEC equilibrium from file ' //           &
     &            TRIM(file_name)
      WRITE (iou,*) ' *** Initializing VMEC equilibrium from file ' //         &
     &              TRIM(file_name)

      ALLOCATE(vmec_construct)
      vmec_construct%ne => ne
      vmec_construct%te => te
      vmec_construct%sxrem => sxrem

      vmec_construct%phi_offset = phi_offset
      vmec_construct%pol_rad_ratio = pol_rad_ratio

!  Initialize VMEC by taking a single step.
      ictrl_array = 0
      ictrl_array(1) = IOR(ictrl_array(1), restart_flag)
      ictrl_array(1) = IOR(ictrl_array(1), readin_flag)
      ictrl_array(1) = IOR(ictrl_array(1), timestep_flag)
      ictrl_array(3) = 1 ! Run for one iteration.
      ictrl_array(4) = 1 ! Limit VMEC to single grid.

      l_v3fit = .true.

      vmec_construct%vmec_file_name = TRIM(file_name)

      IF (wout_file_name .ne. '') THEN
         CALL runvmec(ictrl_array, file_name, .false., wout_file_name)
      ELSE
         CALL runvmec(ictrl_array, file_name, .false.)
      END IF

!  Check for errors in runvmec.
      SELECT CASE (ictrl_array(2))

         CASE (successful_term_flag, norm_term_flag, more_iter_flag)
!  Set the wout file name.
            index_dat = INDEX(file_name,'input.')
            index_end = LEN_TRIM(file_name)
            IF (index_dat .gt. 0) THEN
               vmec_construct%wout_file_name = 'wout_' //                      &
     &            TRIM(file_name(index_dat + 6:index_end)) // '.nc'
            ELSE
               vmec_construct%wout_file_name = 'wout_' //                      &
     &            TRIM(file_name(1:index_end)) // '.nc'
            END IF

         CASE DEFAULT
            CALL err_fatal('vmec_construct runvmec command failed ',           &
     &                     int = ictrl_array(2))

      END SELECT

      vmec_construct%ns_index = MAXLOC(ns_array, 1)

!  Store the inital xc array here. This is needed to reset the xc array for
!  error recovery in the initial equilibirum convergence. The vmec_save_state
!  and vmec_reset_state routines cannot be used here since those assume a wout
!  file was written.
      ALLOCATE(vmec_construct%xc(SIZE(xc)))
      vmec_construct%xc = xc

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref vmec_class object.
!>
!>  Deallocates memory and uninitializes a @ref vmec_class object.
!>
!>  @param[inout] this A @ref vmec_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_destruct(this)
      USE read_wout_mod, only: read_wout_deallocate

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER                :: this

!  local variables
      INTEGER                                   :: i

!  Start of executable code
      IF (ASSOCIATED(this%ne)) THEN
         CALL pprofile_destroy(this%ne)
         DEALLOCATE(this%ne)
         this%ne => null()
      END IF

      IF (ASSOCIATED(this%te)) THEN
         CALL pprofile_destroy(this%te)
         DEALLOCATE(this%te)
         this%te => null()
      END IF

      IF (ASSOCIATED(this%sxrem)) THEN
         DO i = 1, SIZE(this%sxrem)
            IF (ASSOCIATED(this%sxrem(i)%p)) THEN
               CALL pprofile_destroy(this%sxrem(i)%p)
               DEALLOCATE(this%sxrem(i)%p)
               this%sxrem(i)%p => null()
            END IF
         END DO
         DEALLOCATE(this%sxrem)
         this%sxrem => null()
      END IF

      IF (ASSOCIATED(this%magnetic_cache)) THEN
         DEALLOCATE(this%magnetic_cache%rsuv)
         DEALLOCATE(this%magnetic_cache%zsuv)
         DEALLOCATE(this%magnetic_cache%jrsuv)
         DEALLOCATE(this%magnetic_cache%jphisuv)
         DEALLOCATE(this%magnetic_cache%jzsuv)

         this%magnetic_cache%rsuv => null()
         this%magnetic_cache%zsuv => null()
         this%magnetic_cache%jrsuv => null()
         this%magnetic_cache%jphisuv => null()
         this%magnetic_cache%jzsuv => null()

         DEALLOCATE(this%magnetic_cache)
         this%magnetic_cache => null()
      END IF

      IF (ASSOCIATED(this%wout_file)) THEN
         CALL wout_file_cache_destruct(this%wout_file)
         this%wout_file => null()
      END IF

      IF (ASSOCIATED(this%threed1_file)) THEN
         CALL threed1_file_cache_destruct(this%threed1_file)
         this%threed1_file => null()
      END IF

      IF (ASSOCIATED(this%xc)) THEN
         DEALLOCATE(this%xc)
         this%xc => null()
      END IF

!  Deallocate VMEC allocated memory.
      CALL read_wout_deallocate
      CALL free_mem_ns(.true.)
      CALL free_mem_funct3d
      CALL free_mem_nunv

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sets the value of a reconstruction equilibrium parameter.
!>
!>  This method overrides @ref equilibrium::equilibrium_set_param. When a VMEC
!>  parameter is changed, propagate the changes to the VMEC internal state and
!>  inform the caller that the equilibrium needs reconvergence.
!>
!>  @param[inout] this    A @ref vmec_class instance.
!>  @param[in]    id      ID of the parameter.
!>  @param[in]    i_index The ith index of the parameter.
!>  @param[in]    j_index The jth index of the parameter.
!>  @param[in]    value   The value of the parameter.
!>  @returns True if VMEC is left in an unconverged state and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION vmec_set_param(this, id, i_index, j_index, value)
      USE vmec_input, only: rbc, zbs, rbs, zbc,                                &
     &                      ac, ac_aux_s, ac_aux_f,                            &
     &                      ai, ai_aux_s, ai_aux_f,                            &
     &                      am, am_aux_s, am_aux_f,                            &
     &                      extcur, curtor, phiedge, pres_scale,               &
     &                      bloat, mgrid_file, nfp, lfreeb
      USE xstuff, only: xc, xstore, xcdot
      USE vmec_main, only: irzloff, currv, ivac, fsqr, fsqz
      USE stel_constants, only: mu0
      USE mgrid_mod, only: read_mgrid, mgrid_path_old
      USE vacmod0, only: nv

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                          :: vmec_set_param
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(in)              :: id
      INTEGER, INTENT(in)              :: i_index
      INTEGER, INTENT(in)              :: j_index
      REAL (rprec), INTENT(in)         :: value

!  local variables
      REAL (rprec)                     :: phiedge_ratio
      INTEGER                          :: status

!  Start of executable code

      SELECT CASE (id)

         CASE (vmec_rbc_id)
            vmec_set_param = .true.
            rbc(i_index, j_index) = value
            CALL profil1d(xc, xcdot, .false.)
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               fsqr = 1.0
               fsqz = 1.0
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
               CALL profil3d(xstore(1), xstore(1 + irzloff), .true.,           &
     &                       .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_zbs_id)
            vmec_set_param = .true.
            zbs(i_index, j_index) = value
            CALL profil1d(xc, xcdot, .false.)
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               fsqr = 1.0
               fsqz = 1.0
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
               CALL profil3d(xstore(1), xstore(1 + irzloff), .true.,           &
     &                       .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_rbs_id)
            vmec_set_param = .true.
            rbs(i_index, j_index) = value
            CALL profil1d(xc, xcdot, .false.)
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               fsqr = 1.0
               fsqz = 1.0
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
               CALL profil3d(xstore(1), xstore(1 + irzloff), .true.,           &
     &                       .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_zbc_id)
            vmec_set_param = .true.
            zbc(i_index, j_index) = value
            CALL profil1d(xc, xcdot, .false.)
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               fsqr = 1.0
               fsqz = 1.0
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
               CALL profil3d(xstore(1), xstore(1 + irzloff), .true.,           &
     &                       .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_ac_id)
            vmec_set_param = .true.
            ac(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_ac_aux_s_id)
            vmec_set_param = .true.
            ac_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_ac_aux_f_id)
            vmec_set_param = .true.
            ac_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_ai_id)
            vmec_set_param = .true.
            ai(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_ai_aux_s_id)
            vmec_set_param = .true.
            ai_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_ai_aux_f_id)
            vmec_set_param = .true.
            ai_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_am_id)
            vmec_set_param = .true.
            am(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_am_aux_s_id)
            vmec_set_param = .true.
            am_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_am_aux_f_id)
            vmec_set_param = .true.
            am_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_pp_ne_b_id)
            vmec_set_param = .false.
            this%ne%b(i_index) = value

         CASE (vmec_pp_ne_as_id)
            vmec_set_param = .false.
            this%ne%as(i_index) = value

         CASE (vmec_pp_ne_af_id)
            vmec_set_param = .false.
            this%ne%af(i_index) = value

         CASE (vmec_pp_te_b_id)
            vmec_set_param = .false.
            this%te%b(i_index) = value

         CASE (vmec_pp_te_as_id)
            vmec_set_param = .false.
            this%te%as(i_index) = value

         CASE (vmec_pp_te_af_id)
            vmec_set_param = .false.
            this%te%af(i_index) = value

         CASE (vmec_pp_sxrem_b_id)
            vmec_set_param = .false.
            this%sxrem(i_index)%p%b(j_index) = value

         CASE (vmec_pp_sxrem_as_id)
            vmec_set_param = .false.
            this%sxrem(i_index)%p%as(j_index) = value

         CASE (vmec_pp_sxrem_af_id)
            vmec_set_param = .false.
            this%sxrem(i_index)%p%af(j_index) = value

         CASE (vmec_extcur_id)
            vmec_set_param = .true.
            extcur(i_index) = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            mgrid_path_old = ''
            CALL read_mgrid(mgrid_file, extcur, nv, nfp, .false.,              &
     &                      status)
            CALL assert_eq(0, status, 'vmec_set_param: failed to ' //          &
     &                     'mgrid_file')
            ivac = 1
            fsqr = 1.0
            fsqz = 1.0

         CASE (vmec_curtor_id)
            vmec_set_param = .true.
            curtor = value
            currv = mu0*curtor
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_phiedge_id)
            vmec_set_param = .true.
            IF (phiedge .ne. 0.0) THEN
               phiedge_ratio = value/phiedge
            ELSE
               phiedge_ratio = 1.0
            END IF
            phiedge = value
            xc(1 + 2*irzloff:3*irzloff) =                                      &
     &         phiedge_ratio*xc(1 + 2*irzloff:3*irzloff)
            xstore(1 + 2*irzloff:3*irzloff) =                                  &
     &         phiedge_ratio*xstore(1 + 2*irzloff:3*irzloff)
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_pres_scale_id)
            vmec_set_param = .true.
            pres_scale = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_bloat_id)
            vmec_set_param = .true.
            bloat = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

!  For the phi_offset, the magnetic cache needs to be updated. Since the phi
!  offset doesn't alter the equilibrium, the magnetics were not getting updated
!  since the equilibrium was not reconverged. Call vmec_set_magnetic_cache_calc
!  to update the magnetic cache without reconverging VMEC.
         CASE (vmec_phi_offset_id)
            vmec_set_param = .false.
            this%phi_offset = value
            IF (ASSOCIATED(this%magnetic_cache)) THEN
               CALL vmec_set_magnetic_cache(this)
            END IF

         CASE DEFAULT
!>  @todo FIXME: Consider throwing a fatal error here instead of silently doing
!>  nothing.
            vmec_set_param = .true.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Set magnetic cache initialy.
!>
!>  This method overrides @ref equilibrium::equilibrium_set_magnetic_cache. This
!>  allocates a @ref vmec_magnetic_cache structure.
!>
!>  @param[inout] this             A @ref vmec_class instance.
!>  @param[in]    total_tor_planes Total number of toroidal planes.
!>  @param[in]    tor_planes_used  Number of stellarator symmetric toroidal
!>                                 planes.
!>  @param[in]    nfp              Number of field periods the diagnostics were
!>                                 generated with.
!>  @note Stellarator symmetry may be exploited speed up computation and reduce
!>  memory usage. If stellarator symmetry is not used then tor_planes_used
!>  should match total_tor_planes.
!>
!>  @todo This method explicity uses only the first element of the ns_array.
!>  In effect, this excludes support for multigrid VMEC. To support multi grid,
!>  the vmec_set_magnetic_cache_init will need to be reset every time the grid
!>  size changes. Currently the grid are allocated once and never changed.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_magnetic_cache_init(this, total_tor_planes,          &
     &                                        tor_planes_used, nfp)
      USE vmec_input, only: ns_array
      USE stel_constants, only: twopi

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(in)              :: total_tor_planes
      INTEGER, INTENT(in)              :: tor_planes_used
      INTEGER, INTENT(in)              :: nfp

!  local variables
      INTEGER                          :: numU
      INTEGER                          :: ns

!  Start of executable code.
      ALLOCATE(this%magnetic_cache)

      ns = ns_array(this%ns_index)
      numU = MAX(INT(this%pol_rad_ratio*ns), 1)

!  Set the differental volume elements.
      this%magnetic_cache%ds = 1.0/(ns - 1.0)
      this%magnetic_cache%du = twopi/numU
      this%magnetic_cache%dv = (twopi/nfp)/total_tor_planes

!  When stellarator symmetry is used, do not need to store all the toroidal
!  planes in a single field period.
      ALLOCATE(this%magnetic_cache%rsuv(ns, numU, tor_planes_used))
      ALLOCATE(this%magnetic_cache%zsuv(ns, numU, tor_planes_used))
      ALLOCATE(this%magnetic_cache%jrsuv(ns, numU, tor_planes_used))
      ALLOCATE(this%magnetic_cache%jphisuv(ns, numU, tor_planes_used))
      ALLOCATE(this%magnetic_cache%jzsuv(ns, numU, tor_planes_used))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set magnetic cache initialy.
!>
!>  After the equilibrium has been converged calculate the r, z, jr, jphi and
!>  jz on grid of s, u, v points.
!>
!>  @param[inout] this A @ref vmec_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_magnetic_cache_calc(this)
      USE read_wout_mod, only: lasym, xm, xn, mnmax, xm_nyq, xn_nyq,           &
     &                         mnmax_nyq, rmnc, zmns, rmns, zmnc,              &
     &                         currumnc, currvmnc, currumns, currvmns

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this

!  local variables
      INTEGER                                   :: s, u, v
      INTEGER                                   :: numS, numU, numV
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosz, sinz
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: cosu, cosv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: sinu, sinv
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosz_nyq, sinz_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: cosu_nyq, cosv_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: sinu_nyq, sinv_nyq
      REAL (rprec)                              :: ru, rv, zu, zv
      REAL (rprec)                              :: ju, jv
      REAL (rprec)                              :: phi

!  Start of executable code.
!  Get the size of the SUV grid.
      numS = SIZE(this%magnetic_cache%rsuv, 1)
      numU = SIZE(this%magnetic_cache%rsuv, 2)
      numV = SIZE(this%magnetic_cache%rsuv, 3)

      ALLOCATE(cosz(mnmax), sinz(mnmax))
      ALLOCATE(cosu(mnmax, numU), sinu(mnmax, numU),                           &
     &         cosv(mnmax, numV), sinv(mnmax, numV))

      ALLOCATE(cosz_nyq(mnmax_nyq), sinz_nyq(mnmax_nyq))
      ALLOCATE(cosu_nyq(mnmax_nyq, numU), sinu_nyq(mnmax_nyq, numU),           &
     &         cosv_nyq(mnmax_nyq, numV), sinv_nyq(mnmax_nyq, numV))

!  U grid is square with the s grid.
!>  @todo FIXME: There is no reason to recalulate the sin and cosin for the u
!>  componentarrays every time. The v components can change if phi_offset is
!>  being reconstructed. The u components can be cached in the
!>  @ref vmec_magnetic_cache and calculated once in
!>  @ref vmec_set_magnetic_cache_init.
      DO u = 1, numU
         cosu(:,u) = COS(xm*(u - 1)*this%magnetic_cache%du)
         sinu(:,u) = SIN(xm*(u - 1)*this%magnetic_cache%du)
         cosu_nyq(:,u) = COS(xm_nyq*(u - 1)*this%magnetic_cache%du)
         sinu_nyq(:,u) = SIN(xm_nyq*(u - 1)*this%magnetic_cache%du)
      END DO
      DO v = 1, numV
         phi = (v - 1)*this%magnetic_cache%dv + this%phi_offset
         cosv(:,v) = COS(xn*phi)
         sinv(:,v) = SIN(xn*phi)
         cosv_nyq(:,v) = COS(xn_nyq*phi)
         sinv_nyq(:,v) = SIN(xn_nyq*phi)
      END DO

!  Set the suv grids
!
!  ru = dR/du
!  rv = dR/dv   Note: dRHat/dv = PhiHat
!
!  zu = dZ/du
!  zu = dZ/dv
!
!  jr = ju*ru + jv*rv
!  jphi = jv*R
!  jz = ju*zu + jv*zv
      DO u = 1, numU

         DO v = 1, numV

            cosz = cosu(:,u)*cosv(:,v) + sinu(:,u)*sinv(:,v)
            sinz = sinu(:,u)*cosv(:,v) - cosu(:,u)*sinv(:,v)

            cosz_nyq = cosu_nyq(:,u)*cosv_nyq(:,v)                             &
     &               + sinu_nyq(:,u)*sinv_nyq(:,v)
            sinz_nyq = sinu_nyq(:,u)*cosv_nyq(:,v)                             &
     &               - cosu_nyq(:,u)*sinv_nyq(:,v)

            DO s = 1, numS

               this%magnetic_cache%rsuv(s,u,v) = SUM(rmnc(:,s)*cosz)
               this%magnetic_cache%zsuv(s,u,v) = SUM(zmns(:,s)*sinz)
               ru = -SUM(xm*rmnc(:,s)*sinz)
               zu = SUM(xm*zmns(:,s)*cosz)
               rv = SUM(xn*rmnc(:,s)*sinz)
               zv = -SUM(xn*zmns(:,s)*cosz)

               ju = SUM(currumnc(:,s)*cosz_nyq)
               jv = SUM(currvmnc(:,s)*cosz_nyq)

               IF (lasym) THEN
                  this%magnetic_cache%rsuv(s,u,v) =                            &
     &               this%magnetic_cache%rsuv(s,u,v) +                         &
     &               SUM(rmns(:,s)*sinz)
                  this%magnetic_cache%zsuv(s,u,v) =                            &
     &               this%magnetic_cache%zsuv(s,u,v) +                         &
     &               SUM(zmnc(:,s)*cosz)
                  ru = ru + SUM(xm*rmns(:,s)*cosz)
                  zu = zu - SUM(xm*zmnc(:,s)*sinz)
                  rv = rv - SUM(xn*rmns(:,s)*cosz)
                  zv = zv + SUM(xn*zmnc(:,s)*sinz)

                  ju = ju + SUM(currumns(:,s)*sinz_nyq)
                  jv = jv + SUM(currvmns(:,s)*sinz_nyq)
               END IF

               this%magnetic_cache%jrsuv(s,u,v) = ju*ru + jv*rv
               this%magnetic_cache%jphisuv(s,u,v) =                            &
     &            jv*this%magnetic_cache%rsuv(s,u,v)
               this%magnetic_cache%jzsuv(s,u,v) = ju*zu + jv*zv
            END DO
         END DO
      END DO

      DEALLOCATE(cosz, sinz)
      DEALLOCATE(cosu, sinu, cosv, sinv)
      DEALLOCATE(cosz_nyq, sinz_nyq)
      DEALLOCATE(cosu_nyq, sinu_nyq, cosv_nyq, sinv_nyq)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Get the id for a reconstruction parameter.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_param_id.
!>
!>  @param[in] this       A @ref vmec_class instance.
!>  @param[in] param_name Name of a reconstruction parameter.
!>  @returns The id for a reconstruction parameter.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_param_id(this, param_name)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: vmec_get_param_id
      TYPE (vmec_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in) :: param_name

!  Start of executable code
      SELECT CASE (TRIM(param_name))

         CASE ('ac')
            vmec_get_param_id = vmec_ac_id

         CASE ('ac_aux_s')
            vmec_get_param_id = vmec_ac_aux_s_id

         CASE ('ac_aux_f')
            vmec_get_param_id = vmec_ac_aux_f_id

         CASE ('ai')
            vmec_get_param_id = vmec_ai_id

         CASE ('ai_aux_s')
            vmec_get_param_id = vmec_ai_aux_s_id

         CASE ('ai_aux_f')
            vmec_get_param_id = vmec_ai_aux_f_id

         CASE ('am')
            vmec_get_param_id = vmec_am_id

         CASE ('am_aux_s')
            vmec_get_param_id = vmec_am_aux_s_id

         CASE ('am_aux_f')
            vmec_get_param_id = vmec_am_aux_f_id

         CASE ('bloat')
            vmec_get_param_id = vmec_bloat_id

         CASE ('rbc')
            vmec_get_param_id = vmec_rbc_id

         CASE ('zbs')
            vmec_get_param_id = vmec_zbs_id

         CASE ('rbs')
            vmec_get_param_id = vmec_rbs_id

         CASE ('zbc')
            vmec_get_param_id = vmec_zbc_id

         CASE ('extcur')
            vmec_get_param_id = vmec_extcur_id

         CASE ('curtor')
            vmec_get_param_id = vmec_curtor_id

         CASE ('phiedge')
            vmec_get_param_id = vmec_phiedge_id

         CASE ('pres_scale')
            vmec_get_param_id = vmec_pres_scale_id

         CASE ('rmnc')
            vmec_get_param_id = vmec_rmnc_id

         CASE ('zmns')
            vmec_get_param_id = vmec_zmns_id

         CASE ('lmns')
            vmec_get_param_id = vmec_lmns_id

         CASE ('gmnc')
            vmec_get_param_id = vmec_gmnc_id

         CASE ('bsubumnc')
            vmec_get_param_id = vmec_bsubumnc_id

         CASE ('bsubvmnc')
            vmec_get_param_id = vmec_bsubvmnc_id

         CASE ('bsubsmns')
            vmec_get_param_id = vmec_bsubsmns_id

         CASE ('bsupumnc')
            vmec_get_param_id = vmec_bsupumnc_id

         CASE ('bsupvmnc')
            vmec_get_param_id = vmec_bsupvmnc_id

         CASE ('rmns')
            vmec_get_param_id = vmec_rmns_id

         CASE ('zmnc')
            vmec_get_param_id = vmec_zmnc_id

         CASE ('lmnc')
            vmec_get_param_id = vmec_lmnc_id

         CASE ('gmns')
            vmec_get_param_id = vmec_gmns_id

         CASE ('bsubumns')
            vmec_get_param_id = vmec_bsubumns_id

         CASE ('bsubvmns')
            vmec_get_param_id = vmec_bsubvmns_id

         CASE ('bsubsmnc')
            vmec_get_param_id = vmec_bsubsmnc_id

         CASE ('bsupumns')
            vmec_get_param_id = vmec_bsupumns_id

         CASE ('bsupvmns')
            vmec_get_param_id = vmec_bsupvmns_id

         CASE ('phi')
            vmec_get_param_id = vmec_phi_id

         CASE ('iotaf')
            vmec_get_param_id = vmec_iotaf_id

         CASE ('iotas')
            vmec_get_param_id = vmec_iotas_id

         CASE ('pp_ne_b')
            vmec_get_param_id = vmec_pp_ne_b_id

         CASE ('pp_ne_as')
            vmec_get_param_id = vmec_pp_ne_as_id

         CASE ('pp_ne_af')
            vmec_get_param_id = vmec_pp_ne_af_id

         CASE ('pp_sxrem_b_a')
            vmec_get_param_id = vmec_pp_sxrem_b_id

         CASE ('pp_sxrem_as_a')
            vmec_get_param_id = vmec_pp_sxrem_as_id

         CASE ('pp_sxrem_af_a')
            vmec_get_param_id = vmec_pp_sxrem_af_id

         CASE ('pp_te_b')
            vmec_get_param_id = vmec_pp_te_b_id

         CASE ('pp_te_as')
            vmec_get_param_id = vmec_pp_te_as_id

         CASE ('pp_te_af')
            vmec_get_param_id = vmec_pp_te_af_id

         CASE ('vvc_smaleli')
            vmec_get_param_id = vmec_vvc_smaleli_id

         CASE ('vvc_kappa_p')
            vmec_get_param_id = vmec_vvc_kappa_p_id

         CASE ('betatot')
            vmec_get_param_id = vmec_betatot_id

         CASE ('betapol')
            vmec_get_param_id = vmec_betapol_id

         CASE ('betator')
            vmec_get_param_id = vmec_betator_id

         CASE ('betaxis')
            vmec_get_param_id = vmec_betaxis_id

         CASE ('phi_offset')
            vmec_get_param_id = vmec_phi_offset_id

         CASE DEFAULT
            vmec_get_param_id = data_no_id

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the value of a reconstruction VMEC parameter.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_param_value.
!>
!>  @param[in] this    A @ref vmec_class instance.
!>  @param[in] id      ID of the parameter.
!>  @param[in] i_index The ith index of the parameter.
!>  @param[in] j_index The jth index of the parameter.
!>  @returns The value of the parameter.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_param_value(this, id, i_index, j_index)
      USE vmec_input, only: rbc, zbs, rbs, zbc,                                &
     &                      ac, ac_aux_s, ac_aux_f,                            &
     &                      ai, ai_aux_s, ai_aux_f,                            &
     &                      am, am_aux_s, am_aux_f,                            &
     &                      extcur, curtor, phiedge, pres_scale,               &
     &                      bloat
      USE read_wout_mod, only: rmnc, zmns, rmns, zmnc, lmns, lmnc,             &
     &                         gmnc, gmns, bsubumnc, bsubvmnc,                 &
     &                         bsubsmns, bsupumnc, bsupvmnc,                   &
     &                         bsubumns, bsubvmns,                             &
     &                         bsubsmnc, bsupumns, bsupvmns,                   &
     &                         phi, iotaf, iotas,                              &
     &                         betatot, betapol, betator, betaxis
      USE v3f_vmec_comm

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: vmec_get_param_value
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id
      INTEGER, INTENT(in)           :: i_index
      INTEGER, INTENT(in)           :: j_index

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_ac_id)
            vmec_get_param_value = ac(i_index)

         CASE (vmec_ac_aux_s_id)
            vmec_get_param_value = ac_aux_s(i_index)

         CASE (vmec_ac_aux_f_id)
            vmec_get_param_value = ac_aux_f(i_index)

         CASE (vmec_ai_id)
            vmec_get_param_value = ai(i_index)

         CASE (vmec_ai_aux_s_id)
            vmec_get_param_value = ai_aux_s(i_index)

         CASE (vmec_ai_aux_f_id)
            vmec_get_param_value = ai_aux_f(i_index)

         CASE (vmec_am_id)
            vmec_get_param_value = am(i_index)

         CASE (vmec_am_aux_s_id)
            vmec_get_param_value = am_aux_s(i_index)

         CASE (vmec_am_aux_f_id)
            vmec_get_param_value = am_aux_f(i_index)

         CASE (vmec_bloat_id)
            vmec_get_param_value = bloat

         CASE (vmec_rbc_id)
            vmec_get_param_value = rbc(i_index, j_index)

         CASE (vmec_zbs_id)
            vmec_get_param_value = zbs(i_index, j_index)

         CASE (vmec_rbs_id)
            vmec_get_param_value = rbc(i_index, j_index)

         CASE (vmec_zbc_id)
            vmec_get_param_value = zbc(i_index, j_index)

         CASE (vmec_extcur_id)
            vmec_get_param_value = extcur(i_index)

         CASE (vmec_curtor_id)
            vmec_get_param_value = curtor

         CASE (vmec_phiedge_id)
            vmec_get_param_value = phiedge

         CASE (vmec_pres_scale_id)
            vmec_get_param_value = pres_scale

         CASE (vmec_rmnc_id)
            vmec_get_param_value = rmnc(i_index, j_index)

         CASE (vmec_zmns_id)
            vmec_get_param_value = zmns(i_index, j_index)

         CASE (vmec_lmns_id)
            vmec_get_param_value = lmns(i_index, j_index)

         CASE (vmec_gmnc_id)
            vmec_get_param_value = gmnc(i_index, j_index)

         CASE (vmec_bsubumnc_id)
            vmec_get_param_value = bsubumnc(i_index, j_index)

         CASE (vmec_bsubvmnc_id)
            vmec_get_param_value = bsubvmnc(i_index, j_index)

         CASE (vmec_bsubsmns_id)
            vmec_get_param_value = bsubsmns(i_index, j_index)

         CASE (vmec_bsupumnc_id)
            vmec_get_param_value = bsupumnc(i_index, j_index)

         CASE (vmec_bsupvmnc_id)
            vmec_get_param_value = bsupvmnc(i_index, j_index)

         CASE (vmec_rmns_id)
            vmec_get_param_value = rmns(i_index, j_index)

         CASE (vmec_zmnc_id)
            vmec_get_param_value = zmnc(i_index, j_index)

         CASE (vmec_lmnc_id)
            vmec_get_param_value = lmnc(i_index, j_index)

         CASE (vmec_gmns_id)
            vmec_get_param_value = gmns(i_index, j_index)

         CASE (vmec_bsubumns_id)
            vmec_get_param_value = bsubumns(i_index, j_index)

         CASE (vmec_bsubvmns_id)
            vmec_get_param_value = bsubvmns(i_index, j_index)

         CASE (vmec_bsubsmnc_id)
            vmec_get_param_value = bsubsmnc(i_index, j_index)

         CASE (vmec_bsupumns_id)
            vmec_get_param_value = bsupumns(i_index, j_index)

         CASE (vmec_bsupvmns_id)
            vmec_get_param_value = bsupvmns(i_index, j_index)

         CASE (vmec_phi_id)
            vmec_get_param_value = phi(i_index)

         CASE (vmec_iotaf_id)
            vmec_get_param_value = iotaf(i_index)

         CASE (vmec_iotas_id)
            vmec_get_param_value = iotas(i_index)

         CASE (vmec_pp_ne_b_id)
            vmec_get_param_value = this%ne%b(i_index)

         CASE (vmec_pp_ne_as_id)
            vmec_get_param_value = this%ne%as(i_index)

         CASE (vmec_pp_ne_af_id)
            vmec_get_param_value = this%ne%af(i_index)

         CASE (vmec_pp_sxrem_b_id)
            vmec_get_param_value = this%sxrem(i_index)%p%b(j_index)

         CASE (vmec_pp_sxrem_as_id)
            vmec_get_param_value = this%sxrem(i_index)%p%as(j_index)

         CASE (vmec_pp_sxrem_af_id)
            vmec_get_param_value = this%sxrem(i_index)%p%af(j_index)

         CASE (vmec_pp_te_b_id)
            vmec_get_param_value = this%te%b(i_index)

         CASE (vmec_pp_te_as_id)
            vmec_get_param_value = this%te%as(i_index)

         CASE (vmec_pp_te_af_id)
            vmec_get_param_value = this%te%af(i_index)

         CASE (vmec_vvc_smaleli_id)
            vmec_get_param_value = vvc_smaleli

         CASE (vmec_vvc_kappa_p_id)
            vmec_get_param_value = vvc_kappa_p

         CASE (vmec_betatot_id)
            vmec_get_param_value = betatot

         CASE (vmec_betapol_id)
            vmec_get_param_value = betapol

         CASE (vmec_betator_id)
            vmec_get_param_value = betator

         CASE (vmec_betaxis_id)
            vmec_get_param_value = betaxis

         CASE (vmec_phi_offset_id)
            vmec_get_param_value = this%phi_offset

         CASE DEFAULT
            vmec_get_param_value = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the name of a reconstruction VMEC parameter.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_param_name.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns The name of the parameter.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_param_name(this, id)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=data_name_length) :: vmec_get_param_name
      TYPE (vmec_class), INTENT(in)   :: this
      INTEGER, INTENT(in)             :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_ac_id)
            vmec_get_param_name = 'ac'

         CASE (vmec_ac_aux_s_id)
            vmec_get_param_name = 'ac_aux_s'

         CASE (vmec_ac_aux_f_id)
            vmec_get_param_name = 'ac_aux_f'

         CASE (vmec_ai_id)
            vmec_get_param_name = 'ai'

         CASE (vmec_ai_aux_s_id)
            vmec_get_param_name = 'ai_aux_s'

         CASE (vmec_ai_aux_f_id)
            vmec_get_param_name = 'ai_aux_f'

         CASE (vmec_am_id)
            vmec_get_param_name = 'am'

         CASE (vmec_am_aux_s_id)
            vmec_get_param_name = 'am_aux_s'

         CASE (vmec_am_aux_f_id)
            vmec_get_param_name = 'am_aux_f'

         CASE (vmec_bloat_id)
            vmec_get_param_name = 'bloat'

         CASE (vmec_rbc_id)
            vmec_get_param_name = 'rbc'

         CASE (vmec_zbs_id)
            vmec_get_param_name = 'zbs'

         CASE (vmec_rbs_id)
            vmec_get_param_name = 'rbs'

         CASE (vmec_zbc_id)
            vmec_get_param_name = 'zbc'

         CASE (vmec_extcur_id)
            vmec_get_param_name = 'extcur'

         CASE (vmec_curtor_id)
            vmec_get_param_name = 'curtor'

         CASE (vmec_phiedge_id)
            vmec_get_param_name = 'phiedge'

         CASE (vmec_pres_scale_id)
            vmec_get_param_name = 'pres_scale'

         CASE (vmec_rmnc_id)
            vmec_get_param_name = 'rmnc'

         CASE (vmec_zmns_id)
            vmec_get_param_name = 'zmns'

         CASE (vmec_lmns_id)
            vmec_get_param_name = 'lmns'

         CASE (vmec_gmnc_id)
            vmec_get_param_name = 'gmnc'

         CASE (vmec_bsubumnc_id)
            vmec_get_param_name = 'bsubumnc'

         CASE (vmec_bsubvmnc_id)
            vmec_get_param_name = 'bsubvmnc'

         CASE (vmec_bsubsmns_id)
            vmec_get_param_name = 'bsubsmns'

         CASE (vmec_bsupumnc_id)
            vmec_get_param_name = 'bsupumnc'
 
         CASE (vmec_bsupvmnc_id)
            vmec_get_param_name = 'bsupvmnc'

         CASE (vmec_rmns_id)
            vmec_get_param_name = 'rmns'

         CASE (vmec_zmnc_id)
            vmec_get_param_name = 'zmnc'

         CASE (vmec_lmnc_id)
            vmec_get_param_name = 'lmnc'

         CASE (vmec_gmns_id)
            vmec_get_param_name = 'gmns'
 
         CASE (vmec_bsubumns_id)
            vmec_get_param_name = 'bsubumns'

         CASE (vmec_bsubvmns_id)
            vmec_get_param_name = 'bsubvmns'

         CASE (vmec_bsubsmnc_id)
            vmec_get_param_name = 'bsubsmnc'

         CASE (vmec_bsupumns_id)
            vmec_get_param_name = 'bsupumns'

         CASE (vmec_bsupvmns_id)
            vmec_get_param_name = 'bsupvmns'

         CASE (vmec_phi_id)
            vmec_get_param_name = 'phi'

         CASE (vmec_iotaf_id)
            vmec_get_param_name = 'iotaf'

         CASE (vmec_iotas_id)
            vmec_get_param_name = 'iotas'

         CASE (vmec_pp_ne_b_id)
            vmec_get_param_name = 'pp_ne_b'

         CASE (vmec_pp_ne_as_id)
            vmec_get_param_name = 'pp_ne_as'

         CASE (vmec_pp_ne_af_id)
            vmec_get_param_name = 'pp_ne_af'

         CASE (vmec_pp_sxrem_b_id)
            vmec_get_param_name = 'pp_sxrem_b_a'

         CASE (vmec_pp_sxrem_as_id)
            vmec_get_param_name = 'pp_sxrem_as_a'

         CASE (vmec_pp_sxrem_af_id)
            vmec_get_param_name = 'pp_sxrem_af_a'

         CASE (vmec_pp_te_b_id)
            vmec_get_param_name = 'pp_te_b'

         CASE (vmec_pp_te_as_id)
            vmec_get_param_name = 'pp_te_as'

         CASE (vmec_pp_te_af_id)
            vmec_get_param_name = 'pp_te_af'

         CASE (vmec_vvc_smaleli_id)
            vmec_get_param_name = 'vvc_smaleli'

         CASE (vmec_vvc_kappa_p_id)
            vmec_get_param_name = 'vvc_kappa_p'

         CASE (vmec_betatot_id)
            vmec_get_param_name = 'betatot'

         CASE (vmec_betapol_id)
            vmec_get_param_name = 'betapol'

         CASE (vmec_betator_id)
            vmec_get_param_name = 'betator'

         CASE (vmec_betaxis_id)
            vmec_get_param_name = 'betaxis'

         CASE (vmec_phi_offset_id)
            vmec_get_param_name = 'phi_offset'

         CASE DEFAULT
            STOP 'Invalid VMEC parameter id.'

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron density at a position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_ne. If no density
!>  profile was created, return zero.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the electron density at.
!>  @returns The electron density at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_ne(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_ne
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      IF (ASSOCIATED(this%ne)) THEN
         vmec_get_ne = pprofile_evaluate(this%ne,                              &
     &                                   vmec_get_s(this, x_cart))
      ELSE
         vmec_get_ne = 0.0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron temperature at a position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_te. If no
!>  temperature profile was created, return zero.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the electron temperature at.
!>  @returns The electron temperature at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_te(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_te
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      IF (ASSOCIATED(this%te)) THEN
         vmec_get_te = pprofile_evaluate(this%te,                              &
     &                                   vmec_get_s(this, x_cart))
      ELSE
         vmec_get_te = 0.0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the soft x-ray emissivity at a position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_sxrem. If no soft
!>  x-ray emissivity profile was created, return zero.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the soft x-ray emissivity at.
!>  @param[in] index  Index of the soft x-ray emissivity profile.
!>  @returns The soft x-ray emissivity at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_sxrem(this, x_cart, index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_sxrem
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      INTEGER, INTENT(in)                    :: index

!  Start of executable code
      IF (ASSOCIATED(this%sxrem)) THEN
         vmec_get_sxrem = pprofile_evaluate(this%sxrem(index)%p,               &
     &                                      vmec_get_s(this, x_cart))
      ELSE
         vmec_get_sxrem = 0.0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the plasma pressure at a position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_p.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the plasma pressure at.
!>  @returns The plasma pressure at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_p(this, x_cart)
      USE stel_constants, only: mu0

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_p
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Forward declare the interface for pmass. This explicitly defined the type for
!  pmass to avoid a compiler error. If profile_functions.f was implimented as a
!  module this would not be needed.
      INTERFACE
         FUNCTION pmass(s)
         USE stel_kinds
         REAL (rprec)             :: pmass
         REAL (rprec), INTENT(in) :: s
         END FUNCTION
      END INTERFACE

!  Start of executable code
      vmec_get_p = pmass(vmec_get_s(this, x_cart))/mu0

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the magnetic field vector at a position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_B_vec.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the magnetic field vector at.
!>  @returns The magnetic field vector at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_B_vec(this, x_cart)
      USE coordinate_utilities, ONLY : cyl_to_cart_vec, cart_to_cyl
      USE vmec_utils, ONLY : GetBcyl_WOUT

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: vmec_get_B_vec
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  local variables
      REAL (rprec), DIMENSION(3)             :: r_cyl

!  Start of executable code
      r_cyl = cart_to_cyl(x_cart)
      vmec_get_B_vec = 0.0
      r_cyl(2) = r_cyl(2) + this%phi_offset

      CALL GetBcyl_WOUT(r_cyl(1), r_cyl(2), r_cyl(3),                          &
     &                  vmec_get_B_vec(1), vmec_get_B_vec(2),                  &
     &                  vmec_get_B_vec(3))

      vmec_get_B_vec = cyl_to_cart_vec(r_cyl, vmec_get_B_vec)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the loop integrated magnetic field at a position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_Int_B_dphi. This
!>  computes Int[B*dl]
!>
!>  @param[in] this  A @ref vmec_class instance.
!>  @param[in] r     S position to integrate about.
!>  @param[in] theta U angle to integrate about.
!>  @returns The loop integrated magnetic field at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_Int_B_dphi(this, r, theta)
      USE line_segment, only: line_seg
      USE stel_constants, only: twopi
      USE read_wout_mod, only: bsubvmnc, ns

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: vmec_get_Int_B_dphi
      TYPE (vmec_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)      :: r
      REAL (rprec), INTENT(in)      :: theta

!  local variables
      REAL (rprec)                  :: bsubv00c, ds
      INTEGER                       :: i
      REAL (rprec), DIMENSION(ns)   :: s

!  Start of executable code
      ds = 1.0/(ns - 1)
      DO i = 0, ns - 1
         s(i + 1) = i*ds
      END DO

      CALL line_seg(r, bsubv00c, s, bsubvmnc(1,:), ns)

      vmec_get_Int_B_dphi = twopi*bsubv00c

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets VMEC s value at position.
!>
!>  Converts from cartesian coordinates to VMEC's flux coordinates.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the s position at.
!>  @returns The s position at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_s(this, x_cart)
      USE cyl_flux, only: cyl2flx
      USE coordinate_utilities, only: cart_to_cyl
      USE read_wout_mod, only: rzl_local, ns, ntmax, lthreed
      USE vmec_input, only: ntor, mpol, nfp, lasym

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_s
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  local variables
      REAL (rprec), DIMENSION(3)             :: r_cyl
      REAL (rprec), DIMENSION(3)             :: s_flux
      REAL (rprec)                           :: fmin
      INTEGER                                :: info
      INTEGER                                :: nfe

!  Start of executable code
      s_flux = 0.0
      fmin = 0.0
      info = 0
      nfe = 0

      r_cyl = cart_to_cyl(x_cart)
      r_cyl(2) = nfp*(r_cyl(2) + this%phi_offset)

      CALL cyl2flx(rzl_local, r_cyl, s_flux, ns, ntor, mpol, ntmax,            &
     &             lthreed, lasym, info, nfe, fmin)

      vmec_get_s = s_flux(1)

      IF (vmec_get_s .le. 1.0) THEN
         IF (fmin .gt.  1.0E-8_rprec) THEN
            WRITE(*,*) "Warn: VMEC could not converge to flux position."       &
            WRITE(*,*) info, fmin, vmec_get_s
         END IF
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the r and z positions of the outer surface at a toroidal angle.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_plasma_edge.
!>
!>  @param[in]    this A @ref vmec_class instance.
!>  @param[in]    phi  Toroidal angle to determine the outer surface at.
!>  @param[inout] r    The radial postions of the other surface in a single
!>                     toroidal angle.
!>  @param[inout] z    The Z postions of the other surface in a single toroidal
!>                     angle.
!>  @returns The number of elements in the r and z arrays.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_plasma_edge(this, phi, r, z)
      USE vmec_input, only: mpol
      USE read_wout_mod, only: xm, xn, rmnc, zmns, rmns, zmnc, lasym, ns
      USE stel_constants, only: twopi

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: vmec_get_plasma_edge
      TYPE (vmec_class), INTENT(in)                      :: this
      REAL (rprec), INTENT (in)                          :: phi
      REAL (rprec), DIMENSION(:), POINTER                :: r
      REAL (rprec), DIMENSION(:), POINTER                :: z

!  local variables
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: cosv, sinv
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: cosu, sinu
      INTEGER                                 :: i

!  Start of executable code
      vmec_get_plasma_edge = MAX(ns, 2*mpol + 6)

      IF (.not.ASSOCIATED(r)) THEN
         ALLOCATE(r(vmec_get_plasma_edge), z(vmec_get_plasma_edge))
      END IF
      ALLOCATE(cosv(SIZE(xn, 1)), sinv(SIZE(xn, 1)))
      ALLOCATE(cosu(SIZE(xm, 1)), sinu(SIZE(xm, 1)))

      cosv(:) = COS(xn*(phi + this%phi_offset))
      sinv(:) = SIN(xn*(phi + this%phi_offset))
      DO i = 1, vmec_get_plasma_edge
         cosu(:) = COS(xm*(i - 1)*twopi/vmec_get_plasma_edge)
         sinu(:) = SIN(xm*(i - 1)*twopi/vmec_get_plasma_edge)

         r(i) = SUM(rmnc(:, ns)*(cosu(:)*cosv(:) + sinu(:)*sinv(:)))
         z(i) = SUM(zmns(:, ns)*(sinu(:)*cosv(:) - cosu(:)*sinv(:)))

         IF (lasym) THEN
            r(i) = r(i) +                                                      &
     &             SUM(rmns(:, ns)*(sinu(:)*cosv(:) - cosu(:)*sinv(:)))
            z(i) = z(i) +                                                      &
     &             SUM(zmnc(:, ns)*(cosu(:)*cosv(:) + sinu(:)*sinv(:)))
         END IF
      END DO

      DEALLOCATE(cosu, cosv, sinu, sinv)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration radial grid points.
!>
!>  This method overrides
!>  @ref equilibrium::equilibrium_get_magnetic_volume_rgrid.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_magnetic_volume_rgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   vmec_get_magnetic_volume_rgrid
      TYPE (vmec_class), INTENT(in) :: this

!  Start of executable code
      vmec_get_magnetic_volume_rgrid => this%magnetic_cache%rsuv

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration z grid points.
!>
!>  This method overrides
!>  @ref equilibrium::equilibrium_get_magnetic_volume_zgrid.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_magnetic_volume_zgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   vmec_get_magnetic_volume_zgrid
      TYPE (vmec_class), INTENT(in) :: this

!  Start of executable code
      vmec_get_magnetic_volume_zgrid => this%magnetic_cache%zsuv

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration jr grid points.
!>
!>  This method overrides
!>  @ref equilibrium::equilibrium_get_magnetic_volume_jrgrid.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_magnetic_volume_jrgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   vmec_get_magnetic_volume_jrgrid
      TYPE (vmec_class), INTENT(in) :: this

!  Start of executable code
      vmec_get_magnetic_volume_jrgrid => this%magnetic_cache%jrsuv

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration jphi grid points.
!>
!>  This method overrides
!>  @ref equilibrium::equilibrium_get_magnetic_volume_jphigrid.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_magnetic_volume_jphigrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   vmec_get_magnetic_volume_jphigrid
      TYPE (vmec_class), INTENT(in) :: this

!  Start of executable code
      vmec_get_magnetic_volume_jphigrid => this%magnetic_cache%jphisuv

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration jz grid points.
!>
!>  This method overrides
!>  @ref equilibrium::equilibrium_get_magnetic_volume_jzgrid.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_magnetic_volume_jzgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   vmec_get_magnetic_volume_jzgrid
      TYPE (vmec_class), INTENT(in) :: this

!  Start of executable code
      vmec_get_magnetic_volume_jzgrid => this%magnetic_cache%jzsuv

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume integration element.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_volume_int_element.
!>  For vmec this is J*ds*du*dv where J is the jacobian. The jacobian is
!>  integrated into Jr, Jphi and Jz but still need to account for the sign of
!>  the jacobian.
!>
!>  * ds = 1/(ns - 1)
!>  * du = 2*Pi/number of flux surfaces
!>  * dv = 2*Pi/number of field periods/number of toroidal planes
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The volume integration element.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_volume_int_element(this)
      USE read_wout_mod, only: isigng

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: vmec_get_volume_int_element
      TYPE (vmec_class), INTENT(in) :: this

!  Start of executable code
      vmec_get_volume_int_element = isigng*this%magnetic_cache%ds              &
     &                            * this%magnetic_cache%du                     &
     &                            * this%magnetic_cache%dv

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get external current.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_ext_currents. The
!>  array memory containing the external currents is owned by VMEC. Return a
!>  pointer to it. if VMEC is run in fixed boundary mode, return a null()
!>  pointer.
!>
!>  @param[in]  this           A @ref vmec_class instance.
!>  @param[out] scale_currents Informs the caller that currents need to be
!>                             scaled.
!>  @returns The external currents.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_ext_currents(this, scale_currents)
      USE vmec_input, only: lfreeb
      USE mgrid_mod, only: mgrid_mode, nextcur
      USE vmec_input, only: extcur

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:), POINTER :: vmec_get_ext_currents
      TYPE (vmec_class), INTENT(in)       :: this
      LOGICAL, INTENT(out)                :: scale_currents

!  Start of executable code
      IF (lfreeb) THEN
         SELECT CASE (mgrid_mode)

            CASE ('S', 's', 'N', 'n')
               scale_currents = .true.

            CASE DEFAULT
               scale_currents = .false.

         END SELECT
         vmec_get_ext_currents => extcur(1:nextcur)
      ELSE
         vmec_get_ext_currents => null()
      END IF

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Checks if a parameter id is a scaler value.
!>
!>  This method overrides @ref equilibrium::equilibrium_is_scaler_value.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the parameter is a scaler and false if otherwise.
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_scaler_value(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_scaler_value
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_bloat_id, vmec_curtor_id, vmec_phiedge_id,                 &
     &         vmec_pres_scale_id, vmec_vvc_smaleli_id,                        &
     &         vmec_vvc_kappa_p_id, vmec_betatot_id, vmec_betapol_id,          &
     &         vmec_betator_id, vmec_betaxis_id, vmec_phi_offset_id)
            vmec_is_scaler_value = .true.

         CASE DEFAULT
            vmec_is_scaler_value = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if a parameter id is a 1d array.
!>
!>  This method overrides @ref equilibrium::equilibrium_is_1d_array.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the parameter is a 1d array and false if otherwise.
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_1d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_1d_array
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_ac_id, vmec_ac_aux_s_id, vmec_ac_aux_f_id,                 &
     &         vmec_ai_id, vmec_ai_aux_s_id, vmec_ai_aux_f_id,                 &
     &         vmec_am_id, vmec_am_aux_s_id, vmec_am_aux_f_id,                 &
     &         vmec_pp_ne_b_id, vmec_pp_ne_as_id, vmec_pp_ne_af_id,            &
     &         vmec_pp_te_b_id, vmec_pp_te_as_id, vmec_pp_te_af_id,            &
     &         vmec_extcur_id, vmec_phi_id, vmec_iotaf_id,                     &
     &         vmec_iotas_id)
            vmec_is_1d_array = .true.

         CASE DEFAULT
            vmec_is_1d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if a parameter id is a 2d array.
!>
!>  This method overrides @ref equilibrium::equilibrium_is_2d_array.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the parameter is a 2d array and false if otherwise.
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_2d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_2d_array
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_rbc_id, vmec_zbs_id, vmec_rbs_id, vmec_zbc_id,             &
     &         vmec_rmnc_id, vmec_zmns_id, vmec_lmns_id, vmec_gmnc_id,         &
     &         vmec_bsubumnc_id, vmec_bsubvmnc_id,                             &
     &         vmec_bsubsmns_id, vmec_bsupumnc_id, vmec_bsupvmnc_id,           &
     &         vmec_rmns_id, vmec_zmnc_id, vmec_lmnc_id, vmec_gmns_id,         &
     &         vmec_bsubumns_id, vmec_bsubvmns_id,                             &
     &         vmec_bsubsmnc_id, vmec_bsupumns_id, vmec_bsupvmns_id,           &
     &         vmec_pp_sxrem_b_id, vmec_pp_sxrem_as_id,                        &
     &         vmec_pp_sxrem_af_id)
            vmec_is_2d_array = .true.

         CASE DEFAULT
            vmec_is_2d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if a parameter id is a reconstruction parameter.
!>
!>  This method overrides @ref equilibrium::equilibrium_is_recon_param.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the parameter is a reconstruction parameter and false if
!>  otherwise.
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_recon_param(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_recon_param
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_rbc_id, vmec_zbs_id, vmec_rbs_id, vmec_zbc_id,             &
     &         vmec_ac_id, vmec_ac_aux_s_id, vmec_ac_aux_f_id,                 &
     &         vmec_ai_id, vmec_ai_aux_s_id, vmec_ai_aux_f_id,                 &
     &         vmec_am_id, vmec_am_aux_s_id, vmec_am_aux_f_id,                 &
     &         vmec_pp_ne_b_id, vmec_pp_ne_as_id, vmec_pp_ne_af_id,            &
     &         vmec_pp_te_b_id, vmec_pp_te_as_id, vmec_pp_te_af_id,            &
     &         vmec_pp_sxrem_b_id, vmec_pp_sxrem_as_id,                        &
     &         vmec_pp_sxrem_af_id, vmec_extcur_id, vmec_curtor_id,            &
     &         vmec_phiedge_id, vmec_pres_scale_id, vmec_bloat_id,             &
     &         vmec_phi_offset_id)
            vmec_is_recon_param = .true.

         CASE DEFAULT
            vmec_is_recon_param = .false.

      END SELECT

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Solves the VMEC equilibrium.
!>
!>  This method overrides @ref equilibrium::equilibrium_converge. Solves the
!>  vmec equilibrium and loads the resulting wout file in preparation of
!>  computing modeled signals.
!>
!>  @param[inout] this     A @ref vmec_class instance.
!>  @param[inout] num_iter Counter to track the number of iterations.
!>  @param[in]    iou      Input/Output unit of the runlog file.
!>  @returns True if the convergece was sucessful and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION vmec_converge(this, num_iter, iou)
      USE vmec_input, only: delt, niter
      USE vmec_params, only: timestep_flag, output_flag,                       &
     &                       reset_jacdt_flag, jac75_flag,                     &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag
      USE xstuff, only: xcdot, xc
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE vmec_main, only: iterc, iter1, iter2

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                          :: vmec_converge
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(inout)           :: num_iter
      INTEGER, INTENT(in)              :: iou

!  local variables
!  delt_local   Value of delt on entry. Store, for restoration after error
!               recovery.
!  ictrl_array  Array containing various vmec control parameters.
      REAL (rprec)                     :: delt_local
      INTEGER, DIMENSION(5)            :: ictrl_array
      INTEGER                          :: recovery_index
      INTEGER                          :: error, num_iter_init
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: xc_saved

!  Start of executable code
      delt_local = delt
      num_iter_init = iterc
      recovery_index = 0

!  Store the value of the xc array. This xc array is different than the one
!  stored in save/reset state methods. This xc array needs to be stored so that
!  error recovery store the state to before runvmec was run.
      ALLOCATE(xc_saved(SIZE(xc)))
      xc_saved = xc

!  Run vmec until successful or an unrecoverable error has occured.
      DO
         ictrl_array = 0
         ictrl_array(1) = IOR(ictrl_array(1), timestep_flag)
         ictrl_array(1) = IOR(ictrl_array(1), output_flag)
         ictrl_array(1) = IOR(ictrl_array(1), reset_jacdt_flag)
         ictrl_array(3) = niter
         ictrl_array(4) = this%ns_index

!  Update internal vmec counters. These affect how the minimation occurs.
         iter1 = iter2

!  Zero the xcdot array.
         xcdot = 0.0
         CALL runvmec(ictrl_array, this%vmec_file_name, .false.)

!  Set error flag. (Flags declared in vmec_params.f)
!
!  Check error flags for error recovery modes.
!  If a recoverable error has occured, rerun vmec with altered eq_state. Error
!  recovery is attempted a maximum of 2 times.
!  If a non-recoverable error has occured, do not rerun vmec and report error.
!  If vmec terminated sucessfully continue.
         SELECT CASE (ictrl_array(2))
            CASE (jac75_flag)

!  Restore the xc array to the state it was before runvmec was called.
               xc = xc_saved

               IF (recovery_index .le. 2) THEN
                  WRITE (*,*) 'VMEC CONVERGENCE ERROR: jac75_flag'
                  WRITE (*,*) 'delt chnage: old, new:', delt, delt/2.0
                  WRITE (iou,*) 'VMEC CONVERGENCE ERROR: jac75_flag'
                  WRITE (iou,*) 'delt chnage: old, new:', delt, delt/2.0
                  delt = delt/2.0
               ELSE
                  EXIT
               END IF

            CASE DEFAULT
               EXIT

         END SELECT

         WRITE (*,*) 'Attempting error recovery'
         WRITE (iou,*) 'Attempting error recovery'

         recovery_index = recovery_index + 1
      END DO

!  No longer need the xc_saved array.
      DEALLOCATE(xc_saved)

!  Check for errors, if vmec was successful, load vmec data for use with the
!  signals. Other wise report the error.
      vmec_converge = .false.
      SELECT CASE (ictrl_array(2))

         CASE (successful_term_flag, norm_term_flag, more_iter_flag)
            CALL read_wout_file(this%wout_file_name, error)
            IF (error .eq. 0) THEN
               CALL LoadRZL
               vmec_converge = .true.

               IF (ASSOCIATED(this%magnetic_cache)) THEN
                  CALL vmec_set_magnetic_cache(this)
               END IF
            END IF

         CASE DEFAULT
!>  @todo FIXME: Put in error reporting here. For now leave this blank until I
!>  can figure how v3fit decides if vmec has converged.

      END SELECT

!  Increment the number of iterations by the number of iterations used by
!  runvmec.
      num_iter = num_iter + iterc - num_iter_init
      delt = delt_local

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Save the internal state of the equilibrium.
!>
!>  Saves the VMEC xc array and the wout file.
!>
!>  @param[inout] this A @ref vmec_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_save_state(this)
      USE xstuff, only: xc

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%xc)) THEN
         IF (SIZE(xc) .eq. SIZE(this%xc)) THEN
            this%xc = xc
         ELSE
            DEALLOCATE(this%xc)
            ALLOCATE(this%xc(SIZE(xc)))
            this%xc = xc
         END IF
      ELSE
         ALLOCATE(this%xc(SIZE(xc)))
         this%xc = xc
      END IF

      IF(ASSOCIATED(this%wout_file)) THEN
         CALL wout_file_cache_save_state(this%wout_file)
      ELSE
         this%wout_file => wout_file_cache_construct()
      END IF

      IF(ASSOCIATED(this%threed1_file)) THEN
         CALL threed1_file_cache_save_state(this%threed1_file)
      ELSE
         this%threed1_file => threed1_file_cache_construct()
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Reset the internal state of the equilibrium.
!>
!>  Copies the saved VMEC xc array back to VMEC.
!>
!>  @param[inout] this A @ref vmec_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_reset_state(this)
      USE xstuff, only: xc

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%xc) .and. (SIZE(xc) .eq. SIZE(this%xc))) THEN
         xc = this%xc
      END IF

!  Reset the wout file.
      IF (ASSOCIATED(this%wout_file)) THEN
         CALL wout_file_cache_reset_state(this%wout_file)
      END IF
      IF (ASSOCIATED(this%magnetic_cache)) THEN
         CALL vmec_set_magnetic_cache_calc(this)
      END IF
!  Reset the threed1 file.
      IF (ASSOCIATED(this%threed1_file)) THEN
         CALL threed1_file_cache_reset_state(this%threed1_file)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out the equilibrium to an output file.
!>
!>  This method overrides @ref equilibrium::equilibrium_write.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] iou  Input/output unit of the output file.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_write(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: iou

!  local variables
      INTEGER                       :: i

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,*) 'Equilibrium Type : VMEC'
      CALL pprofile_write(this%ne, 'pp_ne', iou)

      DO i = 1, SIZE(this%sxrem)
         CALL pprofile_write(this%sxrem(i)%p, 'pp_sxrem', iou)
      END DO

      CALL pprofile_write(this%te, 'pp_te', iou)

      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!>  @page result_file_vmec VMEC Result File
!>
!>  @tableofcontents
!>  @section result_file_vmec_intro_sec Introduction
!>  This page documents the contents of a result NetCDF file contributed by the
!>  VMEC equilibrium. The remaining parts of the result file are documented in
!>  the @ref result_file_main page.
!>
!>  @section result_file_vmec_dim_sec Dimensions
!>  @header{Dimension, Description, Code Reference}
!>  @begin_table
!>     @item{vmec_nextcur,    Number of vmec external currents.,         mgrid_mod::nextcur}
!>     @item{vmec_boundary_m, Number of poloidal boundary coefficients., vmec_dim::mpol1}
!>     @item{vmec_boundary_n, Number of toroidal boundary coefficients., vmec_input::ntor}
!>  @end_table
!>  @section result_file_vmec_var_sec Variables
!>  @header{Variable(Dimensions), Description, Code Reference}
!>  @begin_table
!>     @item{vmec_rbc (vmec_boundary_n\, vmec_boundary_m), R direction cosine boundary coefficients., vmec_input::rbc}
!>     @item{vmec_rbs (vmec_boundary_n\, vmec_boundary_m), R direction sine boundary coefficients.,   vmec_input::rbs}
!>     @item{vmec_zbc (vmec_boundary_n\, vmec_boundary_m), Z direction cosine boundary coefficients., vmec_input::zbc}
!>     @item{vmec_zbs (vmec_boundary_n\, vmec_boundary_m), Z direction sine boundary coefficients.,   vmec_input::zbs}
!>     @item{vmec_l_sub_i,                                 Internal inducance.,                       v3f_vmec_comm::vvc_smaleli}
!>     @item{vmec_kappa_p,                                 Mean elongation.,                          v3f_vmec_comm::vvc_kappa_p}
!>  @end_table
!-------------------------------------------------------------------------------
!>  @brief Define NetCDF variables for the result file
!>
!>  This method overrides @ref equilibrium::equilibrium_def_result. Defines
!>  dimensions and variables for the VMEC contribution of the result file. Multi
!>  dimensional arrays need to be transposed so arrays appear in the correct
!>  order in non fortran languages.
!>
!>  @param[in] this             A @ref vmec_class instance.
!>  @param[in] result_ncid      NetCDF file id of the result file.
!>  @param[in] maxnsetps_dim_id NetCDF dimension id of the number of steps
!>                              dimension.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_def_result(this, result_ncid, maxnsetps_dim_id)
      USE mgrid_mod, only: nextcur
      USE vmec_dim, only: mpol1
      USE vmec_input, only: ntor
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: result_ncid
      INTEGER, INTENT(in)           :: maxnsetps_dim_id

!  Local variables
      INTEGER                       :: status
      INTEGER                       :: vmec_nextcur_dim_id
      INTEGER                       :: vmec_boundary_m_dim_id
      INTEGER                       :: vmec_boundary_n_dim_id

!  Start of executable code
!  Define dimensions.
      status = nf_def_dim(result_ncid, 'vmec_nextcur', nextcur,                &
     &                    vmec_nextcur_dim_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_dim(result_ncid, 'vmec_boundary_m', mpol1 + 1,           &
     &                    vmec_boundary_m_dim_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_dim(result_ncid, 'vmec_boundary_n', 2*ntor + 1,          &
     &                    vmec_boundary_n_dim_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

!  Define variables.
      status = nf_def_var(result_ncid, 'vmec_extcur', nf_double, 1,            &
     &                    (/ vmec_nextcur_dim_id /),                           &
     &                    this%vmec_extcur_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_curtor', nf_double, 0, 0,         &
     &                    this%vmec_curtor_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_pres_scale', nf_double, 0,        &
     &                    0, this%vmec_pres_scale_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_rbc', nf_double, 2,               &
     &                    (/ vmec_boundary_n_dim_id,                           &
     &                       vmec_boundary_m_dim_id /),                        &
     &                    this%vmec_rbc_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_rbs', nf_double, 2,               &
     &                    (/ vmec_boundary_n_dim_id,                           &
     &                       vmec_boundary_m_dim_id /),                        &
     &                    this%vmec_rbs_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_zbc', nf_double, 2,               &
     &                    (/ vmec_boundary_n_dim_id,                           &
     &                       vmec_boundary_m_dim_id /),                        &
     &                    this%vmec_zbc_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_zbs', nf_double, 2,               &
     &                    (/ vmec_boundary_n_dim_id,                           &
     &                       vmec_boundary_m_dim_id /),                        &
     &                    this%vmec_zbs_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_l_sub_i', nf_double, 1,           &
     &                    (/ maxnsetps_dim_id /), this%vmec_l_sub_i_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_kappa_p', nf_double, 1,           &
     &                    (/ maxnsetps_dim_id /), this%vmec_kappa_p_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write inital data to NetCDF result file
!>
!>  This method overrides @ref equilibrium::equilibrium_write_init_data.
!>
!>  @param[in] this        A @ref vmec_class instance.
!>  @param[in] result_ncid NetCDF file id of the result file.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_write_init_data(this, result_ncid)
      USE mgrid_mod, only: nextcur
      USE vmec_input, only: extcur, curtor, pres_scale, ntor,                  &
     &                      rbc, zbs, rbs, zbc
      USE vmec_dim, only: mpol1
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: result_ncid

!  Local variables
      INTEGER                       :: status
      INTEGER, DIMENSION(2)         :: boundary_size

!  Start of executable code
      status = nf_put_var_double(result_ncid, this%vmec_extcur_id,             &
     &                           extcur(1:nextcur))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_var_double(result_ncid, this%vmec_curtor_id,             &
     &                           curtor)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_var_double(result_ncid, this%vmec_pres_scale_id,         &
     &                           pres_scale)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      boundary_size(1) = 2*ntor + 1
      boundary_size(2) = mpol1 + 1

      status = nf_put_vara_double(result_ncid, this%vmec_rbc_id,               &
     &                            (/ 1, 1 /), boundary_size,                   &
     &                            rbc(-ntor:ntor,0:mpol1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_vara_double(result_ncid, this%vmec_rbs_id,               &
     &                            (/ 1, 1 /), boundary_size,                   &
     &                            rbs(-ntor:ntor,0:mpol1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_vara_double(result_ncid, this%vmec_zbc_id,               &
     &                            (/ 1, 1 /), boundary_size,                   &
     &                            zbc(-ntor:ntor,0:mpol1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_vara_double(result_ncid, this%vmec_zbs_id,               &
     &                            (/ 1, 1 /), boundary_size,                   &
     &                            zbs(-ntor:ntor,0:mpol1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      CALL vmec_write_step_data(this, result_ncid, 1)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write step data to NetCDF result file
!>
!>  This method overrides @ref equilibrium::equilibrium_write_step_data.
!>
!>  @param[in] this         A @ref vmec_class instance.
!>  @param[in] result_ncid  NetCDF file id of the result file.
!>  @param[in] current_step Step index to write variables to.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_write_step_data(this, result_ncid, current_step)
      USE v3f_vmec_comm
      USE ezcdf

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: result_ncid
      INTEGER, INTENT(in)           :: current_step

!  Local variables
      INTEGER                       :: status

!  Start of executable code
      status = nf_put_var1_double(result_ncid, this%vmec_l_sub_i_id,           &
     &                            current_step, vvc_smaleli)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_var1_double(result_ncid, this%vmec_kappa_p_id,           &
     &                            current_step, vvc_kappa_p)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Reset the fixed boundary coefficients.
!>
!>  When changing the boundary coefficients, there is extra work that must be
!>  performed to set the VMEC state.
!>
!>  @param[in] this A @ref vmec_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_reset_boundary(this)
      USE vmec_input, only: lasym, rbc, zbs, rbs, zbc, ntor,                   &
     &                      mfilter_fbdy, nfilter_fbdy, lfreeb
      USE vmec_dim, only: mpol1, ntor1
      USE vmec_main, only: rmn_bdy, zmn_bdy, lthreed, lconm1
      USE vmec_params, only: rcc, rss, rsc, rcs,                               &
     &                       zsc, zcs, zcc, zss, signgs
      USE vparams, only: cp5

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in)        :: this

!  local variables
      REAL (rprec)                            :: delta, sgn
      INTEGER                                 :: m, n, mj, ni
      REAL (rprec)                            :: temp
      REAL (rprec), DIMENSION(:,:), POINTER   :: rbcc, rbss, rbcs, rbsc
      REAL (rprec), DIMENSION(:,:), POINTER   :: zbcc, zbss, zbcs, zbsc
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: temp_array

!  Start of executable code
      IF (lasym) THEN
!  Convert to representation with rbs(m=1) = zbc(m=1)
         delta = ATAN((rbs(0,1) - zbc(0,1))/                                   &
     &                (ABS(rbc(0,1) + ABS(zbs(0,1)))))
         IF (delta .ne. 0.0) THEN
            DO m = 0, mpol1
               DO n = -ntor, ntor
                  temp = rbc(n,m)*COS(m*delta) + rbs(n,m)*SIN(m*delta)
                  rbs(n,m) = rbs(n,m)*COS(m*delta)                             &
     &                     - rbc(n,m)*SIN(m*delta)
                  rbc(n,m) = temp
                  temp = zbc(n,m)*COS(m*delta) + zbs(n,m)*SIN(m*delta)
                  zbs(n,m) = zbs(n,m)*COS(m*delta)                             &
     &                     - zbc(n,m)*SIN(m*delta)
                  zbc(n,m) = temp
               END DO
            END DO
         END IF
      END IF

!  Convert to internal representation of modes
!
!  r = rbcc*cos(m*u)*cos(n*v) + rbss*sin(m*u)*sin(n*v)
!    + rbcs*cos(m*u)*sin(n*v) + rbsc*sin(m*u)*cos(n*v)
!  z = zbcs*cos(m*u)*sin(n*v) + zbsc*sin(m*u)*cos(n*v)
!    + zbcc*cos(m*u)*cos(n*v) + zbsz*sin(m*u)*sin(n*v)

      rbcc => rmn_bdy(:,:,rcc)
      zbsc => zmn_bdy(:,:,zsc)
      IF (lthreed) THEN
         rbss => rmn_bdy(:,:,rss)
         zbcs => zmn_bdy(:,:,zcs)
      END IF

      IF (lasym) THEN
         rbsc => rmn_bdy(:,:,rsc)
         zbcc => rmn_bdy(:,:,zcc)
         IF (lthreed) THEN
            rbcs => rmn_bdy(:,:,rcs)
            zbss => rmn_bdy(:,:,zss)
         END IF
      END IF

      rmn_bdy = 0.0
      zmn_bdy = 0.0

      IF (.not.lfreeb) THEN
         DO m = 0, mpol1
            IF ((mfilter_fbdy .gt. 1) .and. (m .gt. mfilter_fbdy)) THEN
               EXIT
            END IF

            mj = m + LBOUND(rbcc, 2)
            DO n = -ntor, ntor
               IF ((nfilter_fbdy .gt. 0) .and.                                 &
     &             (ABS(n) .gt. nfilter_fbdy)) THEN
                  CYCLE
               END IF

               ni = ABS(n) + LBOUND(rbcc, 1)

               IF (n .eq. 0) THEN
                  sgn = 0.0
               ELSE IF (n .gt. 0) THEN
                  sgn = 1.0
               ELSE
                  sgn = -1.0
               END IF

               rbcc(ni,mj) = rbcc(ni,mj) + rbc(n,m)
               IF (m .gt. 0) THEN
                  zbsc(ni,mj) = zbsc(ni,mj) + zbs(n,m)
               END IF

               IF (lthreed) THEN
                  zbcs(ni,mj) = zbcs(ni,mj) - sgn*zbs(n,m)
                  IF (m .gt. 0) THEN
                     rbss(ni,mj) = rbss(ni,mj) + sgn*rbc(n,m)
                  END IF
               END IF

               IF (lasym) THEN
                  zbcc(ni,mj) = zbcc(ni,mj) + zbc(n,m)
                  IF (m .gt. 0) THEN
                     rbsc(ni,mj) = rbsc(ni,mj) + rbs(n,m)
                  END IF

                  IF (lthreed) THEN
                     rbcs(ni,mj) = rbcs(ni,mj) - sgn*rbs(n,m)
                     IF (m .gt. 0) THEN
                        zbss(ni,mj) = zbss(ni,mj) + sgn*zbc(n,m)
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END IF

!  Check sign of jacobian (should be same as signgs)
      mj = 1 + LBOUND(rbcc, 2)
      signgs = 1.0
      IF (SUM(rbcc(1:ntor1,mj))*SUM(zbsc(1:ntor1,mj)) .gt. 0.0) THEN
         signgs = -1.0
      END IF

!  Convert to internal form for (constrained) m=1 modes
      IF (lconm1 .and. (lthreed .or. lasym)) THEN
         ALLOCATE(temp_array(SIZE(rbcc, 1)))
         IF (lthreed) THEN
            temp_array = rbss(:,mj)
            rbss(:,mj) = cp5*(temp_array + zbcs(:,mj))
            zbcs(:,mj) = cp5*(temp_array - zbcs(:,mj))
         END IF

         IF (lasym) THEN
            temp_array = rbsc(:,mj)
            rbsc(:,mj) = cp5*(temp_array + zbcc(:,mj))
            zbcc(:,mj) = cp5*(temp_array - zbcc(:,mj))
         END IF
         DEALLOCATE(temp_array)
      END IF

      END SUBROUTINE

      END MODULE
