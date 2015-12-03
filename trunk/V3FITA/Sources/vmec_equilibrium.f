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
!>     @item{ac,         1D Array of function current profile parameters.,  @ref vmec_input::ac}
!>     @item{ac_aux_s,   1D Array of segment current profile s poitions.,   @ref vmec_input::ac_aux_s}
!>     @item{ac_aux_f,   1D Array of segment current profile f values.,     @ref vmec_input::ac_aux_f}
!>     @item{ai,         1D Array of function current profile parameters.,  @ref vmec_input::ai}
!>     @item{ai_aux_s,   1D Array of segment iotabar profile s poitions.,   @ref vmec_input::ai_aux_s}
!>     @item{ai_aux_f,   1D Array of segment iotabar profile f values.,     @ref vmec_input::ai_aux_f}
!>     @item{am,         1D Array of function current profile parameters.,  @ref vmec_input::am}
!>     @item{am_aux_s,   1D Array of segment pressure profile s poitions.,  @ref vmec_input::am_aux_s}
!>     @item{am_aux_f,   1D Array of segment pressure profile f values.,    @ref vmec_input::am_aux_f}
!>     @item{bloat,      Extends the field into the vacuum region.,         @ref vmec_input::bloat}
!>     @item{rbc,        2D Array of outer boundary R cosine coeffiecents., @ref vmec_input::rbc}
!>     @item{zbs,        2D Array of outer boundary Z sine coeffiecents.,   @ref vmec_input::zbs}
!>     @item{rbs,        2D Array of outer boundary R sine coeffiecents.,   @ref vmec_input::rbs}
!>     @item{zbc,        2D Array of outer boundary Z cosine coeffiecents., @ref vmec_input::zbc}
!>     @item{extcur,     1D Array of external currents.,                    @ref vmec_input::extcur}
!>     @item{curtor,     Toroidal current.,                                 @ref vmec_input::curtor}
!>     @item{phiedge,    The edge magnetic flux.,                           @ref vmec_input::phiedge}
!>     @item{pres_scale, The scaling parameter for the pressure.,           @ref vmec_input::pres_scale}
!>  @end_table
!>  @subsection vmec_equilibrium_aux_recon_param_sec VMEC Auxiliary Reconstruction Parameters
!>  VMEC Auxiliary reconstruction parameters are parameters that maybe
!>  reconstructed that do not change the equilibirum. These parameters do not
!>  require the equilibrium to reconverge.
!>  @begin_table
!>     @item{pp_ne_b,       1D Array of function density profile parameters.,     @ref pprofile_T::pprofile_class::b}
!>     @item{pp_ne_as,      1D Array of segment density profile s poitions.,      @ref pprofile_T::pprofile_class::as}
!>     @item{pp_ne_af,      1D Array of segment density profile f values.,        @ref pprofile_T::pprofile_class::af}
!>     @item{pp_sxrem_b,    1D Array of function sxrem profile parameters.
!>                          DEPREICATED only use id when using the old profile
!>                          specification.,                                       @ref pprofile_T::pprofile_class::b}
!>     @item{pp_sxrem_as,   1D Array of segment sxrem profile s poitions.
!>                          DEPREICATED only use id when using the old profile
!>                          specification.,                                       @ref pprofile_T::pprofile_class::as}
!>     @item{pp_sxrem_af,   1D Array of segment sxrem profile f values.
!>                          DEPREICATED only use id when using the old profile
!>                          specification.,                                       @ref pprofile_T::pprofile_class::af}
!>     @item{pp_sxrem_b_a,  2D Array of function sxrem profile parameters.
!>                          First index is the profile number.,                   @ref vmec_equilibrium::vmec_class::sxrem}
!>     @item{pp_sxrem_as_a, 2D Array of segment sxrem profile s poitions.
!>                          First index is the profile number.,                   @ref vmec_equilibrium::vmec_class::sxrem}
!>     @item{pp_sxrem_af_a, 2D Array of segment sxrem profile f values.
!>                          First index is the profile number.,                   @ref vmec_equilibrium::vmec_class::sxrem}
!>     @item{pp_te_b,       1D Array of function temperature profile parameters., @ref pprofile_T::pprofile_class::b}
!>     @item{pp_te_as,      1D Array of segment temperature profile s poitions.,  @ref pprofile_T::pprofile_class::as}
!>     @item{pp_te_af,      1D Array of segment temperature profile f values.,    @ref pprofile_T::pprofile_class::af}
!>     @item{phi_offset,    Phi angle offset in radians,                          @ref vmec_equilibrium::vmec_class::phi_offset}
!>     @item{z_offset,      Z offset in meters,                                   @ref vmec_equilibrium::vmec_class::z_offset}
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
      USE vmec_context

      IMPLICIT NONE

!*******************************************************************************
!  vmec equilibrium module parameters
!*******************************************************************************
!  Define id's only for values that can change. All others shouldn't be needed
!  outside of the vmec interface. Parameter id's start at 7 since 0-6 are
!  reserved for non equilibrium model parameters. These numbers will need to be
!  updated if any new model parameters are added.

!  Variable Parameters
!>  1D Array of function current profile parameters. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_id          = 7
!>  1D Array of segment current profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_aux_s_id    = 8
!>  1D Array of segment current profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_aux_f_id    = 9
!>  1D Array of function iotabar profile parameters. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_id          = 10
!>  1D Array of segment iotabar profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_aux_s_id    = 11
!>  1D Array of segment iotabar profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_aux_f_id    = 12
!>  1D Array of function pressure profile parameters. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_id          = 13
!>  1D Array of segment pressure profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_aux_s_id    = 14
!>  1D Array of segment pressure profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_aux_f_id    = 15
!>  Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_bloat_id       = 16
!>  2D Array of outer boundary R cosine coeffiecents. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_rbc_id         = 17
!>  2D Array of outer boundary Z sine coeffiecents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_zbs_id         = 18
!>  2D Array of outer boundary R sine coeffiecents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_rbs_id         = 19
!>  2D Array of outer boundary Z cosine coeffiecents. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_zbc_id         = 20
!>  1D Array of external currents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_extcur_id      = 21
!>  The toroidal current. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_curtor_id      = 22
!>  The edge magnetic flux. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_phiedge_id     = 23
!>  The scaling parameter for the pressure. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_pres_scale_id  = 24

!  Derived Parameters
!>  2D Array of surface R cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_rmnc_id        = 25
!>  2D Array of surface Z sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_zmns_id        = 26
!>  2D Array of surface lambda cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_lmns_id        = 27
!>  2D Array of surface g cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_gmnc_id        = 28
!>  2D Array of surface B_s sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubsmns_id    = 29
!>  2D Array of surface B_u cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubumnc_id    = 30
!>  2D Array of surface B_v cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubvmnc_id    = 31
!>  2D Array of surface B^u cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupumnc_id    = 32
!>  2D Array of surface B^v cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupvmnc_id    = 33
!>  2D Array of surface R sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_rmns_id        = 34
!>  2D Array of surface Z cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_zmnc_id        = 35
!>  2D Array of surface lambda sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_lmnc_id        = 36
!>  2D Array of surface g sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_gmns_id        = 37
!>  2D Array of surface B_s cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubsmnc_id    = 38
!>  2D Array of surface B_u sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubumns_id    = 39
!>  2D Array of surface B_v sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubvmns_id    = 40
!>  2D Array of surface B^u sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupumns_id    = 41
!>  2D Array of surface B^v sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupvmns_id    = 42
!>  Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_phi_id         = 43
!>  1D Array of full grid iotabar. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_iotaf_id       = 44
!>  1D Array of half grid iotabar. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_iotas_id       = 45

!  Auxiliary Model Parameters
!>  1D Array of function density profile parameters. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_ne_b_id     = 46
!>  1D Array of segment density profile s poitions. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_ne_as_id    = 47
!>  1D Array of segment density profile f values. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_ne_af_id    = 48
!>  1D Array of function sxrem profile parameters. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_sxrem_b_id  = 49
!>  1D Array of segment sxrem profile s poitions. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_sxrem_as_id = 50
!>  1D Array of segment sxrem profile f values. Defined in @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_sxrem_af_id = 51
!>  1D Array of function temperature profile parameters. Defined in
!>  @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_te_b_id     = 52
!>  1D Array of segment temperature profile s poitions. Defined in
!>  @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_te_as_id    = 53
!>  1D Array of segment temperature profile f values. Defined in
!>  @ref pprofile_T.
      INTEGER, PARAMETER :: vmec_pp_te_af_id    = 54

!  Derived Parameters
!>  Internal inducance. Defined in @ref v3f_vmec_comm. Computed in @ref eqfor.
      INTEGER, PARAMETER :: vmec_vvc_smaleli_id = 55
!>  Mean elongation. Defined in @ref v3f_vmec_comm. Computed in @ref eqfor.
      INTEGER, PARAMETER :: vmec_vvc_kappa_p_id = 56
!>  Total Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betatot_id     = 57
!>  Poloidal Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betapol_id     = 58
!>  Toroidal Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betator_id     = 59
!>  Beta on axis. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betaxis_id     = 60

!>  Plasma Phi offset. This is a parameter to allow changing the phi angle of a
!>  quasi helical state in an RFP. Defined in @ref vmec_equilibrium.
      INTEGER, PARAMETER :: vmec_phi_offset_id  = 61
!>  Plasma z offset. This is a parameter to allow changing the vertical shift of
!>  the plasma
      INTEGER, PARAMETER :: vmec_z_offset_id    = 62

!>  Virtual casing grid points for vacuum field measurements.
      INTEGER, PARAMETER :: magnetic_cache_vc_min_grid_points = 101
!>  Virtual casing grid dimension.
      REAL (rprec), PARAMETER :: magnetic_cache_vc_grid_dim = 0.01

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
!>  Volume differential element for the s direction.
         REAL (rprec)                            :: ds
!>  Volume differential element for the u direction.
         REAL (rprec)                            :: du
!>  Volume differential element for the v direction.
         REAL (rprec)                            :: dv

!>  Area differential elements for the u direction.
         REAL (rprec)                            :: du_a
!>  Area differential elements for the v direction.
         REAL (rprec)                            :: dv_a

!>  Full surface area differential elements for the u direction.
         REAL (rprec)                            :: du_full
!>  Full surface area differential elements for the v direction.
         REAL (rprec)                            :: dv_full

!>  R positions on the suv gird points.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: rsuv => null()
!>  Z positions on the suv gird points.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: zsuv => null()

!>  Current density in the R direction. Has a 1/Sqrt(g) factor embeded.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jrsuv => null()
!>  Current density in the Phi direction. Has a 1/Sqrt(g) factor embeded.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jphisuv => null()
!>  Current density in the Z direction. Has a 1/Sqrt(g) factor embeded.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jzsuv => null()

!>  Surface current density in the R direction. Has all Jacobian and
!>  normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kruv => null()
!>  Surface current density in the Phi direction. Has all Jacobian and
!>  normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kphiuv => null()
!>  Surface current density in the Z direction. Has all Jacobian and
!>  normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kzuv => null()

!>  X' vector on the full uv grid.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: x_prime => null()
!>  Full surface current density in the R direction. Has all Jacobian and
!>  normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kxuv_full => null()
!>  Full surface current density in the Phi direction. Has all Jacobian and
!>  normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kyuv_full => null()
!>  Full surface current density in the Z direction. Has all Jacobian and
!>  normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kzuv_full => null()

!>  Axisymmtric X' vector on the full uv grid.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: x_axi => null()
!>  Axisymmtric  surface current density in the R direction. Has all Jacobian
!>  and normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kxuv_axi => null()
!>  Full surface current density in the Phi direction. Has all Jacobian and
!>  normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kyuv_axi => null()
!>  Full surface current density in the Z direction. Has all Jacobian and
!>  normalization factors.
         REAL (rprec), DIMENSION(:,:), POINTER :: kzuv_axi => null()
      END TYPE

!-------------------------------------------------------------------------------
!>  Base class representing a vmec_equilibrium.
!>  @par Super Class:
!>  @ref equilibrium
!-------------------------------------------------------------------------------
      TYPE vmec_class
!>  File name of the output of v3fit.
         CHARACTER (len=data_filename_length)   :: wout_file_name
!>  File name of the vmec namelist inout file.
         CHARACTER (len=data_filename_length)   :: vmec_file_name

!>  @ref pprofile_T describing a flux surface constant electron denisty profile.
         TYPE (pprofile_class), POINTER         :: ne => null()
!>  @ref pprofile_T describing a flux surface constant electron temperature
!>  profile.
         TYPE (pprofile_class), POINTER         :: te => null()
!>  @ref pprofile_T describing a flux surface conatsnt soft x-ray emissivity
!>  profile.
         TYPE (pprofile_pointer), DIMENSION(:), POINTER ::                     &
     &      sxrem => null()

!>  Ratio of the poloidal grid points to the poloidal grid points.
         REAL (rprec)                           :: pol_rad_ratio
!>  Instance of a @ref vmec_magnetic_cache object.
         TYPE (vmec_magnetic_cache), POINTER    ::                             &
     &      magnetic_cache => null()

!>  VMEC context.
         TYPE (vmec_context_class), POINTER     ::                             &
     &      vmec_context_save => null()

!  Extra reconstruction parameters
!>  Plasma phi offset.
         REAL (rprec)                           :: phi_offset = 0
!>  Plasma z offset.
         REAL (rprec)                           :: z_offset = 0

!>  Index of the largest ns value.
         INTEGER                                :: ns_index

      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the setting the @ref vmec_magnetic_cache types using
!>  @ref vmec_set_magnetic_cache_init or @ref vmec_set_magnetic_cache_calc. 
!-------------------------------------------------------------------------------
      INTERFACE vmec_set_magnetic_cache
         MODULE PROCEDURE vmec_set_magnetic_cache_responce,                    &
     &                    vmec_set_magnetic_cache_point,                       &
     &                    vmec_set_magnetic_cache_calc,                        &
     &                    vmec_set_magnetic_cache_sync
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

!-------------------------------------------------------------------------------
!>  Interface for the vmec density profile values.
!-------------------------------------------------------------------------------
      INTERFACE vmec_get_ne
         MODULE PROCEDURE vmec_get_ne_cart,                                    &
     &                    vmec_get_ne_radial
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for the vmec temperature profile values.
!-------------------------------------------------------------------------------
      INTERFACE vmec_get_te
         MODULE PROCEDURE vmec_get_te_cart,                                    &
     &                    vmec_get_te_radial
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for the vmec soft x-ray emissivity profile values.
!-------------------------------------------------------------------------------
      INTERFACE vmec_get_sxrem
         MODULE PROCEDURE vmec_get_sxrem_cart,                                 &
     &                    vmec_get_sxrem_radial
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for the vmec pressure profile values.
!-------------------------------------------------------------------------------
      INTERFACE vmec_get_p
         MODULE PROCEDURE vmec_get_p_cart,                                     &
     &                    vmec_get_p_radial
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref vmec_construct types using
!>  @ref vmec_construct_new or @ref vmec_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE vmec_construct
         MODULE PROCEDURE vmec_construct_new,                                  &
     &                    vmec_construct_mpi
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
!>  @param[in] z_offset   Initial Z offset of the plasma relative to the
!>                        machine center.
!>  @param[in] iou        Input/output unit to log messages to.
!>  @returns A pointer to a constructed @ref vmec_class object.
!-------------------------------------------------------------------------------
      FUNCTION vmec_construct_new(file_name, wout_file_name,                   &
     &                            ne, te, sxrem, phi_offset,                   &
     &                            z_offset, pol_rad_ratio, iou)
      USE vmec_params, only: restart_flag, readin_flag, timestep_flag,         &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag
      USE vmec_input, only: l_v3fit, ns_array
      USE read_wout_mod, only: read_wout_file, LoadRZL

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER           :: vmec_construct_new
      CHARACTER (len=*), INTENT(in)        :: file_name
      CHARACTER (len=*), INTENT(in)        :: wout_file_name
      TYPE (pprofile_class), POINTER       :: ne
      TYPE (pprofile_class), POINTER       :: te
      TYPE (pprofile_pointer), DIMENSION(:), POINTER :: sxrem
      REAL (rprec), INTENT(in)             :: phi_offset
      REAL (rprec), INTENT(in)             :: z_offset
      REAL (rprec), INTENT(in)             :: pol_rad_ratio
      INTEGER, INTENT(in)                  :: iou

!  local variables
      INTEGER, DIMENSION(5)                :: ictrl_array
      INTEGER                              :: index_dat, index_end
      INTEGER                              :: error
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (*,*) ' *** Initializing VMEC equilibrium from file ' //           &
     &            TRIM(file_name)
      WRITE (iou,*) ' *** Initializing VMEC equilibrium from file ' //         &
     &              TRIM(file_name)

      ALLOCATE(vmec_construct_new)
      vmec_construct_new%ne => ne
      vmec_construct_new%te => te
      vmec_construct_new%sxrem => sxrem

      vmec_construct_new%phi_offset = phi_offset
      vmec_construct_new%z_offset = z_offset

      vmec_construct_new%pol_rad_ratio = pol_rad_ratio

!  Initialize VMEC by taking a single step.
      ictrl_array = 0
      ictrl_array(1) = IOR(ictrl_array(1), restart_flag)
      ictrl_array(1) = IOR(ictrl_array(1), readin_flag)
      ictrl_array(1) = IOR(ictrl_array(1), timestep_flag)
      ictrl_array(3) = 1 ! Run for one iteration.
      ictrl_array(4) = 1 ! Limit VMEC to single grid.

      l_v3fit = .true.

      vmec_construct_new%vmec_file_name = TRIM(file_name)

      IF (wout_file_name .ne. '') THEN
!  runvmec deletes the string holding the wout file to read. As a result save
!  the wout filename before using it.
         vmec_construct_new%wout_file_name = wout_file_name
         CALL runvmec(ictrl_array, file_name, .false., wout_file_name)

         CALL read_wout_file(vmec_construct_new%wout_file_name, error)
         CALL assert_eq(error, 0, 'Error loading initialized ' //              &
     &                            'wout file.')
         CALL LoadRZL
      ELSE
         CALL runvmec(ictrl_array, file_name, .false.)
      END IF

!  Check for errors in runvmec.
      SELECT CASE (ictrl_array(2))

         CASE (successful_term_flag, norm_term_flag, more_iter_flag)
!  Set the wout file name.
            index_dat = INDEX(file_name, 'input.')
            index_end = LEN_TRIM(file_name)
            IF (index_dat .gt. 0) THEN
               vmec_construct_new%wout_file_name = 'wout_' //                  &
     &            TRIM(file_name(index_dat + 6:index_end)) // '.nc'
            ELSE
               vmec_construct_new%wout_file_name = 'wout_' //                  &
     &            TRIM(file_name(1:index_end)) // '.nc'
            END IF

         CASE DEFAULT
            CALL err_fatal('vmec_construct runvmec command failed ',           &
     &                     int = ictrl_array(2))

      END SELECT

      vmec_construct_new%vmec_context_save => vmec_context_construct()

      vmec_construct_new%ns_index = MAXLOC(ns_array, 1)

      CALL profiler_set_stop_time('vmec_construct', start_time)
      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref vmec_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref vmec_class object. Runs vmec
!>  for one iteration to initialize all the input parameters. Defined in
!>  @ref vmec_input.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref vmec_class object.
!-------------------------------------------------------------------------------
      FUNCTION vmec_construct_mpi(intra_comm)
      USE vmec_params, only: restart_flag, readin_flag, timestep_flag,         &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag
      USE vmec_input, only: l_v3fit, ns_array

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER :: vmec_construct_mpi
      INTEGER, INTENT(in)        :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                    :: error
      INTEGER, DIMENSION(5)      :: ictrl_array
      INTEGER                    :: index_dat, index_end
      INTEGER                    :: i, num_sxrem
      LOGICAL                    :: use_cache

!  Start of executable code
      ALLOCATE(vmec_construct_mpi)

      CALL MPI_BCAST(vmec_construct_mpi%phi_offset, 1, MPI_REAL8, 0,           &
     &               intra_comm, error)
      CALL MPI_BCAST(vmec_construct_mpi%z_offset, 1, MPI_REAL8, 0,             &
     &               intra_comm, error)

      CALL MPI_BCAST(vmec_construct_mpi%pol_rad_ratio, 1, MPI_REAL8, 0,        &
     &               intra_comm, error)

      CALL MPI_BCAST(vmec_construct_mpi%vmec_file_name,                        &
     &               data_filename_length, MPI_CHARACTER, 0,                   &
     &               intra_comm, error)

      vmec_construct_mpi%ne => pprofile_construct(intra_comm)
      vmec_construct_mpi%te => pprofile_construct(intra_comm)

      CALL MPI_BCAST(num_sxrem, 1, MPI_INTEGER, 0, intra_comm, error)
      ALLOCATE(vmec_construct_mpi%sxrem(num_sxrem))
      DO i = 1, num_sxrem
         vmec_construct_mpi%sxrem(i)%p => pprofile_construct(intra_comm)
      END DO

      CALL MPI_BCAST(vmec_construct_mpi%ns_index, 1, MPI_INTEGER, 0,           &
     &               intra_comm, error)

!  Initialize VMEC by taking a single step.
      ictrl_array = 0
      ictrl_array(1) = IOR(ictrl_array(1), restart_flag)
      ictrl_array(1) = IOR(ictrl_array(1), readin_flag)
      ictrl_array(1) = IOR(ictrl_array(1), timestep_flag)
      ictrl_array(3) = 1 ! Run for one iteration.
      ictrl_array(4) = vmec_construct_mpi%ns_index

      l_v3fit = .true.

      CALL runvmec(ictrl_array, vmec_construct_mpi%vmec_file_name,             &
     &             .false.)

!  Check for errors in runvmec.
      SELECT CASE (ictrl_array(2))

         CASE (successful_term_flag, norm_term_flag, more_iter_flag)
!  Set the wout file name.
            index_dat = INDEX(vmec_construct_mpi%vmec_file_name,               &
     &                        'input.')
            index_end = LEN_TRIM(vmec_construct_mpi%vmec_file_name)
            IF (index_dat .gt. 0) THEN
               vmec_construct_mpi%wout_file_name = 'wout_' //                  &
     &            TRIM(vmec_construct_mpi%vmec_file_name(                      &
     &                    index_dat + 6:index_end)) // '.nc'
            ELSE
               vmec_construct_mpi%wout_file_name = 'wout_' //                  &
     &            TRIM(vmec_construct_mpi%vmec_file_name(                      &
     &                    1:index_end)) // '.nc'
            END IF

         CASE DEFAULT
            CALL err_fatal('vmec_construct runvmec command failed ',           &
     &                     int = ictrl_array(2))

      END SELECT

      vmec_construct_mpi%vmec_context_save => vmec_context_construct()

      vmec_construct_mpi%ns_index = MAXLOC(ns_array, 1)

      CALL MPI_BCAST(use_cache, 1, MPI_LOGICAL, 0, intra_comm, error)
      IF (use_cache) THEN
         ALLOCATE(vmec_construct_mpi%magnetic_cache)
         CALL vmec_set_magnetic_cache(vmec_construct_mpi, intra_comm)
      END IF
#endif

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
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER :: this

!  local variables
      INTEGER                    :: i, error

!  Start of executable code
      IF (ASSOCIATED(this%ne)) THEN
         CALL pprofile_destruct(this%ne)
         this%ne => null()
      END IF

      IF (ASSOCIATED(this%te)) THEN
         CALL pprofile_destruct(this%te)
         this%te => null()
      END IF

      IF (ASSOCIATED(this%sxrem)) THEN
         DO i = 1, SIZE(this%sxrem)
            IF (ASSOCIATED(this%sxrem(i)%p)) THEN
               CALL pprofile_destruct(this%sxrem(i)%p)
               this%sxrem(i)%p => null()
            END IF
         END DO
         DEALLOCATE(this%sxrem)
         this%sxrem => null()
      END IF

      IF (ASSOCIATED(this%magnetic_cache)) THEN
         IF (ASSOCIATED(this%magnetic_cache%rsuv)) THEN
            DEALLOCATE(this%magnetic_cache%rsuv)
            this%magnetic_cache%rsuv => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%zsuv)) THEN
            DEALLOCATE(this%magnetic_cache%zsuv)
            this%magnetic_cache%zsuv => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%jrsuv)) THEN
            DEALLOCATE(this%magnetic_cache%jrsuv)
            this%magnetic_cache%jrsuv => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%jphisuv)) THEN
            DEALLOCATE(this%magnetic_cache%jphisuv)
            this%magnetic_cache%jphisuv => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%jzsuv)) THEN
            DEALLOCATE(this%magnetic_cache%jzsuv)
            this%magnetic_cache%jzsuv => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kruv)) THEN
            DEALLOCATE(this%magnetic_cache%kruv)
            this%magnetic_cache%kruv => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kphiuv)) THEN
            DEALLOCATE(this%magnetic_cache%kphiuv)
            this%magnetic_cache%kphiuv => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kzuv)) THEN
            DEALLOCATE(this%magnetic_cache%kzuv)
            this%magnetic_cache%kzuv => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%x_prime)) THEN
            DEALLOCATE(this%magnetic_cache%x_prime)
            this%magnetic_cache%x_prime => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kxuv_full)) THEN
            DEALLOCATE(this%magnetic_cache%kxuv_full)
            this%magnetic_cache%kxuv_full => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kyuv_full)) THEN
            DEALLOCATE(this%magnetic_cache%kyuv_full)
            this%magnetic_cache%kyuv_full => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kzuv_full)) THEN
            DEALLOCATE(this%magnetic_cache%kzuv_full)
            this%magnetic_cache%kzuv_full => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%x_axi)) THEN
            DEALLOCATE(this%magnetic_cache%x_axi)
            this%magnetic_cache%x_axi => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kxuv_axi)) THEN
            DEALLOCATE(this%magnetic_cache%kxuv_axi)
            this%magnetic_cache%kxuv_axi => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kyuv_axi)) THEN
            DEALLOCATE(this%magnetic_cache%kyuv_axi)
            this%magnetic_cache%kyuv_axi => null()
         END IF

         IF (ASSOCIATED(this%magnetic_cache%kzuv_axi)) THEN
            DEALLOCATE(this%magnetic_cache%kzuv_axi)
            this%magnetic_cache%kzuv_axi => null()
         END IF

         DEALLOCATE(this%magnetic_cache)
         this%magnetic_cache => null()
      END IF

!  Delete the cached wout file. Errors here can safely be ignored.
      CALL delete_file(TRIM(this%wout_file_name) // '_cache', error)

      IF (ASSOCIATED(this%vmec_context_save)) THEN
         CALL vmec_context_destruct(this%vmec_context_save)
         this%vmec_context_save => null()
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
     &                      bloat, mgrid_file, nfp, lfreeb, aphi
      USE xstuff, only: xc, xcdot
      USE vmec_main, only: irzloff, currv, ivac
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
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (id)

         CASE (vmec_rbc_id)
            vmec_set_param = .true.
            rbc(i_index, j_index) = value
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_zbs_id)
            vmec_set_param = .true.
            zbs(i_index, j_index) = value
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_rbs_id)
            vmec_set_param = .true.
            rbs(i_index, j_index) = value
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_zbc_id)
            vmec_set_param = .true.
            zbc(i_index, j_index) = value
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_ac_id)
            vmec_set_param = .true.
            ac(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ac_aux_s_id)
            vmec_set_param = .true.
            ac_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ac_aux_f_id)
            vmec_set_param = .true.
            ac_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ai_id)
            vmec_set_param = .true.
            ai(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ai_aux_s_id)
            vmec_set_param = .true.
            ai_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ai_aux_f_id)
            vmec_set_param = .true.
            ai_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_am_id)
            vmec_set_param = .true.
            am(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_am_aux_s_id)
            vmec_set_param = .true.
            am_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_am_aux_f_id)
            vmec_set_param = .true.
            am_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

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
     &                     'read mgrid_file')

         CASE (vmec_curtor_id)
            vmec_set_param = .true.
            curtor = value
            currv = mu0*curtor
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_phiedge_id)
            vmec_set_param = .true.
            IF (phiedge .ne. 0.0) THEN
               phiedge_ratio = value/phiedge
            ELSE
               phiedge_ratio = 1.0
            END IF
            phiedge = value
            IF (aphi(1) .eq. 1.0) THEN
               xc(1 + 2*irzloff:3*irzloff) =                                   &
     &            phiedge_ratio*xc(1 + 2*irzloff:3*irzloff)
            END IF
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE (vmec_pres_scale_id)
            vmec_set_param = .true.
            pres_scale = value
            CALL profil1d(xc, xcdot, .false.)

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

         CASE (vmec_z_offset_id)
            vmec_set_param = .false.
            this%z_offset = value
            IF (ASSOCIATED(this%magnetic_cache)) THEN
               CALL vmec_set_magnetic_cache(this)
            END IF

         CASE DEFAULT
!>  @todo FIXME: Consider throwing a fatal error here instead of silently doing
!>  nothing.
            vmec_set_param = .true.

      END SELECT

      CALL profiler_set_stop_time('vmec_set_param', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Set magnetic cache for plasma resonce.
!>
!>  This method overrides @ref equilibrium::equilibrium_set_magnetic_cache. This
!>  allocates a @ref vmec_magnetic_cache structure. Point measurements require
!>  no array allocations.
!>
!>  @param[inout] this            A @ref vmec_class instance.
!>  @param[in]    response_object A @ref magnetic_response_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_magnetic_cache_responce(this, response_object)
      USE vmec_input, only: ns_array, mpol
      USE stel_constants, only: twopi
      USE magnetic_response

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout)           :: this
      TYPE (magnetic_response_class), INTENT(in) :: response_object

!  local variables
      INTEGER                                    :: numU
      INTEGER                                    :: ns
      REAL (rprec)                               :: start_time

!  Start of executable code.
      start_time = profiler_get_start_time()

      IF (.not.ASSOCIATED(this%magnetic_cache)) THEN
         ALLOCATE(this%magnetic_cache)
      END IF

      IF (magnetic_response_use_plasma(response_object)) THEN
!  In order to sample the fields on a grid, the sampling frequency of the grid
!  must be at least twice the frequency of the highest mode.
         ns = ns_array(this%ns_index)
         numU = MAX(INT(this%pol_rad_ratio*ns), 2*mpol)

!  Set the differental volume elements.
         this%magnetic_cache%ds = 1.0/(ns - 1.0)
         this%magnetic_cache%du = twopi/numU
         this%magnetic_cache%dv = twopi/response_object%n_field_periods        &
     &                          / response_object%num_t

!  When stellarator symmetry is used, do not need to store all the toroidal
!  planes in a single field period.
         ALLOCATE(this%magnetic_cache%rsuv(ns,numU,                            &
     &            SIZE(response_object%a_r)))
         ALLOCATE(this%magnetic_cache%zsuv(ns,numU,                            &
     &            SIZE(response_object%a_r)))
         ALLOCATE(this%magnetic_cache%jrsuv(ns,numU,                           &
     &            SIZE(response_object%a_r)))
         ALLOCATE(this%magnetic_cache%jphisuv(ns,numU,                         &
     &            SIZE(response_object%a_r)))
         ALLOCATE(this%magnetic_cache%jzsuv(ns,numU,                           &
     &            SIZE(response_object%a_r)))
      END IF

!  When conducting shell is used allocate the surface vector arrays.
      IF (magnetic_response_use_shell(response_object)) THEN
!  Set the differental area elements.
         this%magnetic_cache%du_a = twopi/compression_get_dimension1(          &
     &                                       response_object%a_s_r)
         this%magnetic_cache%dv_a                                              &
     &      = (twopi/response_object%n_field_periods)                          &
     &      / response_object%num_t_shell

         ALLOCATE(this%magnetic_cache%kruv(                                    &
     &               compression_get_dimension1(response_object%a_s_r),        &
     &               compression_get_dimension2(response_object%a_s_r)))
         ALLOCATE(this%magnetic_cache%kphiuv(                                  &
     &               compression_get_dimension1(response_object%a_s_r),        &
     &               compression_get_dimension2(response_object%a_s_r)))
         ALLOCATE(this%magnetic_cache%kzuv(                                    &
     &               compression_get_dimension1(response_object%a_s_r),        &
     &               compression_get_dimension2(response_object%a_s_r)))
      END IF

      CALL profiler_set_stop_time('vmec_set_magnetic_cache_responce',          &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set magnetic cache initialy.
!>
!>  This method overrides @ref equilibrium::equilibrium_set_magnetic_cache. This
!>  allocates a @ref vmec_magnetic_cache structure. Point measurements require
!>  no array allocations.
!>
!>  @param[inout] this    A @ref vmec_class instance.
!>  @param[in]    use_axi Magnetics can subtract off axisymmetric components.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_magnetic_cache_point(this, use_axi)
      USE vmec_input, only: rbc, zbs

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this
      LOGICAL, INTENT(in)              :: use_axi

!  local variables
      INTEGER                          :: u_size, v_size
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (.not.ASSOCIATED(this%magnetic_cache)) THEN
         ALLOCATE(this%magnetic_cache)
      END IF

!  Set the grid size based on the size of the rough size of the plasma.
      v_size = MAX(INT(twopi*rbc(0,0)/magnetic_cache_vc_grid_dim),             &
     &             magnetic_cache_vc_min_grid_points)
      IF (MOD(v_size, 2) .eq. 0) THEN
         v_size = v_size + 1
      END IF
      this%magnetic_cache%dv_full = twopi/v_size

      u_size = MAX(INT(twopi*MAX(rbc(0,1), zbs(0,1))/                          &
     &                 magnetic_cache_vc_grid_dim),                            &
     &                magnetic_cache_vc_min_grid_points)
      IF (MOD(u_size, 2) .eq. 0) THEN
         u_size = u_size + 1
      END IF
      this%magnetic_cache%du_full = twopi/u_size

      ALLOCATE(this%magnetic_cache%kxuv_full(u_size,v_size))
      ALLOCATE(this%magnetic_cache%kyuv_full(u_size,v_size))
      ALLOCATE(this%magnetic_cache%kzuv_full(u_size,v_size))

      ALLOCATE(this%magnetic_cache%x_prime(u_size,v_size,3))

      IF (use_axi) THEN
         ALLOCATE(this%magnetic_cache%kxuv_axi(u_size,v_size))
         ALLOCATE(this%magnetic_cache%kyuv_axi(u_size,v_size))
         ALLOCATE(this%magnetic_cache%kzuv_axi(u_size,v_size))

         ALLOCATE(this%magnetic_cache%x_axi(u_size,v_size,3))

      END IF

      CALL profiler_set_stop_time('vmec_set_magnetic_cache_point',             &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Sync the magnetic cache.
!>
!>  This syncronizes the magnetic cache when spawning external mpi processes.
!>
!>  @param[inout] this       A @ref vmec_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_magnetic_cache_sync(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(in)              :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER, DIMENSION(3)            :: dim_lengths
      INTEGER                          :: error

!  Start of executable code.

!  Grid sizes.
      CALL MPI_BCAST(this%magnetic_cache%ds, 1, MPI_REAL8, 0,                  &
     &               intra_comm, error)
      CALL MPI_BCAST(this%magnetic_cache%du, 1, MPI_REAL8, 0,                  &
     &               intra_comm, error)
      CALL MPI_BCAST(this%magnetic_cache%dv, 1, MPI_REAL8, 0,                  &
     &               intra_comm, error)

      CALL MPI_BCAST(this%magnetic_cache%du_a, 1, MPI_REAL8, 0,                &
     &               intra_comm, error)
      CALL MPI_BCAST(this%magnetic_cache%dv_a, 1, MPI_REAL8, 0,                &
     &               intra_comm, error)

      CALL MPI_BCAST(this%magnetic_cache%du_full, 1, MPI_REAL8, 0,             &
     &               intra_comm, error)
      CALL MPI_BCAST(this%magnetic_cache%dv_full, 1, MPI_REAL8, 0,             &
     &               intra_comm, error)

!  Plasma resonse arrays.
      dim_lengths = 0
      IF (ASSOCIATED(this%magnetic_cache%rsuv)) THEN
         dim_lengths = SHAPE(this%magnetic_cache%rsuv)
      END IF
      CALL MPI_BCAST(dim_lengths, 3, MPI_INTEGER, 0, intra_comm, error)
      IF (PRODUCT(dim_lengths) .gt. 0) THEN
         IF (.not.ASSOCIATED(this%magnetic_cache%rsuv)) THEN
            ALLOCATE(this%magnetic_cache%rsuv(dim_lengths(1),                  &
     &                                        dim_lengths(2),                  &
     &                                        dim_lengths(3)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%zsuv)) THEN
            ALLOCATE(this%magnetic_cache%zsuv(dim_lengths(1),                  &
     &                                        dim_lengths(2),                  &
     &                                        dim_lengths(3)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%jrsuv)) THEN
            ALLOCATE(this%magnetic_cache%jrsuv(dim_lengths(1),                 &
     &                                         dim_lengths(2),                 &
     &                                         dim_lengths(3)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%jphisuv)) THEN
            ALLOCATE(this%magnetic_cache%jphisuv(dim_lengths(1),               &
     &                                           dim_lengths(2),               &
     &                                           dim_lengths(3)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%jzsuv)) THEN
            ALLOCATE(this%magnetic_cache%jzsuv(dim_lengths(1),                 &
     &                                         dim_lengths(2),                 &
     &                                         dim_lengths(3)))
         END IF
      END IF

!  Conducting shell arrays.
      dim_lengths = 0
      IF (ASSOCIATED(this%magnetic_cache%kruv)) THEN
         dim_lengths(1:2) = SHAPE(this%magnetic_cache%kruv)
      END IF
      CALL MPI_BCAST(dim_lengths(1:2), 2, MPI_INTEGER, 0, intra_comm,          &
     &               error)
      IF (PRODUCT(dim_lengths(1:2)) .gt. 0) THEN
         IF (.not.ASSOCIATED(this%magnetic_cache%kruv)) THEN
            ALLOCATE(this%magnetic_cache%kruv(dim_lengths(1),                  &
     &                                        dim_lengths(2)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%kphiuv)) THEN
            ALLOCATE(this%magnetic_cache%kphiuv(dim_lengths(1),                &
     &                                          dim_lengths(2)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%kzuv)) THEN
            ALLOCATE(this%magnetic_cache%kzuv(dim_lengths(1),                  &
     &                                        dim_lengths(2)))
         END IF
      END IF

!  Point virtual casing arrays.
      dim_lengths = 0
      IF (ASSOCIATED(this%magnetic_cache%kxuv_full)) THEN
         dim_lengths(1:2) = SHAPE(this%magnetic_cache%kxuv_full)
      END IF
      CALL MPI_BCAST(dim_lengths(1:2), 2, MPI_INTEGER, 0, intra_comm,          &
     &               error)
      IF (PRODUCT(dim_lengths(1:2)) .gt. 0) THEN
         IF (.not.ASSOCIATED(this%magnetic_cache%kxuv_full)) THEN
            ALLOCATE(this%magnetic_cache%kxuv_full(dim_lengths(1),             &
     &                                           dim_lengths(2)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%kyuv_full)) THEN
            ALLOCATE(this%magnetic_cache%kyuv_full(dim_lengths(1),             &
     &                                             dim_lengths(2)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%kzuv_full)) THEN
            ALLOCATE(this%magnetic_cache%kzuv_full(dim_lengths(1),             &
     &                                             dim_lengths(2)))
         END IF
      END IF

      dim_lengths = 0
      IF (ASSOCIATED(this%magnetic_cache%x_prime)) THEN
         dim_lengths = SHAPE(this%magnetic_cache%x_prime)
      END IF
      CALL MPI_BCAST(dim_lengths, 3, MPI_INTEGER, 0, intra_comm, error)
      IF (PRODUCT(dim_lengths) .gt. 0 .and.                                    &
     &    .not.ASSOCIATED(this%magnetic_cache%x_prime)) THEN
         ALLOCATE(this%magnetic_cache%x_prime(dim_lengths(1),                  &
     &                                        dim_lengths(2),                  &
     &                                        dim_lengths(3)))
      END IF

!  Point axisymmtric virtual casing arrays.
      dim_lengths = 0
      IF (ASSOCIATED(this%magnetic_cache%kxuv_axi)) THEN
         dim_lengths(1:2) = SHAPE(this%magnetic_cache%kxuv_axi)
      END IF
      CALL MPI_BCAST(dim_lengths(1:2), 2, MPI_INTEGER, 0, intra_comm,          &
     &               error)
      IF (PRODUCT(dim_lengths(1:2)) .gt. 0) THEN
         IF (.not.ASSOCIATED(this%magnetic_cache%kxuv_axi)) THEN
            ALLOCATE(this%magnetic_cache%kxuv_axi(dim_lengths(1),              &
     &                                            dim_lengths(2)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%kyuv_axi)) THEN
            ALLOCATE(this%magnetic_cache%kyuv_axi(dim_lengths(1),              &
     &                                            dim_lengths(2)))
         END IF
         IF (.not.ASSOCIATED(this%magnetic_cache%kzuv_axi)) THEN
            ALLOCATE(this%magnetic_cache%kzuv_axi(dim_lengths(1),              &
     &                                            dim_lengths(2)))
         END IF
      END IF

      dim_lengths = 0
      IF (ASSOCIATED(this%magnetic_cache%x_axi)) THEN
         dim_lengths = SHAPE(this%magnetic_cache%x_axi)
      END IF
      CALL MPI_BCAST(dim_lengths, 3, MPI_INTEGER, 0, intra_comm, error)
      IF (PRODUCT(dim_lengths) .gt. 0 .and.                                    &
     &    .not.ASSOCIATED(this%magnetic_cache%x_axi)) THEN
         ALLOCATE(this%magnetic_cache%x_axi(dim_lengths(1),                    &
     &                                      dim_lengths(2),                    &
     &                                      dim_lengths(3)))
      END IF

#endif

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
     &                         currumnc, currvmnc, currumns, currvmns,         &
     &                         bsubumnc, bsubumns, bsubvmnc, bsubvmns,         &
     &                         nfp, ns, lwout_opened, lfreeb,                  &
     &                         bsubumnc_sur, bsubumns_sur,                     &
     &                         bsubvmnc_sur, bsubvmns_sur
      USE stel_constants, only: mu0
      USE virtual_casing_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout)          :: this

!  local variables
      INTEGER                                   :: s, u, v, i
      INTEGER                                   :: numU, numV
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosz, sinz
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosu, sinu
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: cosv, sinv
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosphi, sinphi
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosz_nyq, sinz_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosu_nyq, sinu_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: cosv_nyq, sinv_nyq
      REAL (rprec)                              :: ru, rv, zu, zv, r, z
      REAL (rprec)                              :: ju, jv, bu, bv
      REAL (rprec)                              :: kr, kphi
      REAL (rprec)                              :: phi
      REAL (rprec)                              :: start_time
      REAl (rprec), DIMENSION(:), ALLOCATABLE   :: bsubuc, bsubus
      REAl (rprec), DIMENSION(:), ALLOCATABLE   :: bsubvc, bsubvs
      REAL (rprec)                              :: woddlow, woddhigh

!  Start of executable code.
      start_time = profiler_get_start_time()

!  Only generate the magnetic caches if the woutfile is has been opened. This
!  prevents a segfaults when the child processes sync the phi_edge and z_offset
!  parameters before the first reconstruction step when no wout file exists in
!  the child process working directory.
      IF (.not.lwout_opened) THEN
         CALL profiler_set_stop_time('vmec_set_magnetic_cache_calc',           &
     &                               start_time)
         RETURN
      END IF

      IF (ASSOCIATED(this%magnetic_cache%rsuv)) THEN
!  Get the size of the SUV grid.
         numU = SIZE(this%magnetic_cache%rsuv, 2)
         numV = SIZE(this%magnetic_cache%rsuv, 3)

         ALLOCATE(cosv(mnmax, numV), sinv(mnmax, numV))
         ALLOCATE(cosv_nyq(mnmax_nyq, numV), sinv_nyq(mnmax_nyq, numV))

!  U grid is square with the s grid.
!>  @todo FIXME: There is no reason to recalculate the sin and cosin for the u
!>  component arrays every time. The v components can change if phi_offset is
!>  being reconstructed. The u components can be cached in the
!>  @ref vmec_magnetic_cache and calculated once in
!>  @ref vmec_set_magnetic_cache_init.
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(s,u,v,cosz,sinz,cosz_nyq,sinz_nyq,ru,zu,rv,zv,ju,jv,
!$OMP&         cosu,sinu,cosu_nyq,sinu_nyq,phi)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO v = 1, numV
            phi = (v - 1)*this%magnetic_cache%dv + this%phi_offset
            cosv(:,v) = COS(xn*phi)
            sinv(:,v) = SIN(xn*phi)
            cosv_nyq(:,v) = COS(xn_nyq*phi)
            sinv_nyq(:,v) = SIN(xn_nyq*phi)
         END DO
!$OMP END DO

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

!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO u = 1, numU
            ALLOCATE(cosu(mnmax), sinu(mnmax))
            ALLOCATE(cosz(mnmax), sinz(mnmax))
            ALLOCATE(cosu_nyq(mnmax_nyq), sinu_nyq(mnmax_nyq))
            ALLOCATE(cosz_nyq(mnmax_nyq), sinz_nyq(mnmax_nyq))

            cosu = COS(xm*(u - 1)*this%magnetic_cache%du)
            sinu = SIN(xm*(u - 1)*this%magnetic_cache%du)
            cosu_nyq = COS(xm_nyq*(u - 1)*this%magnetic_cache%du)
            sinu_nyq = SIN(xm_nyq*(u - 1)*this%magnetic_cache%du)

            DO v = 1, numV

               cosz = cosu*cosv(:,v) + sinu*sinv(:,v)
               sinz = sinu*cosv(:,v) - cosu*sinv(:,v)

               cosz_nyq = cosu_nyq*cosv_nyq(:,v)                               &
     &                  + sinu_nyq*sinv_nyq(:,v)
               sinz_nyq = sinu_nyq*cosv_nyq(:,v)                               &
     &                  - cosu_nyq*sinv_nyq(:,v)

               DO s = 1, ns

                  this%magnetic_cache%rsuv(s,u,v) = SUM(rmnc(:,s)*cosz)
                  this%magnetic_cache%zsuv(s,u,v) = SUM(zmns(:,s)*sinz)
                  ru = -SUM(xm*rmnc(:,s)*sinz)
                  zu = SUM(xm*zmns(:,s)*cosz)
                  rv = SUM(xn*rmnc(:,s)*sinz)
                  zv = -SUM(xn*zmns(:,s)*cosz)

                  ju = SUM(currumnc(:,s)*cosz_nyq)
                  jv = SUM(currvmnc(:,s)*cosz_nyq)

                  IF (lasym) THEN
                     this%magnetic_cache%rsuv(s,u,v) =                         &
     &                  this%magnetic_cache%rsuv(s,u,v) +                      &
     &                  SUM(rmns(:,s)*sinz)
                     this%magnetic_cache%zsuv(s,u,v) =                         &
     &                  this%magnetic_cache%zsuv(s,u,v) +                      &
     &                  SUM(zmnc(:,s)*cosz)
                     ru = ru + SUM(xm*rmns(:,s)*cosz)
                     zu = zu - SUM(xm*zmnc(:,s)*sinz)
                     rv = rv - SUM(xn*rmns(:,s)*cosz)
                     zv = zv + SUM(xn*zmnc(:,s)*sinz)

                     ju = ju + SUM(currumns(:,s)*sinz_nyq)
                     jv = jv + SUM(currvmns(:,s)*sinz_nyq)
                  END IF

                  this%magnetic_cache%jrsuv(s,u,v) = ju*ru + jv*rv
                  this%magnetic_cache%jphisuv(s,u,v) =                            &
     &               jv*this%magnetic_cache%rsuv(s,u,v)
                  this%magnetic_cache%jzsuv(s,u,v) = ju*zu + jv*zv
               END DO
            END DO

            DEALLOCATE(cosu, sinu)
            DEALLOCATE(cosz, sinz)
            DEALLOCATE(cosu_nyq, sinu_nyq)
            DEALLOCATE(cosz_nyq, sinz_nyq)
         END DO
!$OMP END DO
!$OMP END PARALLEL

         this%magnetic_cache%zsuv = this%magnetic_cache%zsuv                   &
     &                            + this%z_offset

         DEALLOCATE(cosv, sinv)
         DEALLOCATE(cosv_nyq, sinv_nyq)
      END IF

!  Compute the conducting shell if needed. If the kruv array is associated, the
!  the conducting shell is being used. The total conducting shell signal for a
!  magnetic diagnostic is
!
!    S_i = AreaIntegral K(u,v) . R_i(u,v) da                                 (1)
!
!  For now, It will be assumed that the last closed flux surface of the plasma
!  is the conducting shell. This can be refined later. The total magnetic signal
!  in flux coordinates becomes
!
!    S_i = AreaIntegral (-e^s/|e^s| X B/mu0) . R_i * |J||e^s|dudv            (2)
!
!  Using the idenity e^s = 1/J*e_u X e_v equation 2 becomes
!
!    S_i = AreaIntegral -e_u X e_v X B/mu0 . R_i |J|/J dudv                  (3)
!
!  From the triple cross product idenity equation 3 reduces to
!
!    S_i = AreaIntegral -1/mu0(B_u e_v - B_v e_u) . R_i |J|/J dudv           (4)
!
!  Computing the conducting shell surface current is simply.
!
!    K(u,v) = -1/mu0(B_u*e_v - B_v*e_u)                                      (5)
!
!  NOTE: The shell current cannot be used with the Z position shift.
!
!  NOTE: Equation 2.5.49a in "Flux Coordinates and Magnetic Field Structure" by
!        W.D.D'haeseleer, W.N.G.Hitchon, J.D.Callen and J.L.Shohet assumes a
!        positive jacobian. The correct equation should be
!
!            dS(i)=|J||e^i| du^j dv^k
!
!        As a result the signal due to the conducting shell contains a |J|/J
!        term. The sign of the jacobian needs to be taken into account. This is
!        taken care of in vmec_get_area_int_element.
!
!  Compute virtual casing points if needed, If the kruv_full array is
!  associated, the point magnetics are being used. The total signal for a point
!  is
!
!    B(x) = mu0/2Pi AreaIntegral K(u,v) x (x - x')/|x - x'|^3 da             (6)
!
!  The surface used is the last closed flux surface. Simular to the procedure
!  for the conducting shell, the different area element and surface current in
!  equation 2 reduce the point field to.
!
!    B(x) = mu0/2Pi AreaIntegral -1/mu0(B_u*e_v - B_v*e_u)
!         x (x - x')/|x - x'|^3 |J|/J dudv                                   (7)
!
!  Computing the virtual casing surface current is simply.
!
!    K(u,v) = -1/mu0(B_u*e_v - B_v*e_u)
!
!  NOTE: Unlike the shell current, the virtual casing can be used with the Z
!        position shift.

      IF (ASSOCIATED(this%magnetic_cache%kruv) .or.                            &
     &    ASSOCIATED(this%magnetic_cache%kxuv_full)) THEN

         IF (lfreeb) THEN
            bsubuc = bsubumnc_sur
            bsubvc = bsubvmnc_sur
            IF (lasym) THEN
               bsubus = bsubumns_sur
               bsubvs = bsubvmns_sur
            END IF
         ELSE
!  Interpolate the b_u and b_v onto the last closed flux surface. This
!  interpolation follows what is outlined in tosuvspace. For odd modes
!  X ~ SQRT(s) so interpolate Xmn/SQRT(s)
            woddlow = 0.5*SQRT((ns - 1.0)/(ns - 1.5))
            woddhigh = 1.5*SQRT((ns - 1.0)/(ns - 2.5))

            ALLOCATE(bsubuc(mnmax_nyq), bsubvc(mnmax_nyq))
            WHERE (MOD(NINT(xm_nyq), 2) .eq. 0)
               bsubuc = 1.5*bsubumnc(:,ns) - 0.5*bsubumnc(:,ns - 1)
               bsubvc = 1.5*bsubvmnc(:,ns) - 0.5*bsubvmnc(:,ns - 1)
            ELSE WHERE
               bsubuc = woddhigh*bsubumnc(:,ns)                                &
     &                - woddlow*bsubumnc(:,ns - 1)
               bsubvc = woddhigh*bsubvmnc(:,ns)                                &
     &                - woddlow*bsubvmnc(:,ns - 1)
            END WHERE

            IF (lasym) THEN
               ALLOCATE(bsubus(mnmax_nyq), bsubvs(mnmax_nyq))
               WHERE (MOD(NINT(xm_nyq), 2) .eq. 0)
                  bsubus = 1.5*bsubumns(:,ns) - 0.5*bsubumns(:,ns - 1)
                  bsubvs = 1.5*bsubvmns(:,ns) - 0.5*bsubvmns(:,ns - 1)
               ELSE WHERE
                  bsubus = woddhigh*bsubumns(:,ns)                             &
     &                   - woddlow*bsubumns(:,ns - 1)
                  bsubvs = woddhigh*bsubvmns(:,ns)                             &
     &                   - woddlow*bsubvmns(:,ns - 1)
               END WHERE
            END IF
         END IF
      END IF

      IF (ASSOCIATED(this%magnetic_cache%kruv)) THEN
         numU = SIZE(this%magnetic_cache%kruv, 1)
         numV = SIZE(this%magnetic_cache%kruv, 2)

         ALLOCATE(cosv(mnmax, numV), sinv(mnmax, numV))
         ALLOCATE(cosv_nyq(mnmax_nyq, numV), sinv_nyq(mnmax_nyq, numV))

!  U grid is square with the s grid.
!>  @todo FIXME: There is no reason to recalculate the sin and cosin for the u
!>  component arrays every time. The v components can change if phi_offset is
!>  being reconstructed. The u components can be cached in the
!>  @ref vmec_magnetic_cache and calculated once in
!>  @ref vmec_set_magnetic_cache_init.
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(u,v,cosz,sinz,cosz_nyq,sinz_nyq,ru,zu,rv,zv,bu,bv,r,
!$OMP&         cosu,sinu,cosu_nyq,sinu_nyq,phi)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO v = 1, numV
            phi = (v - 1)*this%magnetic_cache%dv_a + this%phi_offset
            cosv(:,v) = COS(xn*phi)
            sinv(:,v) = SIN(xn*phi)
            cosv_nyq(:,v) = COS(xn_nyq*phi)
            sinv_nyq(:,v) = SIN(xn_nyq*phi)
         END DO
!$OMP END DO

!  U grid is square with the v grid.
!>  @todo FIXME: There is no reason to recompute the sine and cosine for the u
!>  component arrays every time.
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO u = 1, numU
!  Arrays need to be allocated within the loop since each thread requires its
!  own allocation.
            ALLOCATE(cosu(mnmax), sinu(mnmax))
            ALLOCATE(cosz(mnmax), sinz(mnmax))
            ALLOCATE(cosu_nyq(mnmax_nyq), sinu_nyq(mnmax_nyq))
            ALLOCATE(cosz_nyq(mnmax_nyq), sinz_nyq(mnmax_nyq))

            cosu = COS(xm*(u - 1)*this%magnetic_cache%du_a)
            sinu = SIN(xm*(u - 1)*this%magnetic_cache%du_a)
            cosu_nyq = COS(xm_nyq*(u - 1)*this%magnetic_cache%du_a)
            sinu_nyq = SIN(xm_nyq*(u - 1)*this%magnetic_cache%du_a)

            DO v = 1, numV

               cosz = cosu*cosv(:,v) + sinu*sinv(:,v)
               sinz = sinu*cosv(:,v) - cosu*sinv(:,v)

               cosz_nyq = cosu_nyq*cosv_nyq(:,v)                               &
     &                  + sinu_nyq*sinv_nyq(:,v)

               ru = -SUM(xm*rmnc(:,ns)*sinz)
               zu = SUM(xm*zmns(:,ns)*cosz)
               rv = SUM(xn*rmnc(:,ns)*sinz)
               zv = -SUM(xn*zmns(:,ns)*cosz)
               r = SUM(rmnc(:,ns)*cosz)

               bu = SUM(bsubuc*cosz_nyq)
               bv = SUM(bsubvc*cosz_nyq)

               IF (lasym) THEN
                  sinz_nyq = sinu_nyq*cosv_nyq(:,v)                            &
     &                     - cosu_nyq*sinv_nyq(:,v)

                  ru = ru + SUM(xm*rmns(:,ns)*cosz)
                  zu = zu - SUM(xm*zmnc(:,ns)*sinz)
                  rv = rv - SUM(xn*rmns(:,ns)*cosz)
                  zv = zv + SUM(xn*zmnc(:,ns)*sinz)
                  r = r + SUM(rmns(:,ns)*sinz)

                  bu = bu + SUM(bsubus*sinz_nyq)
                  bv = bv + SUM(bsubvs*sinz_nyq)
               END IF

               this%magnetic_cache%kruv(u,v) = -1.0/mu0*(bu*rv - bv*ru)
               this%magnetic_cache%kphiuv(u,v) = -1.0/mu0*bu*r
               this%magnetic_cache%kzuv(u,v) = -1.0/mu0*(bu*zv - bv*zu)
            END DO

            DEALLOCATE(cosu, sinu)
            DEALLOCATE(cosz, sinz)
            DEALLOCATE(cosu_nyq, sinu_nyq)
            DEALLOCATE(cosz_nyq, sinz_nyq)
         END DO
!$OMP END DO
!$OMP END PARALLEL

         DEALLOCATE(cosv, sinv)
         DEALLOCATE(cosv_nyq, sinv_nyq)
      END IF

      IF (ASSOCIATED(this%magnetic_cache%kxuv_full)) THEN
         numU = SIZE(this%magnetic_cache%kxuv_full, 1)
         numV = SIZE(this%magnetic_cache%kxuv_full, 2)

         ALLOCATE(cosv(mnmax, numV), sinv(mnmax, numV))
         ALLOCATE(cosphi(numV), sinphi(numV))
         ALLOCATE(cosv_nyq(mnmax_nyq, numV), sinv_nyq(mnmax_nyq, numV))

!>  @todo FIXME: There is no reason to recalculate the sin and cosin for the u
!>  component arrays every time. The v components can change if phi_offset is
!>  being reconstructed. The u components can be cached in the
!>  @ref vmec_magnetic_cache and calculated once in
!>  @ref vmec_set_magnetic_cache_init.
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(u,v,cosz,sinz,cosz_nyq,sinz_nyq,ru,zu,rv,zv,bu,bv,r, z,
!$OMP&         cosu,sinu,cosu_nyq,sinu_nyq,phi,kr,kphi,i)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO v = 1, numV
!  Rotate the signal with respect to a fixed plasma. The phi and z offsets are
!  applied in vmec_get_ext_b_plasma.
            phi = (v - 1)*this%magnetic_cache%dv_full
            cosphi(v) = COS(phi)
            sinphi(v) = SIN(phi)
            cosv(:,v) = COS(xn*phi)
            sinv(:,v) = SIN(xn*phi)
            cosv_nyq(:,v) = COS(xn_nyq*phi)
            sinv_nyq(:,v) = SIN(xn_nyq*phi)
         END DO
!$OMP END DO

!>  @todo FIXME: There is no reason to recompute the sine and cosine for the u
!>  component arrays every time.
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO u = 1, numU
!  Arrays need to be allocated within the loop since each thread requires its
!  own allocation.
            ALLOCATE(cosu(mnmax), sinu(mnmax))
            ALLOCATE(cosz(mnmax), sinz(mnmax))
            ALLOCATE(cosu_nyq(mnmax_nyq), sinu_nyq(mnmax_nyq))
            ALLOCATE(cosz_nyq(mnmax_nyq), sinz_nyq(mnmax_nyq))

            cosu = COS(xm*(u - 1)*this%magnetic_cache%du_full)
            sinu = SIN(xm*(u - 1)*this%magnetic_cache%du_full)
            cosu_nyq = COS(xm_nyq*(u - 1)*this%magnetic_cache%du_full)
            sinu_nyq = SIN(xm_nyq*(u - 1)*this%magnetic_cache%du_full)

            DO v = 1, numV

               cosz = cosu*cosv(:,v) + sinu*sinv(:,v)
               sinz = sinu*cosv(:,v) - cosu*sinv(:,v)

               cosz_nyq = cosu_nyq*cosv_nyq(:,v)                               &
     &                  + sinu_nyq*sinv_nyq(:,v)

               ru = -SUM(xm*rmnc(:,ns)*sinz)
               zu = SUM(xm*zmns(:,ns)*cosz)
               rv = SUM(xn*rmnc(:,ns)*sinz)
               zv = -SUM(xn*zmns(:,ns)*cosz)
               r = SUM(rmnc(:,ns)*cosz)
               z = SUM(zmns(:,ns)*sinz)

               bu = SUM(bsubuc*cosz_nyq)
               bv = SUM(bsubvc*cosz_nyq)

               IF (lasym) THEN
                  sinz_nyq = sinu_nyq*cosv_nyq(:,v)                            &
     &                     - cosu_nyq*sinv_nyq(:,v)

                  ru = ru + SUM(xm*rmns(:,ns)*cosz)
                  zu = zu - SUM(xm*zmnc(:,ns)*sinz)
                  rv = rv - SUM(xn*rmns(:,ns)*cosz)
                  zv = zv + SUM(xn*zmnc(:,ns)*sinz)
                  r = r + SUM(rmns(:,ns)*sinz)
                  z = z + SUM(zmnc(:,ns)*cosz)

                  bu = bu + SUM(bsubus*sinz_nyq)
                  bv = bv + SUM(bsubvs*sinz_nyq)
               END IF

               kr = -1.0/mu0*(bu*rv - bv*ru)
               kphi = -1.0/mu0*bu*r

               this%magnetic_cache%kxuv_full(u,v) = kr*cosphi(v)               &
     &                                            - kphi*sinphi(v)
               this%magnetic_cache%kyuv_full(u,v) = kr*sinphi(v)               &
     &                                            + kphi*cosphi(v)
               this%magnetic_cache%kzuv_full(u,v) = -1.0/mu0*(bu*zv -          &
     &                                                        bv*zu)

               this%magnetic_cache%x_prime(u,v,1) = r*cosphi(v)
               this%magnetic_cache%x_prime(u,v,2) = r*sinphi(v)
               this%magnetic_cache%x_prime(u,v,3) = z

               IF (ASSOCIATED(this%magnetic_cache%kxuv_axi)) THEN
                  r = 0.0
                  z = 0.0
                  ru = 0.0
                  zu = 0.0

                  bu = 0.0
                  bv = 0.0

                  DO i = 1, SIZE(xn)
                     IF (xn(i) .eq. 0.0) THEN
                        ru = ru - xm(i)*rmnc(i,ns)*sinz(i)
                        zu = zu + xm(i)*zmns(i,ns)*cosz(i)
                        r = r + rmnc(i,ns)*cosz(i)
                        z = z + zmns(i,ns)*sinz(i)
                     END IF
                  END DO

                  DO i = 1, SIZE(xn_nyq)
                     IF (xn_nyq(i) .eq. 0.0) THEN
                        bu = bu + bsubuc(i)*cosz_nyq(i)
                        bv = bv + bsubvc(i)*cosz_nyq(i)
                     END IF
                  END DO

                  IF (lasym) THEN
                     DO i = 1, SIZE(xn)
                        IF (xn(i) .eq. 0.0) THEN
                           ru = ru + xm(i)*rmns(i,ns)*cosz(i)
                           zu = zu - xm(i)*zmnc(i,ns)*sinz(i)
                           r = r + rmns(i,ns)*sinz(i)
                           z = z + zmnc(i,ns)*cosz(i)
                        END IF
                     END DO

                     DO i = 1, SIZE(xn_nyq)
                        IF (xn_nyq(i) .eq. 0.0) THEN
                           bu = bu + bsubus(i)*sinz_nyq(i)
                           bv = bv + bsubvs(i)*sinz_nyq(i)
                        END IF
                     END DO
                  END IF

                  kr = 1.0/mu0*bv*ru
                  kphi = -1.0/mu0*bu*r

                  this%magnetic_cache%kxuv_axi(u,v) = kr*cosphi(v)             &
     &                                              - kphi*sinphi(v)
                  this%magnetic_cache%kyuv_axi(u,v) = kr*sinphi(v)             &
     &                                              + kphi*cosphi(v)
                  this%magnetic_cache%kzuv_axi(u,v) = 1.0/mu0*bv*zu

                  this%magnetic_cache%x_axi(u,v,1) = r*cosphi(v)
                  this%magnetic_cache%x_axi(u,v,2) = r*sinphi(v)
                  this%magnetic_cache%x_axi(u,v,3) = z

               END IF
            END DO

            DEALLOCATE(cosu, sinu)
            DEALLOCATE(cosz, sinz)
            DEALLOCATE(cosu_nyq, sinu_nyq)
            DEALLOCATE(cosz_nyq, sinz_nyq)
         END DO
!$OMP END DO
!$OMP END PARALLEL

         DEALLOCATE(cosv, sinv)
         DEALLOCATE(cosphi, sinphi)
         DEALLOCATE(cosv_nyq, sinv_nyq)
      END IF

      IF (ASSOCIATED(this%magnetic_cache%kruv) .or.                            &
     &    ASSOCIATED(this%magnetic_cache%kzuv_full)) THEN
         DEALLOCATE(bsubuc, bsubvc)

         IF (lasym) THEN
            DEALLOCATE(bsubus, bsubvs)
         END IF
      END IF

!  Setup the virtual casing arrays.

      CALL profiler_set_stop_time('vmec_set_magnetic_cache_calc',              &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Sets namelist variables from the solved VMEC equilibrium.
!>
!>  This method updates the boundary and magnetic axis coefficients in the
!>  namelist input file.
!>
!>  @param[in] this A @ref vmec_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_namelist(this)
      USE read_wout_mod, only: rmnc, rmns, zmnc, zmns, xm, xn, mnmax
      USE vmec_input

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this

!  local variables
      INTEGER                       :: i
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Write out the final namelist input file.
!  Update the boundary from the free boundary solution.
      IF (lfreeb) THEN
         DO i = 1, mnmax
            rbc(INT(xn(i)),INT(xm(i))) = rmnc(i,ns_array(this%ns_index))
            zbs(INT(xn(i)),INT(xm(i))) = zmns(i,ns_array(this%ns_index))
            IF (lasym) THEN
               rbs(INT(xn(i)),INT(xm(i))) =                                    &
     &            rmns(i,ns_array(this%ns_index))
               zbc(INT(xn(i)),INT(xm(i))) =                                    &
     &            zmnc(i,ns_array(this%ns_index))
            END IF
         END DO
      END IF

!  Update the magnetic axis.
      raxis_cc(0:ntor) = rmnc(1,1:ntor + 1)
      zaxis_cs(0:ntor) = zmns(1,1:ntor + 1)
      IF (lasym) THEN
         raxis_cs(0:ntor) = rmns(1,1:ntor + 1)
         zaxis_cc(0:ntor) = zmnc(1,1:ntor + 1)
      END IF

      CALL profiler_set_stop_time('vmec_set_namelist', start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

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

         CASE ('z_offset')
            vmec_get_param_id = vmec_z_offset_id

         CASE DEFAULT
            WRITE (*,1000) TRIM(param_name)
            CALL EXIT(1)

      END SELECT

      CALL profiler_set_stop_time('vmec_get_param_id', start_time)

1000  FORMAT('ERROR: ',a,' is not a valid parameter.')

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

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
            vmec_get_param_value = rbs(i_index, j_index)

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

         CASE (vmec_z_offset_id)
            vmec_get_param_value = this%z_offset

         CASE DEFAULT
            vmec_get_param_value = 0.0

      END SELECT

      CALL profiler_set_stop_time('vmec_get_param_value', start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

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

         CASE (vmec_z_offset_id)
            vmec_get_param_name = 'z_offset'

         CASE DEFAULT
            STOP 'Invalid VMEC parameter id.'

      END SELECT

      CALL profiler_set_stop_time('vmec_get_param_name', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron density at a cartesian position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_ne. If no density
!>  profile was created, return zero.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the electron density at.
!>  @returns The electron density at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_ne_cart(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_ne_cart
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_ne_cart = vmec_get_ne(this, vmec_get_s(this, x_cart))

      CALL profiler_set_stop_time('vmec_get_ne_cart', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron density at a s position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_ne. If no density
!>  profile was created, return zero.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] s    Radial position to get the electron density at.
!>  @returns The electron density at s.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_ne_radial(this, s)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: vmec_get_ne_radial
      TYPE (vmec_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)      :: s

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (ASSOCIATED(this%ne)) THEN
         vmec_get_ne_radial = pprofile_get_value(this%ne, s)
      ELSE
         vmec_get_ne_radial = 0.0
      END IF

      CALL profiler_set_stop_time('vmec_get_ne_radial', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron temperature at a cartesian position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_te. If no
!>  temperature profile was created, return zero.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the electron temperature at.
!>  @returns The electron temperature at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_te_cart(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_te_cart
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_te_cart = vmec_get_te(this, vmec_get_s(this, x_cart))

      CALL profiler_set_stop_time('vmec_get_te_cart', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron temperature at a s position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_te. If no
!>  temperature profile was created, return zero.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] s    Radial position to get the electron temperature at.
!>  @returns The electron temperature at s.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_te_radial(this, s)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: vmec_get_te_radial
      TYPE (vmec_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)      :: s

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      IF (ASSOCIATED(this%te)) THEN
         vmec_get_te_radial = pprofile_get_value(this%te, s)
      ELSE
         vmec_get_te_radial = 0.0
      END IF

      CALL profiler_set_stop_time('vmec_get_te_radial', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the soft x-ray emissivity at a cartesian position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_sxrem. If no soft
!>  x-ray emissivity profile was created, return zero.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the soft x-ray emissivity at.
!>  @param[in] index  Index of the soft x-ray emissivity profile.
!>  @returns The soft x-ray emissivity at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_sxrem_cart(this, x_cart, index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_sxrem_cart
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      INTEGER, INTENT(in)                    :: index

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_sxrem_cart = vmec_get_sxrem(this,                               &
     &                                     vmec_get_s(this, x_cart),           &
     &                                     index)

      CALL profiler_set_stop_time('vmec_get_sxrem_cart', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the soft x-ray emissivity at a s position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_sxrem. If no soft
!>  x-ray emissivity profile was created, return zero.
!>
!>  @param[in] this  A @ref vmec_class instance.
!>  @param[in] s     Radial position to get the soft x-ray emissivity at.
!>  @param[in] index Index of the soft x-ray emissivity profile.
!>  @returns The soft x-ray emissivity at s.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_sxrem_radial(this, s, index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: vmec_get_sxrem_radial
      TYPE (vmec_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)      :: s
      INTEGER, INTENT(in)           :: index

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (ASSOCIATED(this%sxrem)) THEN
         vmec_get_sxrem_radial =                                               &
     &      pprofile_get_value(this%sxrem(index)%p, s)
      ELSE
         vmec_get_sxrem_radial = 0.0
      END IF

      CALL profiler_set_stop_time('vmec_get_sxrem_radial', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the plasma pressure at a cartesian position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_p.
!>
!>  @param[in] this   A @ref vmec_class instance.
!>  @param[in] x_cart Cartesian position to get the plasma pressure at.
!>  @returns The plasma pressure at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_p_cart(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_p_cart
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_p_cart = vmec_get_p(this, vmec_get_s(this, x_cart))

      CALL profiler_set_stop_time('vmec_get_p_cart', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the plasma pressure at a s position.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_p.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @param[in] s    Radial position to get the plasma pressure at.
!>  @returns The plasma pressure at s.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_p_radial(this, s)
      USE stel_constants, only: mu0

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: vmec_get_p_radial
      TYPE (vmec_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)      :: s

!  local variables
      REAL (rprec)                  :: start_time

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
      start_time = profiler_get_start_time()
      vmec_get_p_radial = pmass(s)/mu0

      CALL profiler_set_stop_time('vmec_get_p_radial', start_time)

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
      USE mgrid_mod
      USE vmec_input, ONLY: nfp

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: vmec_get_B_vec
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  local variables
      REAL (rprec), DIMENSION(3)             :: r_cyl
      REAL (rprec)                           :: v_norm
      INTEGER                                :: i_low
      INTEGER                                :: j_low
      INTEGER                                :: k_low
      REAL (rprec)                           :: r_low
      REAL (rprec)                           :: phi_low
      REAL (rprec)                           :: z_low
      REAL (rprec)                           :: delphib
      REAL (rprec)                           :: wrl, wrh
      REAL (rprec)                           :: wphil, wphih
      REAL (rprec)                           :: wzl, wzh
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      r_cyl = cart_to_cyl(x_cart)

      vmec_get_B_vec = 0.0
      IF (vmec_get_s(this, x_cart) .le. 1.0) THEN
         r_cyl(2) = r_cyl(2) + this%phi_offset
         r_cyl(3) = r_cyl(3) - this%z_offset

         CALL GetBcyl_WOUT(r_cyl(1), r_cyl(2), r_cyl(3),                       &
     &                     vmec_get_B_vec(1), vmec_get_B_vec(2),               &
     &                     vmec_get_B_vec(3))
      ELSE IF (lfreeb) THEN
         IF (r_cyl(1) .gt. rmaxb .or. r_cyl(1) .lt. rminb) THEN
            STOP 'Radial point outside vacuum grid.'
         END IF
         IF (r_cyl(3) .gt. zmaxb .or. r_cyl(3) .lt. zminb) THEN
            STOP 'Radial point outside vacuum grid.'
         END IF

!  Normalize phi position to fall within 0 - 2Pi.
         v_norm = r_cyl(2)
         DO WHILE(v_norm .gt. twopi .or. v_norm .lt. 0)
            IF (v_norm .gt. twopi) THEN
               v_norm = v_norm - twopi
            ELSE
               v_norm = v_norm + twopi
            END IF
         END DO

!  Normalize phi position to fall within the first field period. We already
!  made sure that v_norm is greater than zero degrees.
         DO WHILE(v_norm .gt. twopi/nfp)
            v_norm = v_norm - twopi/nfp
         END DO

!  Need to interpolate the vacuum field grid to get the vacuum field external to
!  the plasma.
         delphib = twopi/(nfp*np0b)

         i_low = INT((r_cyl(1) - rminb)/delrb) + 1
         j_low = INT(v_norm/delphib) + 1
         k_low = INT((r_cyl(3) - zminb)/delzb) + 1

         r_low = rminb + (i_low - 1.0)*delrb
         phi_low = (j_low - 1.0)*delphib
         z_low = zminb + (k_low - 1.0)*delzb

         wrh = (r_cyl(1) - r_low)/delrb
         wrl = (delrb - r_cyl(1) + r_low)/delrb
         wphih = (v_norm - phi_low)/delphib
         wphil = (delphib - v_norm + phi_low)/delphib
         wzh = (r_cyl(3) - z_low)/delzb
         wzl = (delzb - r_cyl(3) + z_low)/delzb

         vmec_get_B_vec(1)                                                     &
     &      = wrl*wphil*wzl*brvac(i_low,k_low,j_low)                           &
     &      + wrh*wphil*wzl*brvac(i_low + 1,k_low,j_low)                       &
     &      + wrl*wphih*wzl*brvac(i_low,k_low,j_low + 1)                       &
     &      + wrh*wphih*wzl*brvac(i_low + 1,k_low,j_low + 1)                   &
     &      + wrl*wphil*wzh*brvac(i_low,k_low + 1,j_low)                       &
     &      + wrh*wphil*wzh*brvac(i_low + 1,k_low + 1,j_low)                   &
     &      + wrl*wphih*wzh*brvac(i_low,k_low + 1,j_low + 1)                   &
     &      + wrh*wphih*wzh*brvac(i_low + 1,k_low + 1,j_low + 1)

         vmec_get_B_vec(2)                                                     &
     &      = wrl*wphil*wzl*bpvac(i_low,k_low,j_low)                           &
     &      + wrh*wphil*wzl*bpvac(i_low + 1,k_low,j_low)                       &
     &      + wrl*wphih*wzl*bpvac(i_low,k_low,j_low + 1)                       &
     &      + wrh*wphih*wzl*bpvac(i_low + 1,k_low,j_low + 1)                   &
     &      + wrl*wphil*wzh*bpvac(i_low,k_low + 1,j_low)                       &
     &      + wrh*wphil*wzh*bpvac(i_low + 1,k_low + 1,j_low)                   &
     &      + wrl*wphih*wzh*bpvac(i_low,k_low + 1,j_low + 1)                   &
     &      + wrh*wphih*wzh*bpvac(i_low + 1,k_low + 1,j_low + 1)

         vmec_get_B_vec(3)                                                     &
     &      = wrl*wphil*wzl*bzvac(i_low,k_low,j_low)                           &
     &      + wrh*wphil*wzl*bzvac(i_low + 1,k_low,j_low)                       &
     &      + wrl*wphih*wzl*bzvac(i_low,k_low,j_low + 1)                       &
     &      + wrh*wphih*wzl*bzvac(i_low + 1,k_low,j_low + 1)                   &
     &      + wrl*wphil*wzh*bzvac(i_low,k_low + 1,j_low)                       &
     &      + wrh*wphil*wzh*bzvac(i_low + 1,k_low + 1,j_low)                   &
     &      + wrl*wphih*wzh*bzvac(i_low,k_low + 1,j_low + 1)                   &
     &      + wrh*wphih*wzh*bzvac(i_low + 1,k_low + 1,j_low + 1)

         vmec_get_B_vec = vmec_get_ext_b_plasma(this, r_cyl, .false.)          &
     &                  + vmec_get_B_vec
      END IF

      vmec_get_B_vec = cyl_to_cart_vec(r_cyl, vmec_get_B_vec)

      CALL profiler_set_stop_time('vmec_get_B_vec', start_time)

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
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      ds = 1.0/(ns - 1)
      DO i = 0, ns - 1
         s(i + 1) = i*ds
      END DO

      CALL line_seg(r, bsubv00c, s, bsubvmnc(1,:), ns)

      vmec_get_Int_B_dphi = twopi*bsubv00c

      CALL profiler_set_stop_time('vmec_get_B_vec', start_time)

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
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      s_flux = 0.0
      fmin = 0.0
      info = 0
      nfe = 0

      r_cyl = cart_to_cyl(x_cart)
      r_cyl(2) = nfp*(r_cyl(2) + this%phi_offset)
      r_cyl(3) = r_cyl(3) - this%z_offset

      CALL cyl2flx(rzl_local, r_cyl, s_flux, ns, ntor, mpol, ntmax,            &
     &             lthreed, lasym, info, nfe, fmin)

      vmec_get_s = s_flux(1)

      IF (vmec_get_s .le. 1.0) THEN
         IF (fmin .gt.  1.0E-8_rprec) THEN
            WRITE(*,*) "Warn: VMEC could not converge to flux position."       &
            WRITE(*,*) info, fmin, vmec_get_s
         END IF
      END IF

      CALL profiler_set_stop_time('vmec_get_s', start_time)

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
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      vmec_get_plasma_edge = MAX(ns, 2*mpol + 6)

      ALLOCATE(r(vmec_get_plasma_edge), z(vmec_get_plasma_edge))

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
      z = z + this%z_offset
      DEALLOCATE(cosu, cosv, sinu, sinv)

      CALL profiler_set_stop_time('vmec_get_plasma_edge', start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_magnetic_volume_rgrid => this%magnetic_cache%rsuv

      CALL profiler_set_stop_time('vmec_get_magnetic_volume_rgrid',            &
     &                            start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_magnetic_volume_zgrid => this%magnetic_cache%zsuv

      CALL profiler_set_stop_time('vmec_get_magnetic_volume_zgrid',            &
     &                            start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_magnetic_volume_jrgrid => this%magnetic_cache%jrsuv

      CALL profiler_set_stop_time('vmec_get_magnetic_volume_jrgrid',           &
     &                            start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_magnetic_volume_jphigrid => this%magnetic_cache%jphisuv

      CALL profiler_set_stop_time('vmec_get_magnetic_volume_jphigrid',         &
     &                            start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_magnetic_volume_jzgrid => this%magnetic_cache%jzsuv

      CALL profiler_set_stop_time('vmec_get_magnetic_volume_jzgrid',           &
     &                            start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_volume_int_element = isigng*this%magnetic_cache%ds              &
     &                            * this%magnetic_cache%du                     &
     &                            * this%magnetic_cache%dv

      CALL profiler_set_stop_time('vmec_get_volume_int_element',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the conducting surface integration kr grid points.
!>
!>  This method overrides
!>  @ref equilibrium::equilibrium_get_con_surface_krgrid.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_con_surface_krgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:), POINTER ::                                 &
     &   vmec_get_con_surface_krgrid
      TYPE (vmec_class), INTENT(in)         :: this

!  local variables
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_con_surface_krgrid => this%magnetic_cache%kruv

      CALL profiler_set_stop_time('vmec_get_con_surface_krgrid',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the conducting surface integration kphi grid points.
!>
!>  This method overrides
!>  @ref equilibrium::equilibrium_get_con_surface_kphigrid.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_con_surface_kphigrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:), POINTER ::                                 &
     &   vmec_get_con_surface_kphigrid
      TYPE (vmec_class), INTENT(in)         :: this

!  local variables
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_con_surface_kphigrid => this%magnetic_cache%kphiuv

      CALL profiler_set_stop_time('vmec_get_con_surface_kphigrid',             &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the conducting surface integration kz grid points.
!>
!>  This method overrides
!>  @ref equilibrium::equilibrium_get_con_surface_kzgrid.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The radial grid points.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_con_surface_kzgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:), POINTER ::                                 &
     &   vmec_get_con_surface_kzgrid
      TYPE (vmec_class), INTENT(in)         :: this

!  local variables
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_con_surface_kzgrid => this%magnetic_cache%kzuv

      CALL profiler_set_stop_time('vmec_get_con_surface_kzgrid',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get area integration element.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_area_int_element.
!>  For vmec this is |J||e^s|*du*dv where J is the jacobian. A 1/J term is
!>  integrated into kr, kphi and kz. Since due to the resulting |J|/J term, the
!>  the sign of the jacobian must be taken into account.
!>  @see vmec_equilibrium::vmec_set_magnetic_cache_calc equation 3.
!>
!>  * du = 2*Pi/number of u grid
!>  * dv = 2*Pi/number of field periods/number of toroidal planes
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns The area integration element.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_area_int_element(this)
      USE read_wout_mod, only: isigng

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: vmec_get_area_int_element
      TYPE (vmec_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_area_int_element = this%magnetic_cache%du_a
     &                          * this%magnetic_cache%dv_a*isigng

      CALL profiler_set_stop_time('vmec_get_area_int_element',                 &
     &                            start_time)

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
!>  @param[in]  num_currents   Forces the number of currents to return if
!>                             greater than zero.
!>  @param[out] scale_currents Informs the caller that currents need to be
!>                             scaled.
!>  @returns The external currents.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_ext_currents(this, num_currents,                       &
     &                               scale_currents)
      USE vmec_input, only: lfreeb
      USE mgrid_mod, only: mgrid_mode, nextcur
      USE vmec_input, only: extcur

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:), POINTER :: vmec_get_ext_currents
      TYPE (vmec_class), INTENT(in)       :: this
      INTEGER                             :: num_currents
      LOGICAL, INTENT(out)                :: scale_currents

!  local variables
      REAL (rprec)                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (lfreeb) THEN
         SELECT CASE (mgrid_mode)

            CASE ('S', 's', 'N', 'n')
               scale_currents = .true.

            CASE DEFAULT
               scale_currents = .false.

         END SELECT
         vmec_get_ext_currents => extcur(1:nextcur)
      ELSE IF (num_currents .gt. 0) THEN
!  Assume if the currents are being forced that currents are not scaled.
         scale_currents = .false.
         vmec_get_ext_currents => extcur(1:num_currents)
      ELSE
         vmec_get_ext_currents => null()
      END IF

      CALL profiler_set_stop_time('vmec_get_ext_currents', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get external plasma magnetic field.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_ext_b_plasma. The b
!>  field is returned in cyclindical coordinates. This function does not require
!>  shifting the plasma relative to the diagnostics. The plasma has already been
!>  shifted when computing the surface currents and prime position.
!>
!>  @param[in] this     A @ref vmec_class instance.
!>  @param[in] position Position to compute the fields at in cylindrical
!>                      coordinates.
!>  @param[in] axi_only Gives only the axisymmtric component of the magnetic
!>                      field.
!>  @returns The external currents.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_ext_b_plasma(this, position, axi_only)
      USE coordinate_utilities
      USE read_wout_mod, only: isigng

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)              :: vmec_get_ext_b_plasma
      TYPE (vmec_class), INTENT(in)           :: this
      REAL (rprec), DIMENSION(3), INTENT(in)  :: position
      LOGICAL, INTENT(in)                     :: axi_only

!  local variables
      REAL (rprec), DIMENSION(3)              :: x_cart
      REAL (rprec), DIMENSION(3)              :: r_vec
      INTEGER                                 :: u, v
      REAL (rprec)                            :: length
      REAL (rprec), DIMENSION(:,:,:), POINTER :: x_prime
      REAL (rprec), DIMENSION(:,:), POINTER   :: kxuv
      REAL (rprec), DIMENSION(:,:), POINTER   :: kyuv
      REAL (rprec), DIMENSION(:,:), POINTER   :: kzuv
      REAL (rprec), DIMENSION(3)              :: sum_b
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Convert the position to cartesian coordinates. Rotate and shift the signal
!  position with respect to the fixed plasma.
      x_cart = cyl_to_cart(position - (/ 0.0_dp,                               &
     &                                   this%phi_offset,                      &
     &                                   this%z_offset /))

      IF (axi_only) THEN
         CALL assert(ASSOCIATED(this%magnetic_cache%kxuv_axi),                 &
     &               'Cannot get the axisymmtric only ext bfield.')

         x_prime => this%magnetic_cache%x_axi
         kxuv => this%magnetic_cache%kxuv_axi
         kyuv => this%magnetic_cache%kyuv_axi
         kzuv => this%magnetic_cache%kzuv_axi
      ELSE
         CALL assert(ASSOCIATED(this%magnetic_cache%kxuv_full),                &
     &               'Cannot get the ext bfield.')

         x_prime => this%magnetic_cache%x_prime
         kxuv => this%magnetic_cache%kxuv_full
         kyuv => this%magnetic_cache%kyuv_full
         kzuv => this%magnetic_cache%kzuv_full
      END IF

!  The magnetic field produced by a sheet of current is
!
!    B(x) = Int K x (x-x')/|x-x'|^3 da                                       (1)
      sum_b = 0.0
      DO v = 1, SIZE(kxuv,2)
         DO u = 1, SIZE(kxuv,1)
            r_vec = x_cart - x_prime(u,v,:)
            length = SQRT(DOT_PRODUCT(r_vec, r_vec))**3.0

            sum_b(1) = sum_b(1)                                                &
     &               + (kyuv(u,v)*r_vec(3) - kzuv(u,v)*r_vec(2))               &
     &               / length
            sum_b(2) = sum_b(2)                                                &
     &               + (kzuv(u,v)*r_vec(1) - kxuv(u,v)*r_vec(3))               &
     &               / length
            sum_b(3) = sum_b(3)                                                &
     &               + (kxuv(u,v)*r_vec(2) - kyuv(u,v)*r_vec(1))               &
     &               / length
         END DO
      END DO

!  The virtual casing principal acts as if there is a perfectly conducting shell
!  around the last closed flux surface. In this shell, eddy currents induced
!  from the volume currents get canceled out by the shell currents.
!
!    B_v + B_s = 0                                                           (2)
!
!  Thus the shell field is the same magnetude but opposite sign to the volume
!  current.
!
!    B_v = -B_s                                                              (3)
!
!  After the surface integral is completed, negate the result to produce a field
!  that would be the same as a volume integral.
      vmec_get_ext_b_plasma = -mu0/(4.0*pi)*sum_b*isigng                       &
     &                      * this%magnetic_cache%du_full                      &
     &                      * this%magnetic_cache%dv_full

      vmec_get_ext_b_plasma = cart_to_cyl_vec(x_cart,                          &
     &                                        vmec_get_ext_b_plasma)

      CALL profiler_set_stop_time('vmec_get_ext_b_plasma', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get radial grid size.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_grid_size.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns Size of the radial grid.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_grid_size(this)
      USE vmec_input, only: ns_array

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: vmec_get_grid_size
      TYPE (vmec_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_grid_size = ns_array(this%ns_index)
      CALL profiler_set_stop_time('vmec_get_grid_size', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get start of the radial grid.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_grid_start.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns Start of the radial grid.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_grid_start(this)

      IMPLICIT NONE

!  Declare Arguments
      REAl (rprec)                  :: vmec_get_grid_start
      TYPE (vmec_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_grid_start = 0.0
      CALL profiler_set_stop_time('vmec_get_grid_start', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get radial grid size.
!>
!>  This method overrides @ref equilibrium::equilibrium_get_grid_size.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns Size of the radial grid.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_grid_stop(this)

      IMPLICIT NONE

!  Declare Arguments
      REAl (rprec)                  :: vmec_get_grid_stop
      TYPE (vmec_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      vmec_get_grid_stop = 1.0
      CALL profiler_set_stop_time('vmec_get_grid_stop', start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      SELECT CASE (id)

         CASE (vmec_bloat_id, vmec_curtor_id, vmec_phiedge_id,                 &
     &         vmec_pres_scale_id, vmec_vvc_smaleli_id,                        &
     &         vmec_vvc_kappa_p_id, vmec_betatot_id, vmec_betapol_id,          &
     &         vmec_betator_id, vmec_betaxis_id, vmec_phi_offset_id,           &
     &         vmec_z_offset_id)
            vmec_is_scaler_value = .true.

         CASE DEFAULT
            vmec_is_scaler_value = .false.

      END SELECT

      CALL profiler_set_stop_time('vmec_is_scaler_value', start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
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

      CALL profiler_set_stop_time('vmec_is_1d_array', start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
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

      CALL profiler_set_stop_time('vmec_is_2d_array', start_time)

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

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
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
     &         vmec_phi_offset_id, vmec_z_offset_id)
            vmec_is_recon_param = .true.

         CASE DEFAULT
            vmec_is_recon_param = .false.

      END SELECT

      CALL profiler_set_stop_time('vmec_is_recon_param', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if a the point magnetics are being used.
!>
!>  This method overrides @ref equilibrium::equilibrium_is_using_point.
!>
!>  @param[in] this A @ref vmec_class instance.
!>  @returns True if the point magnetic are being used.
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_using_point(this)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_using_point
      TYPE (vmec_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Some compilers will generate code that will early terminate when a first
!  operand to a .and. operator evaluates to false. Others are not so clever.
!  Avoid running the second ASSOCIATED check if the magnetic cache is null.
      vmec_is_using_point = ASSOCIATED(this%magnetic_cache)
      IF (vmec_is_using_point) THEN
         vmec_is_using_point = ASSOCIATED(this%magnetic_cache%kxuv_full)
      END IF

      CALL profiler_set_stop_time('vmec_is_using_point', start_time)

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
      USE vmec_input, only: delt, niter, raxis_cc, raxis_cs,                   &
     &                      zaxis_cc, zaxis_cs
      USE vmec_params, only: timestep_flag, output_flag,                       &
     &                       reset_jacdt_flag, jac75_flag,                     &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag, restart_flag
      USE xstuff, only: xcdot
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE vmec_main, only: iterc, iter1, iter2

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                            :: vmec_converge
      TYPE (vmec_class), INTENT(inout)   :: this
      INTEGER, INTENT(inout)             :: num_iter
      INTEGER, INTENT(in)                :: iou

!  local variables
!  delt_local   Value of delt on entry. Store, for restoration after error
!               recovery.
!  ictrl_array  Array containing various vmec control parameters.
      REAL (rprec)                       :: delt_local
      INTEGER                            :: niter_local
      INTEGER, DIMENSION(5)              :: ictrl_array
      INTEGER                            :: recovery_index
      INTEGER                            :: error, num_iter_init
      TYPE (vmec_context_class), POINTER :: temp_context => null()
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      delt_local = delt
      niter_local = niter
      num_iter_init = iterc
      recovery_index = 0

!  Store the current VMEC context. This context is different than the one stored
!  in save/reset state methods. This context needs to be stored so that error
!  recovery can restore the state to before runvmec was run but after the
!  parameters were changed.
      temp_context => vmec_context_construct()

!  Run vmec until successful or an unrecoverable error has occured.
      DO
         ictrl_array = 0
         ictrl_array(1) = IOR(ictrl_array(1), timestep_flag)
         ictrl_array(1) = IOR(ictrl_array(1), output_flag)
         ictrl_array(1) = IOR(ictrl_array(1), reset_jacdt_flag)

!  Allow multigrid for the inital convergence only.
         IF (num_iter .gt. 1) THEN
            ictrl_array(3) = niter_local
            ictrl_array(4) = this%ns_index
         ELSE
!  This allows the inital convergence to match a stand alone VMEC run.
            ictrl_array(1) = IOR(ictrl_array(1), restart_flag)
         END IF

!  Zero the xcdot array and store the xc array used for timestep resets.
         xcdot = 0.0
         xstore = xc
         CALL reset_params
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

!  Restore the context to the state it was before runvmec was called.
               CALL vmec_context_set_context(temp_context)

               IF (recovery_index .le. 2) THEN
                  WRITE (*,1000) 'jac75_flag'
                  WRITE (*,1001) 'delt', delt, delt/2.0
                  WRITE (iou,1000) 'jac75_flag'
                  WRITE (iou,1001) 'delt', delt, delt/2.0
                  delt = delt/2.0
               ELSE
                  EXIT
               END IF

            CASE (more_iter_flag)

!  Restore the context to the state it was before runvmec was called.
            CALL vmec_context_set_context(temp_context)

            IF (recovery_index .le. 2) THEN
               WRITE (*,1000) 'more_iter_flag'
               WRITE (*,1002) 'niter', niter_local,                            &
     &                        niter_local + niter
               WRITE (iou,1000) 'more_iter_flag'
               WRITE (iou,1002) 'niter', niter_local,                          &
     &                          niter_local + niter

               niter_local = niter_local + niter
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

!  No longer need the local context array.
      DEALLOCATE(temp_context)

!  Check for errors, if vmec was successful, load vmec data for use with the
!  signals. Other wise report the error.
      vmec_converge = .false.
      SELECT CASE (ictrl_array(2))

         CASE (successful_term_flag, norm_term_flag)
            CALL read_wout_file(this%wout_file_name, error)
            CALL assert_eq(error, 0, "Error reading wout file.")
            CALL LoadRZL

            vmec_converge = .true.
            IF (ASSOCIATED(this%magnetic_cache)) THEN
               CALL vmec_set_magnetic_cache(this)
            END IF

         CASE DEFAULT
!>  @todo FIXME: Put in error reporting here. For now leave this blank until I
!>  can figure how v3fit decides if vmec has converged.

      END SELECT

!  Increment the number of iterations by the number of iterations used by
!  runvmec.
      num_iter = num_iter + iterc - num_iter_init
      delt = delt_local

      CALL profiler_set_stop_time('vmec_converge', start_time)

1000  FORMAT('VMEC CONVERGENCE ERROR: ',a)
1001  FORMAT(a,' change: old, new: ',es12.5,' ',es12.5)
1002  FORMAT(a,' change: old, new: ',i12,' ',i12)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Save the internal state of the equilibrium.
!>
!>  Saves the VMEC state and the wout file.
!>
!>  @param[inout] this A @ref vmec_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_save_state(this)
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this

!  local variables
      INTEGER                          :: error
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (ASSOCIATED(this%vmec_context_save)) THEN
         CALL vmec_context_get_context(this%vmec_context_save)
      ELSE
         this%vmec_context_save => vmec_context_construct()
      END IF

!  Cache the wout file by appending _cache to a copy. Use copy here to keep the
!  the orginal file intact.
      CALL copy_file(this%wout_file_name,                                      &
     &               TRIM(this%wout_file_name) // '_cache',                    &
     &               error)
      CALL assert_eq(error, 0, 'Error copying wout file.')

      CALL profiler_set_stop_time('vmec_save_state', start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Reset the internal state of the equilibrium.
!>
!>  Copies the saved VMEC state back to VMEC.
!>
!>  @param[inout] this A @ref vmec_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_reset_state(this)
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this

!  local variables
      INTEGER                          :: error
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (ASSOCIATED(this%vmec_context_save)) THEN
         CALL vmec_context_set_context(this%vmec_context_save)
      END IF

!  Reset the wout file.
      CALL copy_file(TRIM(this%wout_file_name) // '_cache',                    &
     &               TRIM(this%wout_file_name), error)
      CALL assert_eq(error, 0, 'Error moving wout file.')

      CALL read_wout_file(this%wout_file_name, error)
      CALL assert_eq(error, 0, 'Error reading reset wout file.')
      CALL LoadRZL

      IF (ASSOCIATED(this%magnetic_cache)) THEN
         CALL vmec_set_magnetic_cache(this)
      END IF

      CALL profiler_set_stop_time('vmec_reset_state', start_time)

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
      USE safe_open_mod
      USE vmec_input

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: iou

!  local variables
      INTEGER                       :: i
      INTEGER                       :: iou_nl
      INTEGER                       :: status
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()
      WRITE (iou,*)
      WRITE (iou,*) 'Equilibrium Type : VMEC'
      CALL pprofile_write(this%ne, 'pp_ne', iou)

      DO i = 1, SIZE(this%sxrem)
         CALL pprofile_write(this%sxrem(i)%p, 'pp_sxrem', iou)
      END DO

      CALL pprofile_write(this%te, 'pp_te', iou)

!  Update the namelist input variables from the equilibrium solution.
      CALL vmec_set_namelist(this)

      iou_nl = 0
      CALL safe_open(iou_nl, status,                                           &
     &               TRIM(this%vmec_file_name) // '_out',                      &
     &               'replace', 'formatted', delim_in='none')
      WRITE (iou_nl, nml=indata)
      CLOSE (iou_nl, iostat=status)

      CALL profiler_set_stop_time('vmec_write', start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write the current valid input.
!>
!>  This method overrides @ref equilibrium::equilibrium_write_input. The
!>  boundary and other fixed parameters do not get updated as the reconstruction
!>  progresses. Need to update them manually if in free boundary mode.
!>
!>  @param[in] this         A @ref vmec_class instance.
!>  @param[in] current_step Step number to append to input filename.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_write_input(this, current_step)
      USE safe_open_mod
      USE vmec_input
      USE data_parameters
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in)        :: this
      INTEGER, INTENT(in)                  :: current_step

!  local variables
      INTEGER                              :: iou_nl
      INTEGER                              :: status
      CHARACTER (len=data_filename_length) :: filename
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Update the namelist input variables from the equilibrium solution.
      CALL vmec_set_namelist(this)

!  Write out the final namelist input file.
      iou_nl = 0
      status = 0

      WRITE (filename,1000) TRIM(this%vmec_file_name), current_step

      CALL safe_open(iou_nl, status, TRIM(filename), 'replace',                &
     &               'formatted', delim_in='none')
      WRITE (iou_nl, nml=indata)
      CLOSE (iou_nl, iostat=status)
      CALL assert_eq(status, 0, 'Error saving input file.')

!  Save the woutfile by making a copy with a new name.
      WRITE (filename,1000) TRIM(this%wout_file_name), current_step

      CALL copy_file(this%wout_file_name, TRIM(filename), status)
      CALL assert_eq(status, 0, 'Error copying wout file.')

      CALL profiler_set_stop_time('vmec_write_input', start_time)

1000  FORMAT(a,'_',i0.3)

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
!>
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
      INTEGER                       :: vmec_extcur_id
      INTEGER                       :: vmec_curtor_id
      INTEGER                       :: vmec_pres_scale_id
      INTEGER                       :: vmec_rbc_id
      INTEGER                       :: vmec_rbs_id
      INTEGER                       :: vmec_zbc_id
      INTEGER                       :: vmec_zbs_id
      INTEGER                       :: vmec_l_sub_i_id
      INTEGER                       :: vmec_kappa_p_id
      REAL (rprec)                  :: start_time


!  Start of executable code
      start_time = profiler_get_start_time()

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
     &                    (/ vmec_nextcur_dim_id /), vmec_extcur_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_curtor', nf_double, 0, 0,         &
     &                    vmec_curtor_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_pres_scale', nf_double, 0,        &
     &                    0, vmec_pres_scale_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_rbc', nf_double, 2,               &
     &                    (/ vmec_boundary_n_dim_id,                           &
     &                       vmec_boundary_m_dim_id /), vmec_rbc_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_rbs', nf_double, 2,               &
     &                    (/ vmec_boundary_n_dim_id,                           &
     &                       vmec_boundary_m_dim_id /), vmec_rbs_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_zbc', nf_double, 2,               &
     &                    (/ vmec_boundary_n_dim_id,                           &
     &                       vmec_boundary_m_dim_id /), vmec_zbc_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_zbs', nf_double, 2,               &
     &                    (/ vmec_boundary_n_dim_id,                           &
     &                       vmec_boundary_m_dim_id /), vmec_zbs_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_l_sub_i', nf_double, 1,           &
     &                    (/ maxnsetps_dim_id /), vmec_l_sub_i_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_def_var(result_ncid, 'vmec_kappa_p', nf_double, 1,           &
     &                    (/ maxnsetps_dim_id /), vmec_kappa_p_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      CALL profiler_set_stop_time('vmec_def_result', start_time)

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
      INTEGER                       :: vmec_extcur_id
      INTEGER                       :: vmec_curtor_id
      INTEGER                       :: vmec_pres_scale_id
      INTEGER                       :: vmec_rbc_id
      INTEGER                       :: vmec_rbs_id
      INTEGER                       :: vmec_zbc_id
      INTEGER                       :: vmec_zbs_id
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      status = nf_inq_varid(result_ncid, 'vmec_extcur', vmec_extcur_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var_double(result_ncid, vmec_extcur_id,                  &
     &                           extcur(1:nextcur))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'vmec_curtor', vmec_curtor_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var_double(result_ncid, vmec_curtor_id,                  &
     &                           curtor)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'vmec_pres_scale',                    &
     &                      vmec_pres_scale_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var_double(result_ncid, vmec_pres_scale_id,              &
     &                           pres_scale)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      boundary_size(1) = 2*ntor + 1
      boundary_size(2) = mpol1 + 1

      status = nf_inq_varid(result_ncid, 'vmec_rbc', vmec_rbc_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_vara_double(result_ncid, vmec_rbc_id,                    &
     &                            (/ 1, 1 /), boundary_size,                   &
     &                            rbc(-ntor:ntor,0:mpol1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'vmec_rbs', vmec_rbs_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_vara_double(result_ncid, vmec_rbs_id,                    &
     &                            (/ 1, 1 /), boundary_size,                   &
     &                            rbs(-ntor:ntor,0:mpol1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'vmec_zbc', vmec_zbc_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_vara_double(result_ncid, vmec_zbc_id,                    &
     &                            (/ 1, 1 /), boundary_size,                   &
     &                            zbc(-ntor:ntor,0:mpol1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'vmec_zbs', vmec_zbs_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_vara_double(result_ncid, vmec_zbs_id,                    &
     &                            (/ 1, 1 /), boundary_size,                   &
     &                            zbs(-ntor:ntor,0:mpol1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      CALL vmec_write_step_data(this, result_ncid, 1)

      CALL profiler_set_stop_time('vmec_write_init_data', start_time)

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
      INTEGER                       :: vmec_l_sub_i_id
      INTEGER                       :: vmec_kappa_p_id
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      status = nf_inq_varid(result_ncid, 'vmec_l_sub_i',                       &
     &                      vmec_l_sub_i_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var1_double(result_ncid, vmec_l_sub_i_id,                &
     &                            current_step, vvc_smaleli)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_inq_varid(result_ncid, 'vmec_kappa_p',                       &
     &                      vmec_kappa_p_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))
      status = nf_put_var1_double(result_ncid, vmec_kappa_p_id,                &
     &                            current_step, vvc_kappa_p)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      CALL profiler_set_stop_time('vmec_write_step_data', start_time)

      END SUBROUTINE

!*******************************************************************************
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize data to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref vmec_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_sync_children(this, intra_comm)
      USE file_opts
      USE vmec_input, only: lfreeb, mgrid_file

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(in)              :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                          :: error
      INTEGER                          :: i
      INTEGER                          :: num_sxrem
      INTEGER                          :: mpi_rank
      INTEGER                          :: mpi_size
      LOGICAL                          :: sync_cache

!  Start of executable code
      CALL MPI_COMM_SIZE(intra_comm, mpi_size, error)
      CALL MPI_COMM_RANK(intra_comm, mpi_rank, error)

!  Check all relevant files for relative paths and copy them to the sub process
!  directories. This should only be handeled by the main process.
      IF (mpi_rank .eq. 0) THEN
         IF (.not.is_absolute_path(this%vmec_file_name)) THEN
            DO i = 2, mpi_size
               CALL copy_file(this%vmec_file_name, process_dir(i),             &
     &                        error)
            END DO
         END IF

         IF (lfreeb .and. .not.is_absolute_path(mgrid_file)) THEN
            DO i = 2, mpi_size
               CALL copy_file(mgrid_file, process_dir(i), error)
            END DO
         END IF
      END IF

      CALL MPI_BCAST(this%phi_offset, 1, MPI_REAL8, 0, intra_comm,             &
     &               error)
      CALL MPI_BCAST(this%z_offset, 1, MPI_REAL8, 0, intra_comm, error)

      CALL MPI_BCAST(this%pol_rad_ratio, 1, MPI_REAL8, 0, intra_comm,          &
     &               error)

      CALL MPI_BCAST(this%vmec_file_name, data_filename_length,                &
     &               MPI_CHARACTER, 0, intra_comm, error)

      CALL pprofile_sync_children(this%ne, intra_comm)
      CALL pprofile_sync_children(this%te, intra_comm)

      num_sxrem = SIZE(this%sxrem)
      CALL MPI_BCAST(num_sxrem, 1, MPI_INTEGER, 0, intra_comm, error)
      DO i = 1, num_sxrem
         CALL pprofile_sync_children(this%sxrem(i)%p, intra_comm)
      END DO

      CALL MPI_BCAST(this%ns_index, 1, MPI_INTEGER, 0, intra_comm,             &
     &               error)

      sync_cache = ASSOCIATED(this%magnetic_cache)
      CALL MPI_BCAST(sync_cache, 1, MPI_LOGICAL, 0, intra_comm, error)
      IF (sync_cache) THEN
         CALL vmec_set_magnetic_cache(this, intra_comm)
      END IF
#endif

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Syncronize the vmec state to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref vmec_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_sync_state(this, intra_comm)
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
       TYPE (vmec_class), INTENT(inout) :: this
       INTEGER, INTENT(in)              :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                           :: error
      INTEGER                           :: mpi_rank
      INTEGER                           :: mpi_size
      INTEGER                           :: i

!  Start of executable code
      CALL MPI_COMM_SIZE(intra_comm, mpi_size, error)
      CALL MPI_COMM_RANK(intra_comm, mpi_rank, error)

!  If this is the parent process, copy the wout file to the children.
      IF (mpi_rank .eq. 0) THEN
         DO i = 2, mpi_size
            CALL copy_file(TRIM(this%wout_file_name), process_dir(i),          &
     &                     error)
            CALL assert_eq(error, 0, 'Error syncing wout file.')
         END DO
      END IF

      CALL MPI_BARRIER(intra_comm, error)

      IF (mpi_rank .gt. 0) THEN
         CALL read_wout_file(this%wout_file_name, error)
         CALL assert_eq(error, 0, 'Error reading synced wout file.')
         CALL LoadRZL

         IF (ASSOCIATED(this%magnetic_cache)) THEN
            CALL vmec_set_magnetic_cache(this)
         END IF
      END IF

      CALL vmec_context_sync_state(this%vmec_context_save, intra_comm)

#endif
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
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

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

      CALL profiler_set_stop_time('vmec_reset_boundary', start_time)

      END SUBROUTINE

      END MODULE
