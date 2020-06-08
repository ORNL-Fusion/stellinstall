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
!>     @item{jcuru,       Poloidal current density.,                       @ref read_wout_mod::jcuru}
!>     @item{jcurv,       Toroidal current density.,                       @ref read_wout_mod::jcurv}
!>     @item{jdotb,       Current density in direction of b field.,        @ref read_wout_mod::jdotb}
!>     @item{raxis_cc,    Magnetic axis R cosine coeffiecents.,            @ref read_wout_mod::raxis}
!>     @item{raxis_cs,    Magnetic axis R sine coeffiecents.,              @ref read_wout_mod::raxis}
!>     @item{zaxis_cc,    Magnetic axis Z cosine coeffiecents.,            @ref read_wout_mod::zaxis}
!>     @item{zaxis_cs,    Magnetic axis Z sine coeffiecents.,              @ref read_wout_mod::zaxis}
!>     @item{qfact,       Safety factor profile.,                          @ref read_wout_mod::qfact}
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
      USE v3_utilities
      USE vmec_context

      IMPLICIT NONE

!*******************************************************************************
!  vmec equilibrium module parameters
!*******************************************************************************
!  Define id's only for values that can change. All others shouldn't be needed
!  outside of the vmec interface. Parameter id's start at 15 since 0-14 are
!  reserved for non equilibrium model parameters. These numbers will need to be
!  updated if any new model parameters are added.

!  Variable Parameters
!>  1D Array of function current profile parameters. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_id          = 15
!>  1D Array of segment current profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_aux_s_id    = 16
!>  1D Array of segment current profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ac_aux_f_id    = 17
!>  1D Array of function iotabar profile parameters. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_id          = 18
!>  1D Array of segment iotabar profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_aux_s_id    = 19
!>  1D Array of segment iotabar profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_ai_aux_f_id    = 20
!>  1D Array of function pressure profile parameters. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_id          = 21
!>  1D Array of segment pressure profile s poitions. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_aux_s_id    = 22
!>  1D Array of segment pressure profile f values. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_am_aux_f_id    = 23
!>  Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_bloat_id       = 24
!>  2D Array of outer boundary R cosine coeffiecents. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_rbc_id         = 25
!>  2D Array of outer boundary Z sine coeffiecents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_zbs_id         = 26
!>  2D Array of outer boundary R sine coeffiecents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_rbs_id         = 27
!>  2D Array of outer boundary Z cosine coeffiecents. Defined in
!>  @ref vmec_input.
      INTEGER, PARAMETER :: vmec_zbc_id         = 28
!>  1D Array of external currents. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_extcur_id      = 29
!>  The toroidal current. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_curtor_id      = 30
!>  The edge magnetic flux. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_phiedge_id     = 31
!>  The scaling parameter for the pressure. Defined in @ref vmec_input.
      INTEGER, PARAMETER :: vmec_pres_scale_id  = 32

!  Derived Parameters
!>  2D Array of surface R cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_rmnc_id        = 33
!>  2D Array of surface Z sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_zmns_id        = 34
!>  2D Array of surface lambda cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_lmns_id        = 35
!>  2D Array of surface g cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_gmnc_id        = 36
!>  2D Array of surface B_s sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubsmns_id    = 37
!>  2D Array of surface B_u cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubumnc_id    = 38
!>  2D Array of surface B_v cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubvmnc_id    = 39
!>  2D Array of surface B^u cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupumnc_id    = 40
!>  2D Array of surface B^v cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupvmnc_id    = 41
!>  2D Array of surface R sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_rmns_id        = 42
!>  2D Array of surface Z cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_zmnc_id        = 43
!>  2D Array of surface lambda sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_lmnc_id        = 44
!>  2D Array of surface g sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_gmns_id        = 45
!>  2D Array of surface B_s cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubsmnc_id    = 46
!>  2D Array of surface B_u sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubumns_id    = 47
!>  2D Array of surface B_v sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsubvmns_id    = 48
!>  2D Array of surface B^u sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupumns_id    = 49
!>  2D Array of surface B^v sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_bsupvmns_id    = 50
!>  Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_phi_id         = 51
!>  1D Array of full grid iotabar. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_iotaf_id       = 52
!>  1D Array of half grid iotabar. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_iotas_id       = 53
!>  1D Array of magnetic axis R cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_raxis_cc_id    = 54
!>  1D Array of magnetic axis R sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_raxis_cs_id    = 55
!>  1D Array of magnetic axis Z cosine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_zaxis_cc_id    = 56
!>  1D Array of magnetic axis Z sine coeffiecents. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_zaxis_cs_id    = 57
!>  1D Array of full grid safety factor. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_qfact_id       = 58

!  Derived Parameters
!>  Internal inducance. Defined in @ref v3f_vmec_comm. Computed in @ref eqfor.
      INTEGER, PARAMETER :: vmec_vvc_smaleli_id = 74
!>  Mean elongation. Defined in @ref v3f_vmec_comm. Computed in @ref eqfor.
      INTEGER, PARAMETER :: vmec_vvc_kappa_p_id = 75
!>  Total Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betatot_id     = 76
!>  Poloidal Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betapol_id     = 77
!>  Toroidal Beta. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betator_id     = 78
!>  Beta on axis. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_betaxis_id     = 79
!>  Poloidal current density. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_jcuru_id       = 80
!>  Toroidal current density. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_jcurv_id       = 81
!>  Current density in direction of b field. Defined in read_wout_mod.
      INTEGER, PARAMETER :: vmec_jdotb_id       = 82

!-------------------------------------------------------------------------------
!>  Base class representing a vmec_equilibrium.
!>  @par Super Class:
!>  @ref equilibrium
!-------------------------------------------------------------------------------
      TYPE vmec_class
!>  File name of the output of vmec.
         CHARACTER (len=path_length)            :: wout_file_name
!>  File name of the vmec namelist inout file.
         CHARACTER (len=path_length)            :: vmec_file_name

!>  VMEC context.
         TYPE (vmec_context_class), POINTER     ::                             &
     &      vmec_context_save => null()

!>  Index of the largest ns value.
         INTEGER                                :: ns_index

      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************

!-------------------------------------------------------------------------------
!>  Forward declare the interface for @ref runvmec.
!>  @note If @ref runvmec.f was implimented as a module this would not be
!>  necessary.
!-------------------------------------------------------------------------------
      INTERFACE
         SUBROUTINE runvmec(control_array, vmec_input_file,                    &
     &                      use_screen, comm, reset_file_name)
         INTEGER, INTENT(inout), TARGET :: control_array(5)
         LOGICAL, INTENT(in)            :: use_screen
         CHARACTER(LEN=*), INTENT(in)   :: vmec_input_file
         INTEGER, INTENT(in)            :: comm
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
!>  @param[in]    file_name      Filename of the vmec namelist input file.
!>  @param[in]    wout_file_name Filename of the vmec wout input file.
!>  @param[in]    iou            Input/output unit to log messages to.
!>  @param[in]    eq_comm        MPI communicator pool for VMEC.
!>  @param[in]    recon_comm     MPI communicator pool for reconstruction.
!>  @param[inout] state_flags    Bitwise flags to indicate which parts of the
!>                               model changed.
!>  @returns A pointer to a constructed @ref vmec_class object.
!-------------------------------------------------------------------------------
      FUNCTION vmec_construct_eq(file_name, iou,                               &
     &                           eq_comm, recon_comm)
      USE vmec_params, only: restart_flag, readin_flag, timestep_flag,         &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag
      USE vmec_input, only: l_v3fit, ns_array, mgrid_file, lfreeb
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER           :: vmec_construct_eq
      CHARACTER (len=*), INTENT(in)        :: file_name
      INTEGER, INTENT(in)                  :: iou
      INTEGER, INTENT(in)                  :: eq_comm
      INTEGER, INTENT(in)                  :: recon_comm

!  local variables
      INTEGER                              :: i
      INTEGER                              :: eq_rank
      INTEGER                              :: recon_rank
      INTEGER                              :: recon_size
      INTEGER, DIMENSION(5)                :: ictrl_array
      INTEGER                              :: index_dat, index_end
      INTEGER                              :: error
      REAL (rprec)                         :: start_time
      CHARACTER (len=path_length)          :: mgrid_path
      logical :: vmec_lscreen = .true.

!  Start of executable code
      start_time = profiler_get_start_time()

      eq_rank = 0
      recon_rank = 0
      recon_size = 1
#if defined(MPI_OPT)
      CALL MPI_COMM_RANK(eq_comm, eq_rank, error)

      IF (eq_rank .eq. 0) THEN
         CALL MPI_COMM_RANK(recon_comm, recon_rank, error)
         CALL MPI_COMM_SIZE(recon_comm, recon_size, error)
      END IF
      CALL MPI_BCAST(recon_rank, 1, MPI_INTEGER, 0, eq_comm, error)
#endif

      IF (eq_rank .eq. 0) THEN
         WRITE (*,*) ' *** Initializing VMEC equilibrium from file ' //        &
     &               TRIM(file_name)
         WRITE (iou,*) ' *** Initializing VMEC equilibrium from ' //           &
     &                 'file ' // TRIM(file_name)
      END IF

      ALLOCATE(vmec_construct_eq)

!  Initialize VMEC by taking a single step.
      ictrl_array = 0
      ictrl_array(1) = IOR(ictrl_array(1), restart_flag)
      ictrl_array(1) = IOR(ictrl_array(1), readin_flag)
      ictrl_array(1) = IOR(ictrl_array(1), timestep_flag)
      ictrl_array(3) = 1 ! Run for one iteration.
      ictrl_array(4) = 1 ! Limit VMEC to single grid.

      l_v3fit = .true.

      vmec_construct_eq%vmec_file_name = TRIM(file_name)

!  Only the rank zero reconstruction pool runs the inital convergence. Initalize
!  VMEC on that pool then broad cast the size of the largest NS value to the
!  children. This ensures that they are allocated with the correct grid sizes.

      IF (recon_rank .eq. 0) THEN
         write(*,*) '*** Calling runvmec from vmec_construct_eq ...'
         CALL runvmec(ictrl_array, file_name, vmec_lscreen, eq_comm)
         
         vmec_construct_eq%ns_index = MAXLOC(ns_array, 1)
         IF (eq_rank .eq. 0) THEN
            CALL MPI_BCAST(vmec_construct_eq%ns_index, 1, MPI_INTEGER,         &
     &                     0, recon_comm, error)
         END IF
      END IF

!  Check for errors in runvmec.
      SELECT CASE (ictrl_array(2))

         CASE (successful_term_flag, norm_term_flag, more_iter_flag)
!  Set the wout file name.
            index_dat = INDEX(file_name, 'input.')
            index_end = LEN_TRIM(file_name)
            IF (index_dat .gt. 0) THEN
               vmec_construct_eq%wout_file_name = 'wout_' //                   &
     &            TRIM(file_name(index_dat + 6:index_end)) // '.nc'
            ELSE
               vmec_construct_eq%wout_file_name = 'wout_' //                   &
     &            TRIM(file_name(1:index_end)) // '.nc'
            END IF

         CASE DEFAULT
            CALL err_fatal('vmec_construct runvmec command failed ',           &
     &                     int = ictrl_array(2))

      END SELECT

      vmec_construct_eq%vmec_context_save => vmec_context_construct()

      CALL profiler_set_stop_time('vmec_construct_eq', start_time)

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
!>  @param[inout] this        A @ref vmec_class instance.
!>  @param[in]    id          ID of the parameter.
!>  @param[in]    i_index     The ith index of the parameter.
!>  @param[in]    j_index     The jth index of the parameter.
!>  @param[in]    value       The value of the parameter.
!>  @param[in]    eq_comm     MPI communicator for the child equilibrium
!>                            processes.
!>  @param[inout] state_flags Bitwise flags to indicate which parts of the model
!>                            changed.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_param(this, id, i_index, j_index, value,             &
     &                          eq_comm, state_flags)
      USE vmec_input, only: rbc, zbs, rbs, zbc,                                &
     &                      ac, ac_aux_s, ac_aux_f,                            &
     &                      ai, ai_aux_s, ai_aux_f,                            &
     &                      am, am_aux_s, am_aux_f,                            &
     &                      extcur, curtor, phiedge, pres_scale,               &
     &                      bloat, lfreeb, aphi
      USE xstuff, only: xc, xcdot
      USE vmec_main, only: irzloff, currv, ivac
      USE stel_constants, only: mu0
      USE model_state

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(in)              :: id
      INTEGER, INTENT(in)              :: i_index
      INTEGER, INTENT(in)              :: j_index
      REAL (rprec), INTENT(in)         :: value
      INTEGER, INTENT(in)              :: eq_comm
      INTEGER, INTENT(inout)           :: state_flags

!  local variables
      REAL (rprec)                     :: phiedge_ratio
      INTEGER                          :: error
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (id)

         CASE (vmec_rbc_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            rbc(i_index, j_index) = value
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_zbs_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            zbs(i_index, j_index) = value
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_rbs_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            rbs(i_index, j_index) = value
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_zbc_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            zbc(i_index, j_index) = value
            IF (.not.lfreeb) THEN
               CALL vmec_reset_boundary(this)
               CALL profil3d(xc(1), xc(1 + irzloff), .true., .false.)
            ELSE
               CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)
            END IF

         CASE (vmec_ac_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            ac(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ac_aux_s_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            ac_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ac_aux_f_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            ac_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ai_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            ai(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ai_aux_s_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            ai_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_ai_aux_f_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            ai_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_am_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            am(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_am_aux_s_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            am_aux_s(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_am_aux_f_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            am_aux_f(i_index) = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_extcur_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            extcur(i_index) = value
#if defined(MPI_OPT)
            CALL MPI_BCAST(mpi_mgrid_task, 1, MPI_INTEGER, 0, eq_comm,         &
     &                     error)
            CALL MPI_BCAST(i_index, 1, MPI_INTEGER, 0, eq_comm, error)
#endif
            CALL vmec_read_vac_file(this, i_index, eq_comm)

         CASE (vmec_curtor_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            curtor = value
            currv = mu0*curtor
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_phiedge_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
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
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            pres_scale = value
            CALL profil1d(xc, xcdot, .false.)

         CASE (vmec_bloat_id)
            state_flags = IBSET(state_flags, model_state_vmec_flag)
            bloat = value
            CALL profil1d(xc, xcdot, .false.)
            CALL profil3d(xc(1), xc(1 + irzloff), .false., .false.)

         CASE DEFAULT
!>  @todo FIXME: Consider throwing a fatal error here instead of silently doing
!>  nothing.
            state_flags = state_flags

      END SELECT

      CALL profiler_set_stop_time('vmec_set_param', start_time)

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
            rbc(INT(xn(i)/nfp),INT(xm(i))) =                                   &
     &         rmnc(i,ns_array(this%ns_index))
            zbs(INT(xn(i)/nfp),INT(xm(i))) =                                   &
     &         zmns(i,ns_array(this%ns_index))
            IF (lasym) THEN
               rbs(INT(xn(i)/nfp),INT(xm(i))) =                                &
     &            rmns(i,ns_array(this%ns_index))
               zbc(INT(xn(i)/nfp),INT(xm(i))) =                                &
     &            zmnc(i,ns_array(this%ns_index))
            END IF
         END DO
      END IF

!  Update the magnetic axis.
      raxis_cc(0:ntor) = rmnc(1:ntor + 1,1)
      zaxis_cs(0:ntor) = zmns(1:ntor + 1,1)
      IF (lasym) THEN
         raxis_cs(0:ntor) = rmns(1:ntor + 1,1)
         zaxis_cc(0:ntor) = zmnc(1:ntor + 1,1)
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

         CASE ('raxis_cc')
            vmec_get_param_id = vmec_raxis_cc_id

         CASE ('raxis_cs')
            vmec_get_param_id = vmec_raxis_cs_id

         CASE ('zaxis_cc')
            vmec_get_param_id = vmec_zaxis_cc_id

         CASE ('zaxis_cs')
            vmec_get_param_id = vmec_zaxis_cs_id

         CASE ('qfact')
            vmec_get_param_id = vmec_qfact_id

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

         CASE ('jcuru')
            vmec_get_param_id = vmec_jcuru_id

         CASE ('jcurv')
            vmec_get_param_id = vmec_jcurv_id

         CASE ('jdotb')
            vmec_get_param_id = vmec_jdotb_id

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
     &                         betatot, betapol, betator, betaxis,             &
     &                         jcuru, jcurv, jdotb, raxis, zaxis, qfact
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

         CASE (vmec_raxis_cc_id)
            vmec_get_param_value = raxis(i_index, 1)

         CASE (vmec_raxis_cs_id)
            vmec_get_param_value = raxis(i_index, 2)

         CASE (vmec_zaxis_cc_id)
            vmec_get_param_value = zaxis(i_index, 2)

         CASE (vmec_zaxis_cs_id)
            vmec_get_param_value = zaxis(i_index, 1)

         CASE (vmec_qfact_id)
            vmec_get_param_value = qfact(i_index)

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

         CASE (vmec_jcuru_id)
            vmec_get_param_value = jcuru(i_index)

         CASE (vmec_jcurv_id)
            vmec_get_param_value = jcurv(i_index)

         CASE (vmec_jdotb_id)
            vmec_get_param_value = jdotb(i_index)

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

         CASE (vmec_raxis_cc_id)
            vmec_get_param_name = 'raxis_cc'

         CASE (vmec_raxis_cs_id)
            vmec_get_param_name = 'raxis_cs'

         CASE (vmec_zaxis_cc_id)
            vmec_get_param_name = 'zaxis_cc'

         CASE (vmec_zaxis_cs_id)
            vmec_get_param_name = 'zaxis_cs'

         CASE (vmec_qfact_id)
            vmec_get_param_name = 'qfact'

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

         CASE (vmec_jcuru_id)
            vmec_get_param_name = 'jcuru'

         CASE (vmec_jcurv_id)
            vmec_get_param_name = 'jcurv'

         CASE (vmec_jdotb_id)
            vmec_get_param_name = 'jdotb'

         CASE DEFAULT
            WRITE (*,*) id
            STOP 'Invalid VMEC parameter id.'

      END SELECT

      CALL profiler_set_stop_time('vmec_get_param_name', start_time)

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
!>  @param[in]    eq_comm  MPI communicator pool for VMEC.
!>  @returns True if the convergece was sucessful and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION vmec_converge(this, num_iter, iou, eq_comm)
      USE vmec_input, only: delt, niter, nstep
      USE vmec_params, only: timestep_flag, output_flag,                       &
     &                       reset_jacdt_flag, jac75_flag,                     &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag, restart_flag
      USE xstuff, only: xcdot, xc, xstore
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE vmec_main, only: iterc, iter1, iter2

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                            :: vmec_converge
      logical :: vmec_lscreen = .true.
      TYPE (vmec_class), INTENT(inout)   :: this
      INTEGER, INTENT(inout)             :: num_iter
      INTEGER, INTENT(in)                :: iou
      INTEGER, INTENT(in)                :: eq_comm

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
      INTEGER                            :: eq_rank
      INTEGER                            :: nstep_save
      REAL (rprec)                       :: start_time

!  local parameters
      REAL (rprec), PARAMETER            :: delt_factor = 2.0
      INTEGER, PARAMETER                 :: niter_factor = 4

!  Start of executable code
      start_time = profiler_get_start_time()

      delt_local = delt
      niter_local = niter
      num_iter_init = iterc
      recovery_index = 0

      eq_rank = 0
#if defined(MPI_OPT)
      CALL MPI_COMM_RANK(eq_comm, eq_rank, error)
#endif

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

         nstep_save = nstep
         CALL runvmec(ictrl_array, this%vmec_file_name, vmec_lscreen,               &
     &                eq_comm)
         nstep = nstep_save

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
                  IF (eq_rank .eq. 0) THEN
                     WRITE (*,1000) 'jac75_flag'
                     WRITE (*,1001) 'delt', delt, delt/delt_factor
                     WRITE (iou,1000) 'jac75_flag'
                     WRITE (iou,1001) 'delt', delt, delt/delt_factor
                     delt = delt/delt_factor
                  END IF
               ELSE
                  EXIT
               END IF

            CASE (more_iter_flag)

!  Restore the context to the state it was before runvmec was called. Disable
!  this on the inital convergence. Since the niter and niter_array variables in
!  the VMEC namelist input don't change, this has no effect on the initial
!  equilibrium convergence as a consequence of seting the restart_flag.
               CALL vmec_context_set_context(temp_context)

               IF (recovery_index .le. 2 .and. num_iter .gt. 1) THEN
                  IF (eq_rank .eq. 0) THEN
                     WRITE (*,1000) 'more_iter_flag'
                     WRITE (*,1002) 'niter', niter_local,                      &
     &                              niter_local + niter/niter_factor
                     WRITE (iou,1000) 'more_iter_flag'
                     WRITE (iou,1002) 'niter', niter_local,                    &
     &                                niter_local + niter/niter_factor
                  END IF
                  niter_local = niter_local + niter/niter_factor
               ELSE
                  EXIT
               END IF

            CASE DEFAULT
               EXIT

         END SELECT

         IF (eq_rank .eq. 0) THEN
            WRITE (*,*) 'Attempting error recovery'
            WRITE (iou,*) 'Attempting error recovery'
         END IF

         recovery_index = recovery_index + 1
      END DO

!  No longer need the local context array.
      CALL vmec_context_destruct(temp_context)

!  Check for errors, if vmec was successful, load vmec data for use with the
!  signals. Other wise report the error. This should only be performed on the
!  head node.
      IF (eq_rank .eq. 0) THEN
         vmec_converge = .false.
         SELECT CASE (ictrl_array(2))

            CASE (successful_term_flag, norm_term_flag)
               CALL read_wout_file(this%wout_file_name, error)
               CALL assert_eq(error, 0, "Error reading wout file.")
               CALL LoadRZL

               vmec_converge = .true.

            CASE DEFAULT
!>  @todo FIXME: Put in error reporting here. For now leave this blank until I
!>  can figure how v3fit decides if vmec has converged.

         END SELECT
      END IF

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
!>  @brief Loads the vacuum magnetic field file.
!>
!>  Loads the vacuum magnetic field file. This will get on multiple processes to
!>  allow parallel loading of the mgrid file. The extcur array will need to be
!>  broadcast to the child processes.
!>
!>  @param[in] this    A @ref vmec_class instance.
!>  @param[in] index   Index of the changed current.
!>  @param[in] eq_comm MPI communicator pool for VMEC.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_read_vac_file(this, index, eq_comm)
      USE vmec_input, only: extcur, mgrid_file, nfp
      USE mgrid_mod, only: read_mgrid, mgrid_path_old
      USE vacmod0, only: nv

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: index
      INTEGER, INTENT(in)           :: eq_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                       :: error

!  Start of executable code
      mgrid_path_old = ''

      CALL MPI_BCAST(extcur(index), 1, MPI_REAL8, 0, eq_comm, error)
      CALL read_mgrid(mgrid_file, extcur, nv, nfp, .false., error,             &
     &                comm = eq_comm)
      CALL assert_eq(0, error, 'vmec_set_param: failed to read ' //            &
     &               'mgrid_file')
#endif
      END SUBROUTINE

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
      INTEGER                          :: i
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (ASSOCIATED(this%vmec_context_save)) THEN
         CALL vmec_context_get_context(this%vmec_context_save)
      ELSE
         this%vmec_context_save => vmec_context_construct()
      END IF

!  Cache the wout file by appending _cache to a copy. Use copy here to keep the
!  orginal file intact.
      CALL copy_file(TRIM(this%wout_file_name),                                &
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
      INTEGER                          :: i
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (ASSOCIATED(this%vmec_context_save)) THEN
         CALL vmec_context_set_context(this%vmec_context_save)
      END IF

      CALL read_wout_file(this%wout_file_name, error)
      CALL assert_eq(error, 0, 'Error reading reset wout file.')
      CALL LoadRZL

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
      
!  Update the namelist input variables from the equilibrium solution.
      CALL vmec_set_namelist(this)

      iou_nl = 0
      CALL safe_open(iou_nl, status,                                           &
     &               TRIM(this%vmec_file_name) // '_out',                      &
     &               'replace', 'formatted', delim_in='quote')
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
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: current_step

!  local variables
      INTEGER                       :: iou_nl
      INTEGER                       :: status
      CHARACTER (len=path_length)   :: filename
      REAL (rprec)                  :: start_time

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
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize the vmec state to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref vmec_class instance.
!>  @param[in]    recon_comm MPI communicator for the reconstruction processes.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_sync_state(this, recon_comm)
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
       TYPE (vmec_class), INTENT(inout) :: this
       INTEGER, INTENT(in)              :: recon_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                           :: i
      INTEGER                           :: error
      INTEGER                           :: mpi_rank

!  Start of executable code
      CALL MPI_COMM_RANK(recon_comm, mpi_rank, error)

!  The barrier ensures that the wout is ready to be copied from the parent
!  directory. Otherwise a child process can try to copy a wout file that does
!  not exist.
      CALL MPI_BARRIER(recon_comm, error)

!  If this is the child process, load the wout file.
      IF (mpi_rank .gt. 0) THEN
         CALL copy_file(build_path('..', this%wout_file_name),                 &
     &                  TRIM(this%wout_file_name), error)

         CALL read_wout_file(this%wout_file_name, error)
         CALL assert_eq(error, 0, 'Error reading synced wout file.')
         CALL LoadRZL

      END IF

      CALL vmec_context_sync_state(this%vmec_context_save, recon_comm)

#endif
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Syncronize a child vmec state to the parent.
!>
!>  Syncs data between a child and the parent process. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref vmec_class instance.
!>  @param[in]    index      Reconstruction rank to sync.
!>  @param[in]    recon_comm MPI communicator for the reconstruction processes.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_sync_child(this, index, recon_comm)
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE file_opts

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(in)              :: index
      INTEGER, INTENT(in)              :: recon_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                           :: error
      INTEGER                           :: mpi_rank

!  Start of executable code
      CALL MPI_COMM_RANK(recon_comm, mpi_rank, error)

!  If this is the parent process, load the wout file.
      IF (mpi_rank .eq. 0) THEN
         CALL copy_file(build_path(process_dir(index + 1),                     &
     &                             this%wout_file_name),                       &
     &                  TRIM(this%wout_file_name), error)

         CALL read_wout_file(this%wout_file_name, error)
         CALL assert_eq(error, 0, 'Error reading synced wout file.')
         CALL LoadRZL

      END IF

      CALL vmec_context_sync_child(this%vmec_context_save, index,              &
     &                             recon_comm)

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
         zbcc => zmn_bdy(:,:,zcc)
         IF (lthreed) THEN
            rbcs => rmn_bdy(:,:,rcs)
            zbss => zmn_bdy(:,:,zss)
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
