!*******************************************************************************
!  File vmec_equilibrium.f
!  Contains module vmec_equilibrium
!  Defines the base class of the type vmec:
!
!  Super Class: equilibrium
!  Sub Classes: None
!
!*******************************************************************************

      MODULE vmec_equilibrium
      USE stel_kinds
      USE data_parameters
      USE pprofile_T
      USE vparams, ONLY: ndatafmax

!*******************************************************************************
!  vmec equilibrium module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: vmec_profile_type_length = 30

!  Define id's only for values that can change. All others shouldn't be needed
!  outside of the vmec interface.

!  Variable Parameters
      INTEGER, PARAMETER :: vmec_ac_id            =  1
      INTEGER, PARAMETER :: vmec_ac_aux_s_id      =  2
      INTEGER, PARAMETER :: vmec_ac_aux_f_file_id =  3
      INTEGER, PARAMETER :: vmec_ai_id            =  4
      INTEGER, PARAMETER :: vmec_ai_aux_s_id      =  5
      INTEGER, PARAMETER :: vmec_ai_aux_f_id      =  6
      INTEGER, PARAMETER :: vmec_am_id            =  7
      INTEGER, PARAMETER :: vmec_am_aux_s_id      =  8
      INTEGER, PARAMETER :: vmec_am_aux_f_id      =  9
      INTEGER, PARAMETER :: vmec_bloat_id         = 10
      INTEGER, PARAMETER :: vmec_rbc_id           = 11
      INTEGER, PARAMETER :: vmec_zbs_id           = 12
      INTEGER, PARAMETER :: vmec_rbs_id           = 13
      INTEGER, PARAMETER :: vmec_zbc_id           = 14
      INTEGER, PARAMETER :: vmec_extcur_id        = 15
      INTEGER, PARAMETER :: vmec_curtor_id        = 16
      INTEGER, PARAMETER :: vmec_phiedge_id       = 17
      INTEGER, PARAMETER :: vmec_pres_scale_id    = 18

!  Auxilary Parameters
      INTEGER, PARAMETER :: vmec_rmnc_id          = 19
      INTEGER, PARAMETER :: vmec_zmns_id          = 20
      INTEGER, PARAMETER :: vmec_lmns_id          = 21
      INTEGER, PARAMETER :: vmec_gmnc_id          = 22
      INTEGER, PARAMETER :: vmec_bsubumnc_id      = 23
      INTEGER, PARAMETER :: vmec_bsubvmnc_id      = 24
      INTEGER, PARAMETER :: vmec_bsupsmns_id      = 25
      INTEGER, PARAMETER :: vmec_bsupumnc_id      = 26
      INTEGER, PARAMETER :: vmec_bsupvmnc_id      = 27
      INTEGER, PARAMETER :: vmec_rmns_id          = 28
      INTEGER, PARAMETER :: vmec_zmnc_id          = 29
      INTEGER, PARAMETER :: vmec_lmnc_id          = 30
      INTEGER, PARAMETER :: vmec_gmns_id          = 31
      INTEGER, PARAMETER :: vmec_bsubumns_id      = 32
      INTEGER, PARAMETER :: vmec_bsubvmns_id      = 33
      INTEGER, PARAMETER :: vmec_bsupsmnc_id      = 34
      INTEGER, PARAMETER :: vmec_bsupumns_id      = 35
      INTEGER, PARAMETER :: vmec_bsupvmns_id      = 36
      INTEGER, PARAMETER :: vmec_phi_id           = 37
      INTEGER, PARAMETER :: vmec_iotaf_id         = 38
      INTEGER, PARAMETER :: vmec_rzl_id           = 39

!  Auxilary Parameters 2
      INTEGER, PARAMETER :: vmec_bsubu_id         = 40
      INTEGER, PARAMETER :: vmec_bsubv_id         = 41
      INTEGER, PARAMETER :: vmec_rsuv_id          = 42
      INTEGER, PARAMETER :: vmec_zsuv_id          = 43
      INTEGER, PARAMETER :: vmec_gsuv_id          = 44
      INTEGER, PARAMETER :: vmec_currusuv_id      = 45
      INTEGER, PARAMETER :: vmec_currvsuv_id      = 46
      INTEGER, PARAMETER :: vmec_rusuv_id         = 47
      INTEGER, PARAMETER :: vmec_zusuv_id         = 48
      INTEGER, PARAMETER :: vmec_rvsuv_id         = 49
      INTEGER, PARAMETER :: vmec_zvsuv_id         = 50
      INTEGER, PARAMETER :: vmec_phitot_id        = 51

!  Auxilary Model Parameters
      INTEGER, PARAMETER :: vmec_pp_ne_b_id       = 52
      INTEGER, PARAMETER :: vmec_pp_ne_as_id      = 53
      INTEGER, PARAMETER :: vmec_pp_ne_af_id      = 54
      INTEGER, PARAMETER :: vmec_pp_sxrem_b_id    = 55
      INTEGER, PARAMETER :: vmec_pp_sxrem_as_id   = 56
      INTEGER, PARAMETER :: vmec_pp_sxrem_af_id   = 57
      INTEGER, PARAMETER :: vmec_pp_te_b_id       = 58
      INTEGER, PARAMETER :: vmec_pp_te_as_id      = 59
      INTEGER, PARAMETER :: vmec_pp_te_af_id      = 60

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) vmec base class
!  2) vmec auxilary model parameters
!  3) vmec fixed parameters
!  4) vmec variable parameters
!  5) vmec aux 1 parameters
!  6) vmec aux 2 parameters
!
!*******************************************************************************

      TYPE vmec_class
         TYPE (vmec_aux_model)           :: aux_model
         TYPE (vmec_fixed_parameters)    :: fixed_parameters
         TYPE (vmec_variable_parameters) :: aux_model
         TYPE (vmec_aux_1_parameters)    :: aux_1_parameters
         TYPE (vmec_aux_2_parameters)    :: aux_2_parameters
      END TYPE

!-------------------------------------------------------------------------------
!  Define the vmec auxilary model parameters type
!
!  ne      Parameterized profile describing a flux surface constant electron
!          denisty profile.
!
!  te      Parameterized profile describing a flux surface constant electron
!          temperature profile.
!
!  sxrem   Parameterized profile describing a flux surface conatsnt soft x-ray
!          emissivity profile.
!-------------------------------------------------------------------------------
      TYPE vmec_aux_model
         TYPE (pprofile), POINTER :: ne => null()
         TYPE (pprofile), POINTER :: te => null()
         TYPE (pprofile), POINTER :: sxrem => null()
      END TYPE

!-------------------------------------------------------------------------------
!  Define the vmec fixed parameters type
!
!  lasym      logical,   .true.  stellarator asymmetric
!                        .false. stellarator symmetric
!                        vmec_input : vmec_input.f
!
!  lfreeb     logical,   .true.  Do fixed boundary equilibrium
!                        .false. Do free boundary equilibrium
!                        vmec_input : vmec_input.f
!
!  lforbal    logical,   .true.  improved avg force bal
!                        .false. variational avg force bal
!                        vmec_input : vmec_input.f
!
!  lRFP       logical,   .true.  Reversed Field Pinch
!                        .false. Tokamaks, stellarators
!                        vmec_input : vmec_input.f
!
!  mgrid_file character, name for mgrid file.
!                        vmec_input : vmec_input.f
!
!  nfp        integer,   number of toroidal field periods (=1 for Tokamak)
!                        vmec_input : vmec_input.f
!
!  ncurr      integer,   flux conserving (=0) or prescribed toroidal current
!                        (=1)
!                        vmec_input : vmec_input.f
!
!  nvacskip   integer,   iterations skipped between full updates of vacuum
!                        solution
!                        vmec_input : vmec_input.f
!
!  mpol       integer,   number of poloidal modes
!                        vmec_input : vmec_input.f
!
!  ntor       integer,   number of toroidal modes
!                        vmec_input : vmec_input.f
!
!  ntheta     integer,   number of poloidal points, 0 - 2Pi
!                        vmec_input : vmec_input.f
!
!  nzeta      integer,   number of toroifal points per field period
!                        vmec_input : vmec_input.f
!
!  mns        integer,   length of Xmn arrays
!                        vmec_dim : vmec_dim.f
!
!  ns_array   integer,   array of radial (multi-) grids sizes
!                        vmec_dim : vmec_dim.f
!
!  nextcur    integer,   number of external current to use (extcur array)
!                        mgrid_mod.f
!
!  pcurr_type character, specify type of current profile
!                        vmec_input : vmec_input.f
!
!  piota_type character, specify type of iota or q profile
!                        vmec_input : vmec_input.f
!
!  pmass_type character, specify type of pressure profile
!                        vmec_input : vmec_input.f
!
!- Run Control Parameters ------------------------------------------------------
!
!  niter      integer,   maximum number of iterations for convergence (VMEC will
!                        retry once if ftol convergence level not achieved yet)
!                        vmec_dim : vmec_dim.f
!
!  nstep      integer,   iterations between printout of equilibrium data to
!                        threed1 file/screen
!                        vmec_dim : vmec_dim.f
!
!  delt       real,      intial guess (~1) for time step
!                        vmec_dim : vmec_dim.f
!
!  ftol_array real,      array of force convergence tolerances on radial grids
!                        defined by ns_array
!                        vmec_dim : vmec_dim.f
!
!  gamma      real,      value of compressibility index (gamma=0 => pressure
!                        prescribed)
!                        vmec_input : vmec_input.f
!
!  tcon0      real,      optional parameter weight for spectral constraint force
!                        (~1)
!                        vmec_dim : vmec_dim.f
!-------------------------------------------------------------------------------
      TYPE vmec_fixed_parameters
         LOGICAL                                  :: lasym
         LOGICAL                                  :: lfreeb
         LOGICAL                                  :: lforbal
         LOGICAL                                  :: lRFP
         CHARACTER (len=data_filename_length)     :: mgrid_file
         INTEGER                                  :: nfp
         INTEGER                                  :: ncurr
         INTEGER                                  :: nvacskip
         INTEGER                                  :: mpol
         INTEGER                                  :: ntor
         INTEGER                                  :: ntheta
         INTEGER                                  :: nzeta
         INTEGER                                  :: mns
         INTEGER, DIMENSION(100)                  :: ns_array
         INTEGER                                  :: nextcur
         CHARACTER (len=vmec_profile_type_length) :: pcurr_type
         CHARACTER (len=vmec_profile_type_length) :: piota_type
         CHARACTER (len=vmec_profile_type_length) :: pmass_type
         INTEGER                                  :: niter
         INTEGER                                  :: nstep
         REAL (rprec)                             :: delt
         REAL (rprec), DIMENSION(100)             :: ftol_array
         REAL (rprec)                             :: gamma
         REAL (rprec)                             :: tcon0
      END TYPE

!-------------------------------------------------------------------------------
!  Define the vmec variable parameters type
!
!  ns_index   integer,   Index to the ns_array and ftol_array arrays. Called
!                        multi_ns_grid(vmec_main.f) in readin.f and eqsolve.f
!
!  rbc_bounds integer,   Bounds od the rbc, rbs, zbc, and zbs arrays order
!                        assumed is lower_1, upper_1, lower_2, upper_2
!
!  ac         real,      array(0:20) surface averaged toroidal current density
!                        (ncurr=1) expansion coefficients (power series in s)
!                        vmec_input : vmec_input.f
!
!  ac_aux_s   real,      array(ndatafmax[101]) Auxiliary 's' array for ac
!                        specification (splines)
!                        vmec_input : vmec_input.f
!
!  ac_aux_f   real,      array(ndatafmax[101]) Auxiliary 'f' array for ac
!                        specification (splines)
!                        vmec_input : vmec_input.f
!
!  ai         real,      array(0:20) surface dimensionless iota expansion
!                        coefficients (power series in s)
!                        vmec_input : vmec_input.f
!
!  ai_aux_s   real,      array(ndatafmax[101]) Auxiliary 's' array for ai
!                        specification (splines)
!                        vmec_input : vmec_input.f
!
!  ai_aux_f   real,      array(ndatafmax[101]) Auxiliary 'f' array for ai
!                        specification (splines)
!                        vmec_input : vmec_input.f
!
!  am         real,      array(0:20) mass or pressure expansion coefficients
!                        (power series in s) in MKS units (NWT/M**2)
!                        vmec_input : vmec_input.f
!
!  am_aux_s   real,      array(ndatafmax[101]) Auxiliary 's' array for am
!                        specification (splines)
!                        vmec_input : vmec_input.f
!
!  am_aux_f   real,      array(ndatafmax[101]) Auxiliary 'f' array for am
!                        specification (splines)
!                        vmec_input : vmec_input.f
!
!  bloat      real,      used to expand the computational region
!                        vmec_input : vmec_input.f
!
!  rbc        real,      array, (-ntord:ntord,0:mpol1d)
!                        boundary coefficients of COS(m*theta-n*zeta) for R
!                        vmec_input : vmec_input.f
!                        ALLOCATABLE
!
!  zbs        real,      array, (-ntord:ntord,0:mpol1d)
!                        boundary coefficients of SIN(m*theta-n*zeta) for Z
!                        vmec_input : vmec_input.f
!                        ALLOCATABLE
!
!  rbs        real,      array, (-ntord:ntord,0:mpol1d)
!                        boundary coefficients of SIN(m*theta-n*zeta) for R
!                        vmec_input : vmec_input.f
!                        ALLOCATABLE
!
!  zbc        real,      array, (-ntord:ntord,0:mpol1d)
!                        boundary coefficients of COS(m*theta-n*zeta) for Z
!                        vmec_input : vmec_input.f
!                        ALLOCATABLE
!
!  extcur     real,      array(nigroup[nigroup set to 100 in vsvd0.f])
!                        array of currents in each external current group. Used
!                        to multiply Green''s function for fields and loops read
!                        in from MGRID file. Should use real current units (A).
!                        vmec_input : vmec_input.f
!
!  curtor     real,      net toroidal current (used if ncurr=1) in MKS [A]
!                        vmec_dim : vmec_dim.f
!
!  phiedge    real,      value of toroidal flux at plasma edge (defines plasma
!                        volume in free-bdy case)
!                        vmec_input : vmec_input.f
!
!  pres_scale real,      factor used to scale pressure profile (default value =
!                        1) useful so user can fix profile and change beta
!                        without having to change all AM coefficients separately
!                        vmec_input : vmec_input.f
!-------------------------------------------------------------------------------
      TYPE vmec_variable_parameters
         INTEGER                              :: ns_index
         INTEGER                              :: rbc_bounds(4)
         REAL(rprec), DIMENSION(0:20)         :: ac
         REAL(rprec), DIMENSION(ndatafmax)    :: ac_aux_s
         REAL(rprec), DIMENSION(ndatafmax)    :: ac_aux_f
         REAL(rprec), DIMENSION(0:20)         :: ai
         REAL(rprec), DIMENSION(ndatafmax)    :: ai_aux_s
         REAL(rprec), DIMENSION(ndatafmax)    :: ai_aux_f
         REAL(rprec), DIMENSION(0:20)         :: am
         REAL(rprec), DIMENSION(ndatafmax)    :: am_aux_s
         REAL(rprec), DIMENSION(ndatafmax)    :: am_aux_f
         REAL(rprec)                          :: bloat
         REAL(rprec), DIMENSION(:,:), POINTER :: rbc => null()
         REAL(rprec), DIMENSION(:,:), POINTER :: zbs => null()
         REAL(rprec), DIMENSION(:,:), POINTER :: rbs => null()
         REAL(rprec), DIMENSION(:,:), POINTER :: zbc => null()
         REAL(rprec), DIMENSION(100)          :: extcur
         REAL(rprec)                          :: curtor
         REAL(rprec)                          :: phiedge
         REAL(rprec)                          :: pres_scale
      END TYPE

      END MODULE