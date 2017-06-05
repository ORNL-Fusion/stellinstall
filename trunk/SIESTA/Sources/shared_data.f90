!>  \brief Module containing data that is shared with other modules
      MODULE shared_data
      USE stel_kinds
      USE stel_constants

      INTEGER, PARAMETER :: ngmres_steps = 100      !<max number of gmres steps (10-100) should scale with ns?
      INTEGER, PARAMETER :: PREDIAG=1               !<diagonal preconditioning
      INTEGER, PARAMETER :: PREBLOCK=2              !<block preconditioning
      INTEGER            :: neqs                    !<No. elements in xc array
      INTEGER            :: ndims                   !<No. independent (real) variables, =3 (stell sym) or 6
      INTEGER            :: mblk_size               !<mpol X ntor X ndims linear block size in Fourier space
      INTEGER            :: nprecon
      INTEGER            :: nprecon_type
      INTEGER            :: niter
      INTEGER            :: ngmres_type = 2         !<=1, no "peek" improvement, =2 "peek" performance improvement
      INTEGER            :: iortho = 3              !controls orthogonalization in GMRES
                                                    !3: ICGS (recommended); 2: CGS; 1: IMGS; 0: MGS
      INTEGER            :: hessPass_Test=-1        !for dumping block and data files (testing)
      REAL(dp)           :: mupar_test

      REAL(dp),PARAMETER   :: fsq_res = 1.E-16_dp   !<threshold force for turning off resistive perturbations
      REAL(dp), PARAMETER  :: levm_ped=1.E-08_dp, mu_ped=1.E-6_dp  !<pedestal values of levenberg/mu|| (levm_ped between 10^-5 and 10^-10)
      REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET :: gc  !<1D Array of Fourier mode MHD force components
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE     :: buv_res  !<Resonant magnetic field perturbation
      REAL(dp), DIMENSION(:), ALLOCATABLE         :: xc  !<1D array of Fourier mode displacement components
      REAL(dp), DIMENSION(:), ALLOCATABLE         :: xc0, gc0
      REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE   :: col_scale  !<column scaling factor

      REAL(dp)        :: fsq_total, fsq_total1, fsq_gmres, fsq_lin, etak_tol
      REAL(dp)        :: levm_scale = 1
      REAL(dp)        :: wtotal, wtotal0, xmax_lin
      REAL(dp)        :: delta_t
      REAL(dp)        :: fsqvs, fsqvu, fsqvv
      REAL(dp)        :: ste(4)                    !<Spectral Truncation RMS Error    
      REAL(dp)        :: bs0(12), bu0(12)
      REAL(dp)        :: bsbu_ratio_s, jsju_ratio_s
      REAL(dp)        :: bsbu_ratio_a, jsju_ratio_a
      REAL(dp)        :: scale_s, scale_u          !Remove this once col_scale implemented

      
!Logicals for evolving origin and edge
      LOGICAL         :: l_push_edge                !<=FALSE, u,v components of forces NOT solved at s=1
      LOGICAL         :: l_push_s, l_push_u, l_push_v !<=FALSE, s,u,v components of forces NOT solved at s=0
      LOGICAL         :: l_natural=.TRUE.           !<=TRUE use natural bc at s=0 (preferred: maintains symmetry)
      LOGICAL         :: l_linearize, l_conjgrad
      LOGICAL         :: l_getwmhd                  !<=TRUE to compute MHD energy
      LOGICAL         :: l_getfsq                   !<=TRUE to compute integrated MHD force
      LOGICAL         :: l_ApplyPrecon              !<controls when preconditioner is applied
      LOGICAL         :: l_PrintOriginForces=.FALSE.
      LOGICAL         :: l_init_state               !M=<Internal flag to storing initial field/pressure state
      LOGICAL         :: l_update_state=.FALSE.     !<=TRUE before funct_island called to get ste array
      LOGICAL         :: l_par_state                !<=TRUE if parallel allocation state (in init_state_par)
      LOGICAL         :: lcolscale                  !<=TRUE to apply col-scaling to Hessian
      LOGICAL         :: lasym = .FALSE.            !<=FALSE, symmetric; =TRUE, asymmetric plasma
      LOGICAL         :: lrecon = .FALSE.

      END MODULE shared_data
