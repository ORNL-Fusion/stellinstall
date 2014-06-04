      MODULE parambs
      USE vmec0
      USE bootsj_input
      REAL(rprec), PARAMETER :: pi = 3.141592653589793238462643_dp
      REAL(rprec), PARAMETER :: mu0 = 4.0e-7_dp*pi
      INTEGER, PARAMETER :: nlambda = 101       !this now ONLY governs the trapped calcualtion
      INTEGER :: nthetah, nzetah                !poloidal, toridal grid parameters
      INTEGER, PARAMETER :: nthetahm = 32       !poloidal, toroidal  grid refinement
      INTEGER, PARAMETER :: nzetahm = 16        !for bmax
      INTEGER, PARAMETER :: iotasign = 1
      REAL(rprec), PARAMETER :: zetasign = -1
      INTEGER :: irdim, irup
      INTEGER, DIMENSION(:), ALLOCATABLE :: jlist, jlist_idx
      REAL(rprec)  delta_rho_1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: flux,
     1    aiogar, aipsi, gpsi, pres1,
     2    betar, dense, densi, tempe1, tempi1
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1  bfield, gsqrt_b,  b2obm, omb32, bfieldm,    !  arrays to store quantities on the
     2  sinmi, sinnj, cosmi, cosnj,                 !  theta phi grid
     2  sinmim, sinnjm, cosmim, cosnjm              !  quantities that are on the theata

      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1    dibs, aibs, dibst, aibst, phip,
     2    bsdense, bsdensi, bstempe, bstempi, bsdenste, bsdensti,
     3    bstempte, bstempti, qsafety, capr, caps, h2,
     4    ftrapped, fpassing, epsttok, fttok, b2avg,
     5    gbsnorm, aiterm1, other1, ajBbs,
     6    rhoar, bsnorm, fptok, amain, d_rho,
     7    bmax1, thetamax, zetahmax
      REAL(rprec)   alphae, alphai , psimax
      REAL(rprec)   temperho1, tempirho1, densrho1
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: amnfit
      REAL(rprec)   periods
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: theta, zetah
      REAL(rprec), DIMENSION(nthetahm) :: thetam
      REAL(rprec), DIMENSION(nzetahm) :: zetahm
      complex(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    dmn, fmn, alpha1mn
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rfmn
      REAL(rprec)   avgbobm2, sum_gsqrt_b
      REAL(rprec), DIMENSION(136) :: w1, alamb_h
CCCCfix needed--put nlambda variables into the trapped fraction calculation
      REAL(rprec)   dlambda, drho
      REAL(rprec)  sign_jacobian
      LOGICAL, DIMENSION(:), ALLOCATABLE :: lsurf
      LOGICAL l_boot_all, lscreen
      END MODULE parambs
