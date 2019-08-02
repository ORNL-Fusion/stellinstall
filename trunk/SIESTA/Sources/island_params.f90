      MODULE island_params
      USE stel_kinds
!
!     CONTAINS VARIABLE DECLARATIONS FOR (ISLAND) CODE
!     THAT WILL BE SHARED THROUGHOUT THE PROJECT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
	  INTEGER :: ns_i, nu_i, nv_i, nuv_i                                  ! size of s, u, v dimensions of (island) metric tensors
      INTEGER :: mpol_i, ntor_i
      INTEGER :: ns_v, mpol_v, ntor_v
      INTEGER :: nfp_i, mnmax_i
      INTEGER :: nsh                                                      ! Number of points in half mesh
      INTEGER, ALLOCATABLE :: irefu(:), jrefv(:)                          ! -u, -v indices
!>  Radial grid spacing. hs = s_i+1 - s_i
      REAL (dp) :: hs_i
!>  Radial grid derivative factor.
!>    d X/ ds = (x_i+1 - x_i)/(s_i+1 - s_i)                                  (1)
!>  Where ohs = 1/(s_i+1 - s_i) = ns - 1
      REAL (dp) :: ohs_i
!>  Differential area element.
!>    dudv 1/(ntheta*nzeta)                                                  (1)
      REAL (dp) :: dnorm_i
!>  
      REAL (dp) :: gnorm_i
!>  Magnetic energy. 2Pi*wb @ref read_wout_mod::wb.
      REAL (dp) :: wb_i
!>  Thermal energy. 2Pi*wp @ref read_wout_mod::wp.
      REAL (dp) :: wp_i
!>  Equilibrium volume @ref read_wout_mod::volume.
      REAL (dp) :: volume_i
!>  Major radius.
      REAL (dp) :: rmajor_i
      REAL(dp), PARAMETER :: gamma = 5._dp/3._dp                          ! Adiabatic constant
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: cosmu, cosmum,            &
          cosmui, sinmu, sinmum, sinmui, cosnv, cosnvn, sinnv, sinnvn
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: phipf_i, chipf_i,         &
                                               presf_i, vp_f
!-----------------------------------------------

      END MODULE island_params
