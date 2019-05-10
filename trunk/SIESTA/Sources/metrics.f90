!*******************************************************************************
!>  @file metrics.f90
!>  @brief Contains module @ref metrics.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module for reading in VMEC input and computing metric elements on the SIESTA
!>  (sqrt-flux) radial grid. Computes and stores the real-space metric elements,
!>  Jacobian based on square root of fluc splined vmec coordinate system.
!*******************************************************************************

      MODULE metrics
      USE stel_kinds
      USE stel_constants
      USE island_params
      USE shared_data, ONLY: lasym, r1_i, z1_i, ru_i, zu_i, rv_i, zv_i,  &
                             lverbose, jsupvdotA, nsp
      USE read_wout_mod, ns_vmec=>ns, mpol_vmec=>mpol, ntor_vmec=>ntor,  &
          rmnc_vmec=>rmnc, zmns_vmec=>zmns, lmns_vmec=>lmns,             &
          xm_vmec=>xm, xn_vmec=>xn, chipf_vmec=>chipf,                   &  ! MRC 4/1/2016
          rmns_vmec=>rmns, zmnc_vmec=>zmnc, lmnc_vmec=>lmnc,             &  ! MRC 12/1/2016
          phipf_vmec=>phipf, presf_vmec=>presf, nfp_vmec=>nfp,           &
          wb_vmec=>wb, wp_vmec=>wp, gamma_vmec=>gamma,                   &
          volume_vmec=>volume, raxis_vmec=>raxis, lasym_vmec=>lasym,     &
          iasym_vmec=>iasym
      USE descriptor_mod, ONLY: iam
      USE grid_extension, ONLY: grid_extender, read_vessel_file
      USE timer_mod
      USE utilities, ONLY: to_full_mesh
      USE siesta_namelist, ONLY: l_vessel
      USE vmec_info

      IMPLICIT NONE

!*******************************************************************************
!  metrics module parameters
!*******************************************************************************
!>  Flag to switch between using Sqrt(Phi) or Phi directly.
      LOGICAL, PARAMETER, PRIVATE :: iScale = .true.
!>  Symmetric parity flag.
      INTEGER, PARAMETER, PRIVATE :: isym = 0
!>  Asymmetric parity flag.
      INTEGER, PARAMETER, PRIVATE :: iasym = 1
!>  Spline lower boundary condition. Uses last four points.
      REAL(dp), PARAMETER         :: lower_b = -1.e30_dp
!>  Spline upper boundary condition. Uses last four points.
      REAL(dp), PARAMETER         :: upper_b = -1.e30_dp

!>     SET iflipj -> -1 FOR POSITIVE JACOBIAN SYSTEM
      INTEGER, PARAMETER          :: iflipj = 1

!*******************************************************************************
!  metrics module variables
!*******************************************************************************
!>  Jacobian on half grid
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: sqrtg
!>  Lower metric tensor half grid. e_s . e_s
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gss
!>  Lower metric tensor half grid. e_s . e_u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gsu
!>  Lower metric tensor half grid. e_s . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gsv
!>  Lower metric tensor half grid. e_u . e_u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: guu
!>  Lower metric tensor half grid. e_u . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: guv
!>  Lower metric tensor half grid. e_v . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gvv
!>  Upper metric tensor half grid. e^s . e^s
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: hss
!>  Upper metric tensor. e^s . e^u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: hsu
!>  Upper metric tensor. e^s . e^v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: hsv
!>  Upper metric tensor. e^u . e^u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: huu
!>  Upper metric tensor. e^u . e^v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: huv
!>  Upper metric tensor. e^v . e^v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: hvv
!>  Lower metric tensor full grid. e_s . e_s
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gssf
!>  Lower metric tensor full grid. e_s . e_u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gsuf
!>  Lower metric tensor full grid. e_s . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gsvf
!>  Lower metric tensor full grid. e_u . e_u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: guuf
!>  Lower metric tensor full grid. e_u . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: guvf
!>  Lower metric tensor full grid. e_v . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gvvf

!>  Maximum of the grid R inside the last closed flux surface.
      REAL(dp)                              :: rmax
!>  Minimum of the grid R inside the last closed flux surface.
      REAL(dp)                              :: rmin
!>  Maximum of the grid Z inside the last closed flux surface.
      REAL(dp)                              :: zmax
!>  Minimum of the grid Z inside the last closed flux surface.
      REAL(dp)                              :: zmin

!>  Start timer
      REAL(dp)                              :: skston
!>  Stop timer
      REAL(dp)                              :: skstoff

!*******************************************************************************
!  Interface Blocks
!*******************************************************************************
      PRIVATE spline_fourier_modes, loadrzl_vmec, add_ghost_points

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Initialize the metric elements.
!>
!>  Loads values on to the SIESTA meshes. This involves splining and
!>  interpolating the VMEC quantities from the VMEC radial grid to the SIESTA
!>  radial grid.
!>
!>  @param[in] ns_in     Number of SIESTA radial surfaces.
!>  @param[in] mpol_in   Number of SIESTA poloidal modes.
!>  @param[in] ntor_in   Number of SIESTA toroidal modes.
!>  @param[in] nfp_in    Number of SIESTA field periods.
!>  @param[in] wout_file VMEC woutfile.
!-------------------------------------------------------------------------------
      SUBROUTINE init_metric_elements(ns_in, mpol_in, ntor_in,           &
                                      nfp_in, wout_file)
      USE timer_mod, ONLY: init_timers
      USE shared_data, ONLY: torflux, polflux
      USE island_params, hs=>hs_i, ns=>ns_i
      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(IN)          :: ns_in
      INTEGER, INTENT(IN)          :: mpol_in
      INTEGER, INTENT(IN)          :: ntor_in
      INTEGER, INTENT(IN)          :: nfp_in
      CHARACTER(len=*), INTENT(in) :: wout_file

!  Local variables
      INTEGER                      :: istat
      INTEGER                      :: js
      REAL(dp)                     :: temp

!  Start of executable code

#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL init_timers
#if defined(SKS) 
      CALL second0(skstoff)
      init_timers_time=init_timers_time+(skstoff-skston)
#endif

      IF (.not.iScale .AND. IAM.EQ.0)                                   &
        PRINT *,'ISCALE=False using VMEC grid.'

!  READ DATA FROM VMEC WOUT FILE (LIBSTELL ROUTINE)
#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL read_wout_file(wout_file, istat)
#if defined(SKS) 
      CALL second0(skstoff)
      read_wout_file_time=read_wout_file_time+(skstoff-skston)
#endif
      CALL ASSERT(istat.EQ.0, 'Read-wout error in INIT_METRIC_ELEMENTS')

      IF (nfp_in .lt. 1) THEN
         nfp_i = nfp_vmec
      ELSE
         nfp_i = nfp_in
      END IF
      wb_i  = (twopi*twopi)*wb_vmec
      wp_i  = (twopi*twopi)*wp_vmec
      volume_i = volume_vmec

      CALL ASSERT(nfp_i.NE.0.AND.ALLOCATED(raxis_vmec),    &
         'ERROR INITIALIZING FROM WOUT_FILE IN METRICS!')

      rmajor_i = raxis_vmec(0,1)

      CALL ASSERT(wb_i.gt.zero,'wb_vmec = 0!')
      gnorm_i = wb_i + wp_i/(gamma-1)
      gnorm_i = 1._dp/gnorm_i

      ns_i = ns_in
      nsh  = ns_i-1
      mpol_i = mpol_in
      ntor_i = ntor_in

!  Set number of points == number of modes for now! (May want mid-points for
!  flux conservation)
      nu_i = mpol_i+2  
      nv_i = 2*ntor_i+2

!USE 3/2 (ORSZAG) RULE FOR ANTI-ALIASING OF EVOLUTION EQNS
!Suppresses RADIAL grid separation in pressure
      nu_i = 3*nu_i/2
!SPH051617      nu_i = nu_i+MOD(nu_i,2)
      IF (lasym) nu_i = 2*nu_i

!SPH051617
      IF (.NOT.lasym .AND. MOD(nu_i,2).NE.1) nu_i = nu_i+1

      nv_i = 3*nv_i/2
      nv_i = nv_i+MOD(nv_i,2)

      IF (ntor_i .EQ. 0) nv_i = 1
      nuv_i = nu_i*nv_i
      mnmax_i = (mpol_i + 1)*(2*ntor_i + 1) ! Added RS. Contains total number of modes.

!  READ VACUUM VESSEL DATA (before calling init_fourier) (points on the vessel
!  surface are already re-distributed by VMESH code).
      IF (l_vessel) THEN
         CALL read_vessel_file(wout_file, istat)
         IF (istat .EQ. 0 .and. iam .EQ. 0) PRINT *,'Vacuum Vessel data successfully read.'
      END IF

!  UPDATED WOUT TO WRITE OUT CHIPF (SPH)
!
!    iotaf_vmec = iotaf_vmec*phipf_vmec                 MRC 4/1/2016
!
!  SPLINE R, Z. L FOURIER COMPONENTS in s FROM ORIGINAL VMEC MESH (s ~ phi,
!  ns_vmec points) TO A "POLAR" MESH [s ~ sqrt(phi), ns_i POINTS] WITH BETTER
!  AXIS RESOLUTION
      CALL vmec_info_construct_spline(mnmax, ns_i, lasym)
      CALL GetRZLSplines(rmnc_vmec, zmns_vmec, lmns_vmec,               &
                         currumnc, currvmnc,                            &
                         rmnc_spline, zmns_spline, lmns_spline,         &
                         currumnc_spline, currvmnc_spline, isym)
      IF (lasym)                                                        &
      CALL GetRZLSplines(rmns_vmec, zmnc_vmec, lmnc_vmec,               &
                         currumns, currvmns,                            &
                         rmns_spline, zmnc_spline, lmnc_spline,         &
                         currumns_spline, currvmns_spline, iasym)

!  Spline 1-D arrays: careful -> convert phipf VMEC and multiply by
!  ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island
      ALLOCATE(phipf_i(ns_i), chipf_i(ns_i), presf_i(ns_i), stat=istat)
      CALL Spline_OneD_Array (chipf_vmec, chipf_i, istat)
      CALL Spline_OneD_Array (phipf_vmec, phipf_i, istat)
      presf_vmec = mu0*presf_vmec
      CALL Spline_OneD_Array (presf_vmec, presf_i, istat)

!  Scale phipf_i and convert to sqrt(flux) mesh by multiplying by 2*s
      phipf_i = phipf_i/(2*pi)
      chipf_i = chipf_i/(2*pi)

!  Mapping s_vmec = (s_siesta)^2, so d(s_vmec)/d(s_siesta) = 2 s_siesta, where
!  s_siesta = hs_i*(js-1)
      IF (iScale) THEN
         DO js = 1, ns_i
            phipf_i(js) = 2*hs_i*(js - 1)*phipf_i(js)
            chipf_i(js) = 2*hs_i*(js - 1)*chipf_i(js)
         END DO
      END IF

      IF (iam .EQ. 0) THEN
         temp = twopi*hs_i*(SUM(phipf_i(2:ns_i)) + SUM(phipf_i(1:ns_i-1)))/2
         WRITE (33, 1000) volume_i, 1.E-6_dp*(wb_i+wp_i/(gamma-1))/mu0,  &
                          1.E-6_dp*wp_i/mu0, temp, wp_i/wb_i,            &
                          chipf_i(2)/phipf_i(2),                         &
                          chipf_i(ns_i)/phipf_i(ns_i), gamma
      END IF

! Added SKS
      ALLOCATE(torflux(ns_i), polflux(ns_i))
      torflux(1) = 0
      polflux(1) = 0;
      DO js = 2, ns_i
         polflux(js) = polflux(js - 1) + hs*(chipf_i(js) + chipf_i(js - 1))/2
         torflux(js) = torflux(js - 1) + hs*(phipf_i(js) + phipf_i(js - 1))/2
      END DO

!  CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON SQRT(FLUX) - "POLAR" - MESH AND
!  COMPUTE METRIC ELEMENTS AND JACOBIAN
      CALL LoadRZL_VMEC(istat)
      CALL ASSERT(istat.EQ.0,'LoadRZL error in INIT_METRIC_ELEMENTS')

      IF (l_vessel) THEN
        CALL grid_extender(wout_file, 'quad', istat)
      END IF

!  COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN
      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
!  COMPUTE FULL-MESH METRIC ELEMENTS
      CALL full_mesh_metrics

1000  FORMAT(/,' INITIAL PARAMETERS (FROM VMEC)',/,                      &
               ' PHYSICS QUANTITIES',/,                                  &
               ' PLASMA VOLUME (M^3): ',1pE12.4,                         &
               ' TOTAL MHD ENERGY (MJ):   ', 1pE16.8,/,                  &
               ' THERMAL ENERGY (MJ): ', 1pE12.4,                        &
               ' EDGE TOROIDAL FLUX (Wb): ', 1pE12.4,/,                  &
               ' DIMENSIONLESS QUANTITIES',/,                            &
               ' <BETA>: ',1pE12.4, ' IOTA(0): ',1pE12.4,                &
               ' IOTA(1) ',1pE12.4, ' GAMMA: ',1pE12.4,/, 21('-'),/)

      END SUBROUTINE init_metric_elements

!-------------------------------------------------------------------------------
!>  @brief Surface average a quantity.
!>
!>  @param[out] average Surface average
!>  @param[in]  q3d     3D quantity is real space.
!>  @param[in]  nsmin   Minimum radial index.
!>  @param[in]  nsmax   Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE SurfAvg(average, q3d, nsmin, nsmax)
      USE island_params, ONLY: cosmui

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nsmin:nsmax), INTENT(out)          :: average
      REAL(dp), DIMENSION(nu_i,nv_i,nsmin:nsmax), INTENT(IN) :: q3d
      INTEGER, INTENT(IN)                                    :: nsmin
      INTEGER, INTENT(IN)                                    :: nsmax

!  Local variables
      INTEGER                                                :: l
      INTEGER                                                :: js
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)                :: wint

!  Start of executable code

      ALLOCATE(wint(nu_i,nv_i,nsmin:nsmax))
      DO l = 1, nu_i
         wint(l,:,:) = cosmui(l,0)
      END DO
      DO js = nsmin, nsmax
         average(js) = SUM(q3d(:,:,js)*wint(:,:,js))
      END DO
      DEALLOCATE(wint)

      END SUBROUTINE SurfAvg
 
!-------------------------------------------------------------------------------
!>  @brief Splines the VMEC quantities.
!>
!>  Recasts the VMEC quantities to the SIESTA mesh.
!>
!>  @param[in]    rmn_vmec       VMEC R Fourier amplitudes.
!>  @param[in]    zmn_vmec       VMEC Z Fourier amplitudes.
!>  @param[in]    lmn_vmec       VMEC lambda Fourier amplitudes.
!>  @param[in]    currumn_vmec   VMEC J^u Fourier amplitudes.
!>  @param[in]    currvmn_vmec   VMEC J^v Fourier amplitudes.
!>  @param[in]    rmn_spline     Spline R Fourier amplitudes.
!>  @param[in]    zmn_spline     Spline Z Fourier amplitudes.
!>  @param[in]    lmn_spline     Spline lambda Fourier amplitudes.
!>  @param[in]    currumn_spline Spline J^u Fourier amplitudes.
!>  @param[in]    currvmn_spline Spline J^v Fourier amplitudes.
!>  @param[in]    parity         Parity of the Fourier amplitudes.
!-------------------------------------------------------------------------------
      SUBROUTINE GetRZLSplines(rmn_vmec, zmn_vmec, lmn_vmec,            &
                               currumn_vmec, currvmn_vmec,              &
                               rmn_spline, zmn_spline, lmn_spline,      &
                               currumn_spline, currvmn_spline, parity)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(:,:), INTENT(in)  :: rmn_vmec
      REAL(dp), DIMENSION(:,:), INTENT(in)  :: zmn_vmec
      REAL(dp), DIMENSION(:,:), INTENT(in)  :: lmn_vmec
      REAL(dp), DIMENSION(:,:), INTENT(in)  :: currumn_vmec
      REAL(dp), DIMENSION(:,:), INTENT(in)  :: currvmn_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rmn_spline
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: zmn_spline
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: lmn_spline
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: currumn_spline
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: currvmn_spline
      INTEGER, INTENT(in)                   :: parity

!  Local variables
      INTEGER                               :: istat
      INTEGER                               :: modes
      INTEGER                               :: modes_nyq
      INTEGER                               :: n_off
      INTEGER                               :: n_max
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tempmn1
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tempmn2

!  Start of executable code
      IF (parity.EQ.iasym .AND. .NOT.lasym_vmec) THEN
         rmn_spline = 0
         zmn_spline = 0
         lmn_spline = 0
         currumn_spline = 0
         currvmn_spline = 0
         RETURN
      END IF

      istat = 0
      CALL second0(skston)
      CALL Spline_Fourier_Modes(rmn_vmec, rmn_spline, xm_vmec, istat)
      CALL second0(skstoff)
      CALL ASSERT(istat.EQ.0,'Error splining rmn')

      CALL second0(skston)
      CALL Spline_Fourier_Modes(zmn_vmec, zmn_spline, xm_vmec, istat)
      CALL second0(skstoff)
      CALL ASSERT(istat.EQ.0,'Error splining zmn')

      ALLOCATE(tempmn1(mnmax, ns_vmec + 1), stat=istat)
      CALL Add_Ghost_Points(lmn_vmec, tempmn1)
      CALL second0(skston)
      CALL Spline_Fourier_Modes(tempmn1, lmn_spline, xm_vmec, istat)
      CALL second0(skstoff)
      CALL ASSERT(istat.EQ.0,'Error splining lmn')
      DEALLOCATE(tempmn1)

!  The currents use the _nyq spectrum so we need to truncae down to just the
!  non-nyquest modes.
      ALLOCATE (tempmn1(mnmax,ns_vmec), tempmn2(mnmax,ns_vmec), stat=istat)
      tempmn1 = 0
      tempmn2 = 0

!  Assuming here that _nyq modes will be larger than the regular modes.
      CALL ASSERT(mnmax_nyq .ge. mnmax, 'ERROR: More modes than _nyq modes')

      n_max = INT(MAXVAL(xn_vmec))
      n_off = 2*INT((MAXVAL(xn_nyq) - MAXVAL(xn_vmec))/nfp_vmec) + 1

      modes_nyq = 1
      DO modes = 1, mnmax
         tempmn1(modes,:) = mu0*currumn_vmec(modes_nyq,:)
         tempmn2(modes,:) = mu0*currvmn_vmec(modes_nyq,:)

         IF (INT(xn_vmec(modes)) .eq. n_max) THEN
            modes_nyq = modes_nyq + n_off
         ELSE
            modes_nyq = modes_nyq + 1
         END IF
      END DO
      IF (parity.EQ.isym) THEN
         jsupvdotA = mu0*(SUM(tempmn2(1,2:ns_vmec - 1))          &
                   + (tempmn2(1, 1) + tempmn2(1, ns_vmec))/2)
         jsupvdotA = jsupvdotA/(ns_vmec - 1)
      END IF

      CALL Spline_Fourier_Modes(tempmn1, currumn_spline, xm_vmec, istat)
      CALL Spline_Fourier_Modes(tempmn2, currvmn_spline, xm_vmec, istat)

      DEALLOCATE(tempmn1, tempmn2)
     
      END SUBROUTINE GetRZLSplines
 
!-------------------------------------------------------------------------------
!>  @brief Spline reinterpolate a VMEC Fourier quantity.
!>
!>  Redistributes the VMEC quantity from the VMEC radial grid to he SIESTA
!>  radial grid.
!>
!>  @param[in]    y_vmec   Orginal VMEC Fourier quantity.
!>  @param[out]   y_spline Splined VMEC Fourier quantity.
!>  @param[in]    xm_in    Poloidal mode numbers.
!>  @param[inout] istat    Error status.
!>
!>  @note
!>  When called using the lambda quantities, this is one element larger in the
!>  radial.
!-------------------------------------------------------------------------------
      SUBROUTINE Spline_Fourier_Modes(ymn_vmec, ymn_spline, xm_in, istat)
      USE island_params, hs=>hs_i, ohs=>ohs_i

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(:,:), INTENT(in)  :: ymn_vmec
      REAL(dp), DIMENSION(:,:), INTENT(out) :: ymn_spline
      REAL(dp), DIMENSION(:), INTENT(in)    :: xm_in
      INTEGER, INTENT(inout)                :: istat

!  Local variables
      INTEGER                               :: js
      INTEGER                               :: modes
      INTEGER                               :: mp
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: snodes_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: y2_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: fac1
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: y_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: snodes
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: fac2
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: y_spline
      REAL(dp)                              :: hsv
      REAL(dp)                              :: expm
      REAL(dp)                              :: offset

!  Local parameters
      REAL(dp), PARAMETER                   :: one=1

!  Start of executable code

!  Check array sizes.
      IF (SIZE(ymn_vmec, 2) .le. 1) THEN
         istat = 1
         RETURN
      END IF

      IF (ns_i .LE. 1) THEN
         istat = 2
         RETURN
      END IF

      ALLOCATE(snodes_vmec(SIZE(ymn_vmec, 2)))
      ALLOCATE(y2_vmec(SIZE(ymn_vmec, 2)))
      ALLOCATE(fac1(SIZE(ymn_vmec, 2)))
      ALLOCATE(y_vmec(SIZE(ymn_vmec, 2)))

      ALLOCATE(snodes(SIZE(ymn_spline, 2)))
      ALLOCATE(fac2(SIZE(ymn_spline, 2)))
      ALLOCATE(y_spline(SIZE(ymn_spline, 2)))

      IF (SIZE(snodes_vmec) .eq. ns_vmec) THEN
         offset = 1.0_dp
      ELSE
         offset = 1.5_dp
         snodes_vmec(SIZE(snodes_vmec)) = 1
         fac1(SIZE(snodes_vmec)) = 1;
         fac2(SIZE(snodes)) = 1
         snodes(SIZE(snodes)) = 1
      END IF

!  Set up "knots" on initial (vmec, svmec ~ phi) mesh and factor for taking out
!  (putting back) [sqrt(s)]**m factor for odd m.

      hsv = one/(ns_vmec - 1)
      snodes_vmec(1) = 0;  fac1(1) = 1;
      DO js = 2, ns_vmec
         snodes_vmec(js) = hsv*(js - offset)
         fac1(js) = one/SQRT(snodes_vmec(js))     !regularization factor: sqrt(FLUX) for odd modes
      END DO
      fac1(2:ns_vmec) = one/SQRT(snodes_vmec(2:ns_vmec))

!  Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!  if s_siesta ~ sfinal, this is then the original vmec mesh.
      ohs = ns_i - 1
      hs = one/ohs
      snodes(1) = 0
      DO js = 2, ns_i
         snodes(js) = hs*(js - offset)
      END DO
      IF (iscale) THEN
         snodes = snodes*snodes !SIESTA s==fac2 ~ sqrt(s-vmec) mesh
      END IF
      fac2 = SQRT(snodes)

      DO modes = 1, SIZE(xm_in)
         expm = 0  ! MRC Check this. Orginal code looked like this would always
                   ! be non zero if the first loop iteration reset it.

         y_vmec = ymn_vmec(modes,:)
         mp = xm_in(modes)

         IF (istat.EQ.0 .and. mp.GT.0) THEN
            IF (MOD(mp,2) .eq. 1) THEN 
               expm = 1
            ELSE
               expm = 2
            END IF
            y_vmec = y_vmec*(fac1**expm)
            IF (mp .le. 2) y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
         END IF

!  Initialize spline for each mode amplitude (factor out sqrt(s) factor for
!  odd-m) (spline routines in LIBSTELL Miscel folder).
         CALL spline(snodes_vmec, y_vmec, SIZE(snodes_vmec),            &
                     lower_b, upper_b, y2_vmec)

!  Interpolate onto snodes mesh.
         DO js = 1, SIZE(snodes)
            CALL splint(snodes_vmec, y_vmec, y2_vmec,                   &
                        SIZE(snodes_vmec), snodes(js), y_spline(js))
         END DO

         IF (istat.eq.0 .and. mp.gt.0 .and. expm.ne.0) THEN
            y_spline = y_spline*(fac2**expm)
         END IF
         ymn_spline(modes,:) = y_spline(:)

!  IDENTITY TEST THAT INTERPOLATION WORKS
         IF (lverbose .and.                                             &
             .not.iscale .and.                                          &
             SIZE(snodes).EQ.SIZE(snodes_vmec) .and.                    &
             ANY(ABS(ymn_vmec(modes,2:) - ymn_spline(modes,2:)).GT.1.E-10_dp)) THEN
            PRINT *,'MODE: ',modes,'y_vmec .ne. y_spline'
         END IF

      END DO

      istat = 0

      DEALLOCATE(snodes_vmec)
      DEALLOCATE(y2_vmec)
      DEALLOCATE(fac1)
      DEALLOCATE(y_vmec)

      DEALLOCATE(snodes)
      DEALLOCATE(fac2)
      DEALLOCATE(y_spline)

      END SUBROUTINE Spline_Fourier_Modes
      
!-------------------------------------------------------------------------------
!>  @brief Spline reinterpolate a VMEC quantity.
!>
!>  Redistributes the VMEC quantity from the VMEC radial grid to he SIESTA
!>  radial grid.
!>
!>  @param[in]  y_vmec   Orginal VMEC quantity.
!>  @param[out] y_spline Splined VMEC quantity.
!>  @param[out] istat    Error status.
!-------------------------------------------------------------------------------
      SUBROUTINE Spline_OneD_Array (y_vmec, y_spline, istat)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(ns_vmec), INTENT(in)  :: y_vmec
      REAL(dp), DIMENSION(ns_i)   , INTENT(out) :: y_spline
      INTEGER, INTENT(OUT)                      :: istat

!  Local variables
      INTEGER                                   :: js
      REAL(dp), DIMENSION(ns_vmec)              :: snodes_vmec
      REAL(dp), DIMENSION(ns_vmec)              :: y2_vmec
      REAL(dp), DIMENSION(ns_i)                 :: snodes
      REAL(dp)                                  :: hs_vmec

!  Local parameters
      REAL(dp), PARAMETER                       :: one = 1

!  Start of executable code

!  Check grid dimensions.
      IF (ns_vmec .le. 1) THEN
         istat = 1
         RETURN
      END IF

      IF (ns_i .LE. 1) THEN
         istat = 2
         RETURN
      END IF

!  Set up "knots" on initial (vmec, svmec ~ phi) mesh.
      hs_vmec = one/(ns_vmec-1)
      DO js = 1, ns_vmec
         snodes_vmec(js) = hs_vmec*(js-1)
      END DO

!  Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!  if s_siesta ~ sfinal, this is then the original vmec mesh.
      ohs_i = ns_i - 1
      hs_i = one/ohs_i
      DO js = 1, ns_i
         snodes(js) = hs_i*(js-1) !vmec mesh   s-siesta~s-vmec
         IF (iScale) THEN
            snodes(js) = snodes(js)*snodes(js) !polar mesh, s-siesta~s-vmec**2
         END IF
      END DO

!  Initialize spline for each mode amplitude (factor out sqrt(s) factor for
!  odd-m)
      CALL spline(snodes_vmec, y_vmec, ns_vmec, lower_b, upper_b, y2_vmec)

!  Interpolate onto snodes mesh
      DO js = 1, ns_i
         CALL splint(snodes_vmec, y_vmec, y2_vmec, ns_vmec,              &
                     snodes(js), y_spline(js))
      END DO

      END SUBROUTINE Spline_OneD_Array

!-------------------------------------------------------------------------------
!>  @brief Add ghost points to at js=1 (s=0) and js=ns+1 (s=1).
!>
!>  Extends the VMEC lmn to the ns+1 grid.
!>
!>  @param[in]    VMEC Fourier quantity.
!>  @param[inout] Ghost Fourier quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE Add_Ghost_Points(ymn_vmec, ymn_ghost)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(:,:), INTENT(in)  :: ymn_vmec
      REAL(dp), DIMENSION(:,:), INTENT(out) :: ymn_ghost

!  Local variables
      INTEGER                      :: modes

!  Start of executable code
      ymn_ghost(:,2:ns_vmec) = ymn_vmec(:,2:ns_vmec)
      ymn_ghost(:,ns_vmec+1) = 2*ymn_vmec(:,ns_vmec) - ymn_vmec(:,ns_vmec-1)

      WHERE (xm_vmec .eq. 0)
         ymn_ghost(:,1) = 2*ymn_vmec(:,2) - ymn_vmec(:,3)
      ELSE WHERE
         ymn_ghost(:,1) = 0
      END WHERE

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Load the R, Z and lambda arrays from VMEC.
!>
!>  R, Z and Lambda are recast onto the SIESTA grid.
!>
!>  @param[out] istat Error status.
!-------------------------------------------------------------------------------
      SUBROUTINE LoadRZL_VMEC (istat)
      USE fourier, ONLY: toijsp, init_fourier, f_cos, f_sin,             &
                         f_none, f_sum, f_du, f_dv

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(out) :: istat

!  Local variables
      INTEGER              :: m
      INTEGER              :: n

!  Start executable code.

!  LOAD FOURIER FACTORS USING VMEC mpol,ntor AND ISLAND nu,nv
!  IDEA: USE VMEC mpol, ntor VALUES AND CALL fixarray, FLIP VMEC n -> -n in
!  REPACK
      
      CALL init_fourier !  MRC: FIXME: move this some where else. This will be bypassed in the new restart.

!  COMPUTE R, Z (Metric elements, jacobian), LAMBDA (For initializing B)
!
!  FIRST, REPACK SPLINED ARRAYS FROM VMEC-ORDERING TO SIESTA-ORDERING AND SWAP
!  THE N->-N MODES TO BE CONSISTENT WITH mu+nv SIESTA ARGUMENT
      CALL vmec_info_construct_island(mpol_i, ntor_i, ns_i, lasym)
      CALL repack(rmnc_i, zmns_i, lmns_i, jcurrumnc, jcurrvmnc,         &
                  rmnc_spline, zmns_spline, lmns_spline,                &
                  currumnc_spline, currvmnc_spline, isym)
      IF (lasym)                                                        &
      CALL repack(rmns_i, zmnc_i, lmnc_i, jcurrumns, jcurrvmns,         &
                  rmns_spline, zmnc_spline, lmnc_spline,                &
                  currumns_spline, currvmns_spline, isym)
      CALL vmec_info_destruct_spline

!  COMPUTE AND STORE R, Z AND THEIR ANGULAR DERIVATIVES, AND LAMBDA (NEED FOR
!  VECTOR POT)
      ALLOCATE (r1_i(nu_i, nv_i, ns_i), z1_i(nu_i, nv_i, ns_i),         &
                ru_i(nu_i, nv_i, ns_i), zu_i(nu_i, nv_i, ns_i),         &
                rv_i(nu_i, nv_i, ns_i), zv_i(nu_i, nv_i, ns_i),         &
                stat = istat)

!      IF (istat .NE. 0) WRITE (*,*) 'Allocation failed in LoadRZL_VMEC'

      CALL toijsp (rmnc_i, r1_i, f_none, f_cos)
      CALL toijsp (rmnc_i, ru_i, f_du,   f_cos)
      CALL toijsp (rmnc_i, rv_i, f_dv,   f_cos)

      CALL toijsp (zmns_i, z1_i, f_none, f_sin)
      CALL toijsp (zmns_i, zu_i, f_du,   f_sin)
      CALL toijsp (zmns_i, zv_i, f_dv,   f_sin)

      IF (lasym) THEN
         m = IOR(f_sum, f_du)
         n = IOR(f_sum, f_dv)
         CALL toijsp (rmns_i, r1_i, f_sum, f_sin)
         CALL toijsp (rmns_i, ru_i, m,     f_sin)
         CALL toijsp (rmns_i, rv_i, n,     f_sin)

         CALL toijsp (zmnc_i, z1_i, f_sum, f_cos)
         CALL toijsp (zmnc_i, zu_i, m,     f_cos)
         CALL toijsp (zmnc_i, zv_i, n,     f_cos)
      END IF

      rmax = MAXVAL(r1_i);      rmin = MINVAL(r1_i)
      zmax = MAXVAL(z1_i);      zmin = MINVAL(z1_i)
      IF (.NOT.lasym) zmin = MIN(zmin, -zmax)                           !ONLY TOP HALF

#if defined(NOSKS)
!  COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN
      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)

!  COMPUTE FULL-MESH METRIC ELEMENTS
      CALL full_mesh_metrics
#endif

      END SUBROUTINE LoadRZL_VMEC

!-------------------------------------------------------------------------------
!>  @brief Resets VMEC Fourier quanties to use SIESTA convension.
!>
!>  The splined arrays (rmnc_spl, etc) are all in VMEC ordering.
!>
!>    rmnc_spl(1:mnmax,ns_i)
!>
!>  Now pack them so they can be used by ISLAND Fourier routines.
!>
!>    rmn(0:mpol, -ntor:ntor, ns_i)
!>
!>  Arrays passed into the spline arguments get deallocated.
!>
!>  @param[inout] isl    Repacked Fourier coefficents.
!>  @param[in]    spl    Splined Fourier coefficents.
!>  @param[in]    parity True if using stellarator symmetry.
!>
!>  @note
!>  mpol_i == mpol_vmec, ntor_i == ntor_vmec here. Also, must flip the sign of
!>  'n' terms since VMEC uses mu-nv and SIESTA uses mu+nv for args of cos and
!>  sin.
!>
!>  Keep lmn on half mesh and add extrapolation points at s=1 (rho=0)(ns+1 -> 1)
!>
!>  Rountines here assume that *mn have 0:mpol_i and -ntor_i:ntor_i array
!>  bounds.
!-------------------------------------------------------------------------------
      SUBROUTINE repack (rmn, zmn, lmn, currumn, currvmn,               &
                         rmn_spl, zmn_spl, lmn_spl,                     &
                         currumn_spl, currvmn_spl, parity)
      USE fourier, ONLY: orthonorm

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rmn
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: zmn
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: lmn
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: currumn
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: currvmn
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: rmn_spl
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: zmn_spl
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: lmn_spl
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: currumn_spl
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: currvmn_spl
      INTEGER, INTENT(in)                     :: parity

!  Local variables
      INTEGER                                 :: modes
      INTEGER                                 :: js
      INTEGER                                 :: m
      INTEGER                                 :: n
      INTEGER                                 :: n1
      INTEGER                                 :: s1
      INTEGER                                 :: s2
      INTEGER                                 :: istat

!  Start executable code.
      CALL ASSERT(SIZE(rmn_spl,1).EQ.mnmax,'rmn_spl wrong size1 in REPACK')
      CALL ASSERT(SIZE(rmn_spl,2).EQ.ns_i,'rmn_spl wrong size2 in REPACK')
      CALL ASSERT(SIZE(lmn_spl,2).EQ.ns_i+1,'lmn_spl wrong size2 in REPACK')

      IF (parity .EQ. isym) THEN
         s1 = 1
      ELSE
         s2 = 1
      END IF

      IF (iam.EQ.0 .AND. lverbose .AND. parity.EQ.isym) THEN
         m = INT(MAXVAL(xm_vmec))
         n = INT(MAXVAL(ABS(xn_vmec)))/nfp_i
         IF (m.GT.mpol_i .OR. n.GT.ntor_i) WRITE(*,1000) m, n
         CALL ASSERT(MOD(nfp_vmec, nfp_i) .EQ. 0, 'nfpin should an even ' //   &
                                                  'divisor of the VMEC nfp.')
      ENDIF

!  FLIP VMEC n -> SIESTA -n SO TRIG ARG IS mu+nv IN SIESTA
      DO modes = 1,mnmax
         m = xm_vmec(modes)
         n = xn_vmec(modes)/nfp_i
         n1 = ABS(n)
         IF (m.GT.mpol_i .OR. n1.GT.ntor_i) CYCLE
!
!     LOAD n>=0 ONLY FOR M=0
!
         IF (m .EQ. 0) THEN
            IF (parity .EQ. isym)  s2 = -SIGN(1,n)
            IF (parity .EQ. iasym) s1 = -SIGN(1,n)
            rmn(m,n1,:) = rmn(m,n1,:) + s1*rmn_spl(modes,:)
            zmn(m,n1,:) = zmn(m,n1,:) + s2*zmn_spl(modes,:)
            lmn(m,n1,:) = lmn(m,n1,:) + s2*lmn_spl(modes,:)*iflipj
            currumn(m,n1,:) = currumn(m,n1,:) + s1*currumn_spl(modes,:)
            currvmn(m,n1,:) = currvmn(m,n1,:) + s1*currvmn_spl(modes,:)

         ELSE
            rmn(m,-n*iflipj,:) = rmn_spl(modes,:)
            zmn(m,-n*iflipj,:) = zmn_spl(modes,:)*iflipj
            lmn(m,-n*iflipj,:) = lmn_spl(modes,:)
            currumn(m,-n*iflipj,:) = currumn_spl(modes,:)
            currvmn(m,-n*iflipj,:) = currvmn_spl(modes,:)
         ENDIF
      END DO

!
! RS (10/03/06)  Divide out "orthonorm" factor used in Fourier FFT routines
!                mu0 already multiplied into currents (in GetRZLSpline)
! 
      DO js = 1, ns_i
         lmn(:,:,js) = lmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         rmn(:,:,js) = rmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         zmn(:,:,js) = zmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         currumn(:,:,js) = 2*hs_i*(js - 1)*currumn(:,:,js)      &
                         / orthonorm(0:mpol_i,-ntor_i:ntor_i)
         currvmn(:,:,js) = 2*hs_i*(js - 1)*currvmn(:,:,js)      &
                         / orthonorm(0:mpol_i,-ntor_i:ntor_i)
      END DO
      lmn(:,:,ns_i + 1) = lmn(:,:,ns_i + 1)                     &
                        / orthonorm(0:mpol_i,-ntor_i:ntor_i)

1000  FORMAT(' It is recommended to increase the number of modes to at' &
             ' least VMEC values mpol=',i4,' ntor=',i4)
      END SUBROUTINE repack

!-------------------------------------------------------------------------------
!>  @brief Compute the metric elements on the half mesh.
!>
!>  @param[in] r1_i R on the full mesh.
!>  @param[in] ru_i dRdu on the full mesh.
!>  @param[in] rv_i dRdv on the full mesh.
!>  @param[in] z1_i Z on the full mesh.
!>  @param[in] zu_i dZdu on the full mesh.
!>  @param[in] zv_i dZdv on the full mesh.
!-------------------------------------------------------------------------------
      SUBROUTINE half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
      USE island_params, nuv=>nuv_i, nu=>nu_i, nv=>nv_i!, ns=>ns_i

      IMPLICIT NONE  

!  Declare Arguments
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: r1_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: ru_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: rv_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: z1_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: zu_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: zv_i

!  Local variables
      INTEGER                                   :: js
      INTEGER                                   :: js1
      INTEGER                                   :: istat
      INTEGER                                   :: imin(2)
      INTEGER                                   :: imax(2)
      REAL(dp)                                  :: mintest
      REAL(dp)                                  :: maxtest
      REAL(dp), DIMENSION(nuv)                  :: r12
      REAL(dp), DIMENSION(nuv)                  :: rs12
      REAL(dp), DIMENSION(nuv)                  :: ru12
      REAL(dp), DIMENSION(nuv)                  :: rv12
      REAL(dp), DIMENSION(nuv)                  :: zs12
      REAL(dp), DIMENSION(nuv)                  :: zu12
      REAL(dp), DIMENSION(nuv)                  :: zv12

!  Local parameters.
      REAL(dp), PARAMETER                     :: p5 = 1.d0/2.d0
      REAL(dp), PARAMETER                     :: zero = 0

!  Start executable code.

!  ALLOCATE METRIC ELEMENT ARRAYS
      ALLOCATE(sqrtg(nuv,ns_i),                                         &
               gss(nuv,ns_i), gsu(nuv,ns_i), gsv(nuv,ns_i),             &
               guu(nuv,ns_i), guv(nuv,ns_i), gvv(nuv,ns_i), stat=istat)
!      IF (istat .ne. 0) CALL ErrorAbort('ALLOCATION ERROR1 IN HALF_MESH_METRICS')

!  COMPUTE ALL ON THE HALF MESH
      sqrtg(:,1) = 0

      DO js = 2, ns_i
         js1  = js - 1
         r12  = p5*(r1_i(:,js) + r1_i(:,js1))
         rs12 =    (r1_i(:,js) - r1_i(:,js1))*ohs_i
         ru12 = p5*(ru_i(:,js) + ru_i(:,js1))
         rv12 = p5*(rv_i(:,js) + rv_i(:,js1))
         zs12 =    (z1_i(:,js) - z1_i(:,js1))*ohs_i
         zu12 = p5*(zu_i(:,js) + zu_i(:,js1))
         zv12 = p5*(zv_i(:,js) + zv_i(:,js1))
         guu(:,js) = ru12*ru12 + zu12*zu12
         guv(:,js) = ru12*rv12 + zu12*zv12
         gvv(:,js) = r12*r12 + rv12*rv12 + zv12*zv12
         gsu(:,js) = rs12*ru12 + zs12*zu12
         gsv(:,js) = rs12*rv12 + zs12*zv12
         gss(:,js) = rs12*rs12 + zs12*zs12
         sqrtg(:,js) = r12*(ru12*zs12 - rs12*zu12)
      END DO

      sqrtg(:,1) = sqrtg(:,2)

!NEED THESE FOR CURRENT CALCULATION AT ORIGIN
      gss(:,1) = gss(:,2);  gsu(:,1) = 0;  gsv(:,1) = gsv(:,2)
      guu(:,1) = 0;         guv(:,1) = 0;  gvv(:,1) = gvv(:,2)

      mintest = MINVAL(sqrtg(:,2:))
      maxtest = MAXVAL(sqrtg(:,2:))
      
      IF (mintest*maxtest .LE. zero) THEN
         imin = MINLOC(sqrtg(:,2:)); imax = MAXLOC(sqrtg(:,2:))
         IF (iam .EQ. 0) WRITE(*,1000) ' MIN-G: ', mintest, imin(2),    &
                                       ' MAX-G: ', maxtest, imax(2)
         CALL ASSERT(mintest*maxtest.GT.zero,' Jacobian changed sign in half_mesh_metrics')
      END IF

1000  FORMAT(2(a,1p,e12.4,' AT JS: ',i5))
      END SUBROUTINE half_mesh_metrics

!-------------------------------------------------------------------------------
!>  @brief Gets the full grid metric elements.
!>
!>  This subroutine gets the lower metric elements on the full mesh. Preserves
!>  positive definiteness of metric tensor.
!-------------------------------------------------------------------------------
      SUBROUTINE full_mesh_metrics
      USE island_params, nuv=>nuv_i, nzeta=>nv_i, ntheta=>nu_i!, ns=>ns_i

      IMPLICIT NONE  

!  Local variables
      INTEGER                                 :: istat
      INTEGER                                 :: js
      REAL(dp), DIMENSION(nuv)                :: det
#if defined(DEBUG)
      REAL(dp)                                :: eps
      REAL(dp)                                :: maxtest
      REAL(dp)                                :: r1
      REAL(dp)                                :: s1
      REAL(dp)                                :: t1
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: temps
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: tempu
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: tempv
#endif

!  Start executable code.
      ALLOCATE(gssf(nuv,ns_i), guuf(nuv,ns_i), gvvf(nuv,ns_i),           &
               gsuf(nuv,ns_i), gsvf(nuv,ns_i), guvf(nuv,ns_i),           &
               hss(nuv,ns_i),  huu(nuv,ns_i),  hvv(nuv,ns_i),            &
               hsu(nuv,ns_i),  hsv(nuv,ns_i),  huv(nuv,ns_i), stat=istat)

      CALL to_full_mesh(gss,gssf)
      CALL to_full_mesh(guu,guuf)
      CALL to_full_mesh(gvv,gvvf)
      CALL to_full_mesh(gsu,gsuf)
      CALL to_full_mesh(gsv,gsvf)
      CALL to_full_mesh(guv,guvf)
      guuf(:,1) = 0;  gsuf(:,1) = 0;  guvf(:,1) = 0

!  Compute upper metric elements (inverse of lower matrix) on full mesh.
!  NOTE: DET == 1/det|Gij| == 1/sqrtg(full)**2
      hss(:,1) = 0;  hsu(:,1) = 0;  hsv(:,1) = 0
      huu(:,1) = 0;  huv(:,1) = 0;  hvv(:,1) = 0

      DO js = 2, ns_i
         det(:) =                                                        &
                gssf(:,js)*(guuf(:,js)*gvvf(:,js)-guvf(:,js)*guvf(:,js)) &
              + gsuf(:,js)*(guvf(:,js)*gsvf(:,js)-gsuf(:,js)*gvvf(:,js)) &
              + gsvf(:,js)*(gsuf(:,js)*guvf(:,js)-gsvf(:,js)*guuf(:,js))
         IF (ANY(det .LE. zero)) THEN
            IF (iam .EQ. 0) PRINT *,'Determinant |gijf| <= 0 at js = ',js
            IF (js .EQ. ns_i) THEN
               det(:) = sqrtg(:,ns_i)**2
               gssf(:,ns_i) = gss(:,ns_i)
               guuf(:,ns_i) = guu(:,ns_i)
               gvvf(:,ns_i) = gvv(:,ns_i)
               gsuf(:,ns_i) = gsu(:,ns_i)
               guvf(:,ns_i) = guv(:,ns_i)
               gsvf(:,ns_i) = gsv(:,ns_i)
            ELSE
               det(:) = ((sqrtg(:,js) + sqrtg(:,js + 1))/2)**2
!               CALL ErrorAbort('Full_mesh_metrics gijf')
            END IF
         END IF
         det = one/det
         hss(:,js) = det(:)*                                             &
                    (guuf(:,js)*gvvf(:,js) - guvf(:,js)*guvf(:,js))
         hsu(:,js) = det(:)*                                             &
                    (guvf(:,js)*gsvf(:,js) - gsuf(:,js)*gvvf(:,js))
         hsv(:,js) = det(:)*                                             &
                    (gsuf(:,js)*guvf(:,js) - gsvf(:,js)*guuf(:,js))
         huu(:,js) = det(:)*                                             &
                    (gssf(:,js)*gvvf(:,js) - gsvf(:,js)*gsvf(:,js))
         huv(:,js) = det(:)*                                             &
                    (gsvf(:,js)*gsuf(:,js) - gssf(:,js)*guvf(:,js))
         hvv(:,js) = det(:)*                                             &
                    (gssf(:,js)*guuf(:,js) - gsuf(:,js)*gsuf(:,js))

      END DO

#if defined(DEBUG)
!  CONFIRM ACCURACY OF INVERSE
      eps = 0.01_dp*SQRT(EPSILON(eps))

      ALLOCATE(temps(ntheta,nzeta,ns_i),                                 &
               tempu(ntheta,nzeta,ns_i),                                 &
               tempv(ntheta,nzeta,ns_i), stat=istat)

!  gijf * hj1   (first column of hij matrix)
      CALL tolowerf(hss, hsu, hsv, temps, tempu, tempv, 1, ns_i)
      r1 = MAXVAL(ABS(temps(:,:,2:) - 1))
      s1 = MAXVAL(ABS(tempu(:,:,2:)))
      t1 = MAXVAL(ABS(tempv(:,:,2:)))
      maxtest = MAX(r1,s1,t1)
!      IF (maxtest.GT.eps .AND. iam.EQ.0) CALL ErrorAbort('ERROR1 IN FULL_MESH_METRICS')

!  gijf * hj2   (second column of hij matrix)
      CALL tolowerf(hsu, huu, huv, temps, tempu, tempv, 1, ns_i)
      r1 = MAXVAL(ABS(temps(:,:,2:)))
      s1 = MAXVAL(ABS(tempu(:,:,2:) - 1))
      t1 = MAXVAL(ABS(tempv(:,:,2:)))
      maxtest = MAX(r1,s1,t1)
!      IF (maxtest.GT.eps .AND. iam.EQ.0) CALL ErrorAbort('ERROR2 IN FULL_MESH_METRICS')

!  gijf * hj3   (third column of hij matrix)
      CALL tolowerf(hsv, huv, hvv, temps, tempu, tempv, 1, ns_i)
      r1 = MAXVAL(ABS(temps(:,:,2:)))
      s1 = MAXVAL(ABS(tempu(:,:,2:)))
      t1 = MAXVAL(ABS(tempv(:,:,2:) - 1))
      maxtest = MAX(r1,s1,t1)
!      IF (maxtest.GT.eps .AND. iam.EQ.0) CALL ErrorAbort('ERROR3 IN FULL_MESH_METRICS')

      DEALLOCATE(temps, tempu, tempv)
#endif
      END SUBROUTINE full_mesh_metrics

!-------------------------------------------------------------------------------
!>  @brief Converts to contravariant component full grid.
!>
!>  Computes contravariant (upper) from covariant (lower)  components on the
!>  full mesh.
!>
!>    xsupsij = hss*xsubsij + hsu*xsupuij + hsv*xsubvij
!>    xsupuij = hsu*xsubsij + huu*xsupuij + huv*xsubvij
!>    xsupvij = hsv*xsubsij + huv*xsupuij + hvv*xsubvij
!>
!>  @param[in]  xsubsij Covariant s component.
!>  @param[in]  xsubuij Covariant u component.
!>  @param[in]  xsubvij Covariant v component.
!>  @param[out] xsupsij Contravariant s component.
!>  @param[out] xsupuij Contravariant u component.
!>  @param[out] xsupvij Contravariant v component.
!>  @param[in]  nsmin   Min radial index.
!>  @param[in]  nsmax   Max radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE toupper(xsubsij, xsubuij, xsubvij,                     &
                         xsupsij, xsupuij, xsupvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, nuv=>nuv_i

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsubsij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsubuij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsubvij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsupsij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsupuij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsupvij
      INTEGER, INTENT(in)                               :: nsmin
      INTEGER, INTENT(in)                               :: nsmax

!  Local variables
        INTEGER                                         :: js

!  Start executable code.
      DO js = nsmin, nsmax
         xsupsij(:,js) = hss(:,js)*xsubsij(:,js)                        &
                       + hsu(:,js)*xsubuij(:,js)                        &
                       + hsv(:,js)*xsubvij(:,js)
         xsupuij(:,js) = hsu(:,js)*xsubsij(:,js)                        &
                       + huu(:,js)*xsubuij(:,js)                        &
                       + huv(:,js)*xsubvij(:,js)
         xsupvij(:,js) = hsv(:,js)*xsubsij(:,js)                        &
                       + huv(:,js)*xsubuij(:,js)                        &
                       + hvv(:,js)*xsubvij(:,js)
      END DO

      END SUBROUTINE toupper

!-------------------------------------------------------------------------------
!>  @brief Converts to covariant component half grid.
!>
!>  Computes covariant (lower) from contravariant (upper) magnetic field
!>  components on the half mesh.
!>
!>    xsubsij = gss*xsupsij + gsu*xsupuij + gsv*xsupvij
!>    xsubuij = gsu*xsupsij + guu*xsupuij + guv*xsupvij
!>    xsubvij = gsv*xsupsij + guv*xsupuij + gvv*xsupvij
!>
!>  @param[in]  xsupsij Contravariant s component.
!>  @param[in]  xsupuij Contravariant u component.
!>  @param[in]  xsupvij Contravariant v component.
!>  @param[out] xsubsij Covariant s component.
!>  @param[out] xsubuij Covariant u component.
!>  @param[out] xsubvij Covariant v component.
!>  @param[in]  nsmin   Min radial index.
!>  @param[in]  nsmax   Max radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE tolowerh(xsupsij, xsupuij, xsupvij,                      &
                          xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: nuv=>nuv_i

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsupsij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsupuij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsupvij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsubsij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsubuij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsubvij
      INTEGER, INTENT(IN)                               :: nsmin
      INTEGER, INTENT(IN)                               :: nsmax

!  Start executable code.
      xsubsij = gss(:,nsmin:nsmax)*xsupsij                               &
              + gsu(:,nsmin:nsmax)*xsupuij                               &
              + gsv(:,nsmin:nsmax)*xsupvij

      xsubuij = gsu(:,nsmin:nsmax)*xsupsij                               &
              + guu(:,nsmin:nsmax)*xsupuij                               &
              + guv(:,nsmin:nsmax)*xsupvij
      
      xsubvij = gsv(:,nsmin:nsmax)*xsupsij                               &
              + guv(:,nsmin:nsmax)*xsupuij                               &
              + gvv(:,nsmin:nsmax)*xsupvij

      END SUBROUTINE tolowerh

!-------------------------------------------------------------------------------
!>  @brief Converts to covariant component full grid.
!>
!>  Computes covariant (lower) from contravariant (upper) magnetic field
!>  components on the full mesh.
!>
!>    xsubsij = gssf*xsupsij + gsuf*xsupuij + gsvf*xsupvij
!>    xsubuij = gsuf*xsupsij + guuf*xsupuij + guvf*xsupvij
!>    xsubvij = gsvf*xsupsij + guvf*xsupuij + gvvf*xsupvij
!>
!>  @param[in]  xsupsij Contravariant s component.
!>  @param[in]  xsupuij Contravariant u component.
!>  @param[in]  xsupvij Contravariant v component.
!>  @param[out] xsubsij Covariant s component.
!>  @param[out] xsubuij Covariant u component.
!>  @param[out] xsubvij Covariant v component.
!>  @param[in]  nsmin   Min radial index.
!>  @param[in]  nsmax   Max radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE tolowerf(xsupsij, xsupuij, xsupvij,                     &
                          xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: ns=>ns_i, nuv=>nuv_i

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nuv,ns), INTENT(in)  :: xsupsij  !  MRC: Why ns and not nsmin:nsmax here?
      REAL(dp), DIMENSION(nuv,ns), INTENT(in)  :: xsupuij
      REAL(dp), DIMENSION(nuv,ns), INTENT(in)  :: xsupvij
      REAL(dp), DIMENSION(nuv,ns), INTENT(out) :: xsubsij
      REAL(dp), DIMENSION(nuv,ns), INTENT(out) :: xsubuij
      REAL(dp), DIMENSION(nuv,ns), INTENT(out) :: xsubvij
      INTEGER, INTENT(in)                        :: nsmin
      INTEGER, INTENT(in)                        :: nsmax

!  Local variables
      INTEGER                                    :: js
      INTEGER                                    :: js1

!  Start executable code.
      DO js = nsmin, nsmax
         js1 = js - nsmin + 1
         xsubsij(:,js1) = gssf(:,js)*xsupsij(:,js1)                      &
                        + gsuf(:,js)*xsupuij(:,js1)                      &
                        + gsvf(:,js)*xsupvij(:,js1)

         xsubuij(:,js1) = gsuf(:,js)*xsupsij(:,js1)                      &
                        + guuf(:,js)*xsupuij(:,js1)                      &
                        + guvf(:,js)*xsupvij(:,js1)

         xsubvij(:,js1) = gsvf(:,js)*xsupsij(:,js1)                      &
                        + guvf(:,js)*xsupuij(:,js1)                      &
                        + gvvf(:,js)*xsupvij(:,js1)
      END DO

      END SUBROUTINE tolowerf

!-------------------------------------------------------------------------------
!>  @brief Deallocate memory containing metric elements on the half mesh.
!>
!>  Also deallocates the grid variables.
!>
!>  @note
!>  sqrtg is deallocated in @ref init_bcovar and stored in @ref jacob variable
!-------------------------------------------------------------------------------
      SUBROUTINE cleanup_metric_elements
      USE fourier, ONLY: dealloc_fourier
      IMPLICIT NONE

!  Local variables
      INTEGER :: istat

!  Start executable code.
      DEALLOCATE(gss, gsu, gsv, guu, guv, gvv,                           &
                 hss, hsu, hsv, huu, huv, hvv,                           &
                 phipf_i, chipf_i, presf_i, stat=istat)
      
      DEALLOCATE (r1_i, z1_i, ru_i, zu_i, rv_i, zv_i)

      CALL dealloc_fourier
      CALL vmec_info_destruct_island

      END SUBROUTINE cleanup_metric_elements

!-------------------------------------------------------------------------------
!>  @brief Deallocate memory containing metric elements on the full mesh.
!-------------------------------------------------------------------------------
      SUBROUTINE dealloc_full_lower_metrics
      USE stel_kinds

      IMPLICIT NONE

!  Local variables
      INTEGER :: istat

!  Start executable code.
      DEALLOCATE(gssf, guuf, gvvf, gsuf, gsvf, guvf, stat = istat)
!      IF (istat .NE. 0) CALL ErrorAbort('Problem dealloc. in LOWER_METRIC_FULL')

      END SUBROUTINE dealloc_full_lower_metrics

      END MODULE metrics
