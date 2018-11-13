!>
!!  \brief Module for reading in VMEC input and computing metric elements
!!   on the SIESTA (sqrt-flux) radial grid.
!<
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
!     
!     WRITTEN 06-27-06 BY S. P. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: COMPUTES AND STORES THE REAL-SPACE METRIC ELEMENTS, JACOBIAN BASED 
!     ON A SQRT(FLUX) SPLINED VMEC - COORDINATE SYSTEM
!

!     VARIABLE DECLARATIONS (3D for now; will convert to 2D or 1D as needed)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), DIMENSION(:,:), ALLOCATABLE ::                           &
                   sqrtg,                                                & !sqrt(g): Jacobian on half grid
                   gss, gsu, gsv, guu, guv, gvv,                         & !elements of lower metric tensor (half mesh)
                   hss, hsu, hsv, huu, huv, hvv                            !elements of upper metric tensor (full mesh)

      REAL(dp) :: rmax, rmin, zmax, zmin

      REAL(dp), DIMENSION(:,:), ALLOCATABLE:: gssf, guuf, gvvf,          & !Lower Metric elements (full)
                   gsuf, gsvf, guvf   

      INTEGER :: iflipj=1                                                  !=-1 to flip sign of vmec jac
!
!     LOCAL (PRIVATE) HELPER ROUTINES (NOT ACCESSIBLE FROM OUTSIDE THIS MODULE)
!
      PRIVATE spline_fourier_modes, loadrzl_vmec, add_ghost_points
      INTEGER, PARAMETER, PRIVATE :: iScale=0                              !=0, s ~ SQRT(PHI) - default; = 1, s ~ PHI (VMEC)
      INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1

      REAL(dp) :: skston, skstoff

      CONTAINS

      SUBROUTINE init_metric_elements(ns_in, mpol_in, ntor_in, wout_file)
      USE timer_mod, ONLY: init_timers
      USE shared_data, ONLY: torflux, polflux
      USE island_params, hs=>hs_i, ns=>ns_i
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*)            :: wout_file
      INTEGER, INTENT(IN)      :: ns_in, mpol_in, ntor_in
      INTEGER                  :: istat, ntype, imesh, js
      REAL(dp)                 :: t1, t2
      INTEGER                  :: m, n, mn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp),ALLOCATABLE     :: lmns_half(:,:), lmnc_half(:,:)
!-----------------------------------------------
!
!     LOADS VALUES FOR MESHES TO BE USED IN VMECPP ISLAND SOLVER
!     (GENERALLY DIFFERENT FROM VMEC MESHES. USE SUBSCRIPT _i FOR ISLAND VARIABLES)
!     READS wout_file FROM VMEC TO GET FOURIER COMPONENTS ON VMEC MESH
!     SPLINES THE VMEC COMPONENTS TO RADIAL ISLAND MESH (PHI->SQRT(PHI))
!     CALLS Fourier MODULE fixarray TO LOAD TRIG ARRAYS FOR COMPUTING ISLAND METRICS
!     COMPUTES gij (sub/sup) AND JACOBIAN ON ISLAND MESHES IN REAL SPACE
!

#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL init_timers
#if defined(SKS) 
      CALL second0(skstoff)
      init_timers_time=init_timers_time+(skstoff-skston)
#endif

      IF (iScale.NE.0 .AND. IAM.EQ.0)                                   &
        PRINT *,'RESET ISCALE=0 AFTER TESTING!'

!
!     READ DATA FROM VMEC WOUT FILE (LIBSTELL ROUTINE)
!
#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL read_wout_file(wout_file, istat)
#if defined(SKS) 
      CALL second0(skstoff)
      read_wout_file_time=read_wout_file_time+(skstoff-skston)
#endif
      CALL ASSERT(istat.EQ.0, 'Read-wout error in INIT_METRIC_ELEMENTS')

!  Asymmetric terms may be overriden in the namelist input file.
!      lasym = lasym .OR. lasym_vmec

      nfp_i = nfp_vmec
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


!Set number of points == number of modes for now! (May want mid-points for flux conservation)
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
      mnmax_i = (mpol_i + 1)*(2*ntor_i + 1)             ! Added RS. Contains total number of modes.

!
!     READ VACUUM VESSEL DATA (before calling init_fourier)
!     (points on the vessel surface are already re-distributed by VMESH code).  
   
      IF  (l_vessel) THEN
        CALL read_vessel_file(wout_file, istat)
        IF (istat .EQ. 0) THEN
          IF (iam .EQ. 0) PRINT *,'Vacuum Vessel data successfully read.'
        END IF
      END IF

!     UPDATED WOUT TO WRITE OUT CHIPF (SPH) 
!
!      iotaf_vmec = iotaf_vmec*phipf_vmec                 MRC 4/1/2016
!
!     SPLINE R, Z. L FOURIER COMPONENTS in s FROM ORIGINAL VMEC MESH (s ~ phi, ns_vmec points) 
!     TO A "POLAR" MESH [s ~ sqrt(phi), ns_i POINTS] WITH BETTER AXIS RESOLUTION
!
      CALL GetRZLSplines(rmnc_vmec, zmns_vmec, lmns_vmec,               &
                         currumnc, currvmnc,                            &
                         rmnc_spline, zmns_spline, lmns_spline,         & 
                         currumnc_spline, currvmnc_spline, isym)
      IF (lasym)                                                        &
      CALL GetRZLSplines(rmns_vmec, zmnc_vmec, lmnc_vmec,               &
                         currumns, currvmns,                            &
                         rmns_spline, zmnc_spline, lmnc_spline,         &
                         currumns_spline, currvmns_spline, iasym)

!
!     Spline 1-D arrays: careful -> convert phipf VMEC and multiply
!     by ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island
!
#if defined(SKS) 
      CALL second0(skston)
#endif
      ALLOCATE(phipf_i(ns_i), chipf_i(ns_i), presf_i(ns_i), stat=istat)

!      CALL Spline_OneD_Array (iotaf_vmec, chipf_i, istat)                       MRC 4/1/2016
      CALL Spline_OneD_Array (chipf_vmec, chipf_i, istat)                      ! MRC 4/1/2016
      CALL Spline_OneD_Array (phipf_vmec, phipf_i, istat)
      presf_vmec = mu0 * presf_vmec
      CALL Spline_OneD_Array (presf_vmec, presf_i, istat)
#if defined(SKS) 
      CALL second0(skstoff)
      Spline_OneD_Array_time=Spline_OneD_Array_time+(skstoff-skston)
#endif
!
!     Scale phipf_i and convert to sqrt(flux) mesh by multiplying by 2*s
!
      phipf_i = phipf_i / (2*pi)
      chipf_i = chipf_i / (2*pi)

!Mapping s_vmec = (s_siesta)^2, so d(s_vmec)/d(s_siesta) = 2 s_siesta, where
!s_siesta = hs_i*(js-1)

      IF (iScale == 0) THEN
         DO js = 1, ns_i
            phipf_i(js) = 2 * hs_i*(js-1) * phipf_i(js)
            chipf_i(js) = 2 * hs_i*(js-1) * chipf_i(js)
         END DO
      END IF
121   FORMAT(a,i4, 4(a,1p,e12.3))

      t1 = twopi*hs_i*(SUM(phipf_i(2:ns_i)) + SUM(phipf_i(1:ns_i-1)))/2
      IF (iam .EQ. 0) THEN
         WRITE (33, 100) volume_i, 1.E-6_dp*(wb_i+wp_i/(gamma-1))/mu0,   &
         1.E-6_dp*wp_i/mu0,                                              &
         t1,wp_i/wb_i,chipf_i(2)/phipf_i(2),                             &
         chipf_i(ns_i)/phipf_i(ns_i), gamma
      ENDIF

 100  FORMAT(/,' INITIAL PARAMETERS (FROM VMEC)',/,                      &
               ' PHYSICS QUANTITIES',/,                                  &
               ' PLASMA VOLUME (M^3): ',1pE12.4,                         &
               ' TOTAL MHD ENERGY (MJ):   ', 1pE16.8,/,                  &
               ' THERMAL ENERGY (MJ): ', 1pE12.4,                        &
               ' EDGE TOROIDAL FLUX (Wb): ', 1pE12.4,/,                  &
               ' DIMENSIONLESS QUANTITIES',/,                            &
               ' <BETA>: ',1pE12.4, ' IOTA(0): ',1pE12.4,                &
               ' IOTA(1) ',1pE12.4, ' GAMMA: ',1pE12.4,/, 21('-'),/)

! Added SKS
             ALLOCATE(torflux(ns), polflux(ns))
             torflux(1) = 0;  polflux(1) = 0;
             DO js = 2, ns
               polflux(js) = polflux(js-1) + hs*(chipf_i(js)+chipf_i(js-1))/2
               torflux(js) = torflux(js-1) + hs*(phipf_i(js)+phipf_i(js-1))/2
             END DO
! End Added SKS
!
!     CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON SQRT(FLUX) - "POLAR" - MESH
!     AND COMPUTE METRIC ELEMENTS AND JACOBIAN
!
#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL LoadRZL_VMEC(istat)
#if defined(SKS) 
      CALL second0(skstoff)
      LoadRZL_VMEC_time=LoadRZL_VMEC_time+(skstoff-skston)
#endif
      CALL ASSERT(istat.EQ.0,'LoadRZL error in INIT_METRIC_ELEMENTS')

      IF (l_vessel) THEN
!        IF (iam .EQ. 0) print *, "ns BEFORE grid_extender call: ", ns, ", ", __FILE__, __LINE__
        CALL grid_extender(wout_file, 'quad', istat)
#if defined(NOSKS)
        ns=nsp
        IF (iam .EQ. 0) THEN
        print *," "
        print *," !!! WARNING!!! Setting ns=nsp for debugging purposes !!!"
        print *," !!! Remember to unset this assignment for production runs !!!"
        print *," "
        END IF
#endif
!        IF (iam .EQ. 0) print *, "ns  AFTER grid_extender call: ", ns, ", ", __FILE__, __LINE__
      END IF
!
!     COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN
!
      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)

!     COMPUTE FULL-MESH METRIC ELEMENTS
      CALL full_mesh_metrics

!      IF (lwout_opened) CALL read_wout_deallocate  !  Don't deallocate the wout file here.

      END SUBROUTINE init_metric_elements

      SUBROUTINE SurfAvg(average, q3d, nsmin, nsmax)
        USE island_params, ONLY: cosmui 
        IMPLICIT NONE

        INTEGER, INTENT(IN)   :: nsmin, nsmax
        REAL(dp), INTENT(OUT) :: average(nsmin:nsmax)
        REAL(dp), INTENT(IN)  :: q3d(nu_i,nv_i,nsmin:nsmax)

        INTEGER :: l, js
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: wint

        ALLOCATE(wint(nu_i,nv_i,nsmin:nsmax))
        DO l = 1, nu_i
          wint(l,:,:) = cosmui(l,0)
        END DO         
        DO js = nsmin, nsmax
          average(js) = SUM(q3d(:,:,js)*wint(:,:,js))
        END DO

        DEALLOCATE(wint)
      END SUBROUTINE SurfAvg
 

      SUBROUTINE GetRZLSplines(rmn_vmec,   zmn_vmec,   lmn_vmec,        &
                               currumn_vmec, currvmn_vmec,              &
                               rmn_spline, zmn_spline, lmn_spline,      &
                               currumn_spline, currvmn_spline, ipar)
      USE island_params, ns=>ns_i
!      USE read_wout_mod, mnmax_vmec=>mnmax, ns_vmec=>ns  !SCREWS UP iasym!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar
      REAL(dp), DIMENSION(mnmax, ns_vmec), INTENT(IN) ::                &
                rmn_vmec, zmn_vmec, lmn_vmec
      REAL(dp), DIMENSION(mnmax_nyq, ns_vmec), INTENT(IN)  ::           &
                currumn_vmec, currvmn_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rmn_spline, zmn_spline,  &
                              lmn_spline, currumn_spline, currvmn_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: num_ntypes = 5
      INTEGER :: istat, ntype, modes, m2, m, n, mode_check
      REAL(dp), ALLOCATABLE     :: lmn_half(:,:), curr(:,:)
!-----------------------------------------------
      
      ALLOCATE(rmn_spline(mnmax,ns), zmn_spline(mnmax,ns),              &
               lmn_spline(mnmax,ns+1), currumn_spline(mnmax,ns),        &
               currvmn_spline(mnmax,ns), stat=istat)
      CALL ASSERT(istat.EQ.0, 'Allocation error in GetRZLSplines')

      IF (ipar.EQ.iasym .AND. .NOT.lasym_vmec) THEN
           rmn_spline = 0; zmn_spline = 0; lmn_spline = 0
           currumn_spline = 0; currvmn_spline = 0
         RETURN
      END IF

      ALLOCATE (curr(mnmax,ns_vmec), stat=istat)

       DO ntype = 1, num_ntypes
         IF (ntype .EQ. 1) THEN
           istat = 0
           CALL second0(skston)
           CALL Spline_Fourier_Modes(rmn_vmec, rmn_spline, xm_vmec,    &
             ns_vmec, ns, mnmax, istat)
           CALL second0(skstoff)
           CALL ASSERT(istat.EQ.0,'Error splining rmn')
         ELSE IF (ntype .EQ. 2) THEN
           istat = 0
           CALL second0(skston)
           CALL Spline_Fourier_Modes(zmn_vmec, zmn_spline, xm_vmec,    &
             ns_vmec, ns, mnmax, istat)   
           CALL second0(skstoff)
           CALL ASSERT(istat.EQ.0,'Error splining zmn')
         ELSE IF (ntype .EQ. 3) THEN
           ALLOCATE(lmn_half(mnmax,ns_vmec+1), stat=istat)
           CALL Add_Ghost_Points(lmn_vmec, lmn_half)
           istat = 0
           CALL second0(skston)
           CALL Spline_Fourier_Modes(lmn_half, lmn_spline, xm_vmec,     &
             ns_vmec+1, ns+1, mnmax, istat)   
           CALL second0(skstoff)
           CALL ASSERT(istat.EQ.0,'Error splining lmn')
           DEALLOCATE(lmn_half)
         ELSE IF (ntype.EQ.4 .OR. ntype.EQ.5) THEN
            curr = 0;  mode_check = 0
            DO modes = 1,mnmax
               m = xm_vmec(modes)
               n = xn_vmec(modes)
               DO m2 = 1,mnmax_nyq
                  IF (INT(xm_nyq(m2)).EQ.m .AND.                        &
                      INT(xn_nyq(m2)).EQ.n) EXIT
               END DO
               IF (m2 .GT. mnmax_nyq) CYCLE
               mode_check = mode_check+1
               IF (ntype.EQ.4) curr(modes,:) = mu0*currumn_vmec(m2,:)
               IF (ntype.EQ.5) THEN
                 curr(modes,:) = mu0*currvmn_vmec(m2,:)
                 IF (ipar.EQ.isym .AND. m.EQ.0 .AND. n.EQ.0) THEN
                   jsupvdotA = mu0*(SUM(curr(modes,2:ns_vmec-1))       &
                     + (curr(modes,1) + curr(modes,ns_vmec))/2)
                   jsupvdotA = jsupvdotA/(ns_vmec-1)
                 END IF
               END IF
           END DO

           CALL ASSERT(mode_check .EQ. mnmax,' MODE_CHECK != mnmax')

           IF (ntype.EQ.4)                                              &
             CALL Spline_Fourier_Modes(curr, currumn_spline, xm_vmec,   &
             ns_vmec, ns, mnmax, istat)   
           IF (ntype.EQ.5)                                              &
             CALL Spline_Fourier_Modes(curr, currvmn_spline, xm_vmec,   &
             ns_vmec, ns, mnmax, istat)   
         END IF
       END DO

      DEALLOCATE (curr)
     
      END SUBROUTINE GetRZLSplines
 

      ! Added SKS
      SUBROUTINE Spline_Fourier_Modes(ymn_vmec, ymn_spline, xm_in,      &
                                      nsv, nsi, mnmax_in, istat)
      USE island_params, hs=>hs_i, ohs=>ohs_i, ns=>ns_i
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)        :: nsv, nsi, mnmax_in
      INTEGER, INTENT(INOUT)     :: istat
      REAL(dp), INTENT(IN)       :: xm_in(mnmax_in)
      REAL(dp), DIMENSION(mnmax_in, nsv), TARGET,                       &
                   INTENT(IN)    :: ymn_vmec
      REAL(dp), DIMENSION(mnmax_in, nsi), TARGET,                       &
                   INTENT(OUT)   :: ymn_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER          :: one=1
      INTEGER                      :: js, modes, ntype, mp
      REAL(dp), DIMENSION(nsv)     :: snodes_vmec, y2_vmec, y_save, fac1, y_vmec
      REAL(dp), DIMENSION(nsi)     :: snodes, fac2, y_spline
      REAL(dp)                     :: hsv, yp1, ypn, expm=0
      LOGICAL                      :: lfull, lscale
!-----------------------------------------------
!
!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
!

!
!     1. Set up "knots" on initial (vmec, svmec ~ phi) mesh 
!        and factor for taking out (putting back) [sqrt(s)]**m factor for odd m
!

      IF (nsv .le. 1) THEN
         istat = 1
         RETURN
      END IF

      lfull = (nsv .eq. ns_vmec)
      lscale = (iscale .eq. 0)

      hsv = one/(ns_vmec-1)
      snodes_vmec(1) = 0;  fac1(1) = 1;
      DO js = 2, ns_vmec
         IF (lfull) THEN
            snodes_vmec(js) = hsv*(js-1)
         ELSE
            snodes_vmec(js) = hsv*(js-1.5_dp)
         END IF
         fac1(js) = one/SQRT(snodes_vmec(js))     !regularization factor: sqrt(FLUX) for odd modes
      END DO

      IF (.not.lfull) THEN
         snodes_vmec(nsv) = 1
         fac1(nsv) = 1;
      END IF

!
!     2. Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!        if s_siesta ~ sfinal, this is then the original vmec mesh
!
      IF (ns .LE. 1) THEN
         istat = 2
         RETURN
      END IF

      ohs = ns-1
      hs = one/ohs
      fac2(1) = 0;  snodes(1) = 0
      DO js = 2, ns
         IF (lfull) THEN
            fac2(js) = hs*(js-1)
         ELSE
            fac2(js) = hs*(js-1.5_dp)
         END IF
         IF (lscale) THEN
            snodes(js) = fac2(js)*fac2(js)                                !SIESTA s==fac2 ~ sqrt(s-vmec) mesh
         ELSE
            snodes(js) = fac2(js)                                         !SIESTA s==fac2 ~ s-vmec mesh
         END IF
         fac2(js) = SQRT(snodes(js))
      END DO

      IF (.not.lfull) THEN
         fac2(nsi) = 1
         snodes(nsi) = 1
      END IF

      DO modes = 1, mnmax_in
         y_vmec = ymn_vmec(modes,:)
         mp = xm_in(modes)

         IF (istat.EQ.0 .and. mp.GT.0) THEN
            IF (MOD(mp,2) .eq. 1) THEN 
               expm = 1
            ELSE
               expm = 2
            END IF
            IF (expm .ne. 0) y_vmec = y_vmec*(fac1**expm)
            IF (mp .le. 2) y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
         END IF

!
!      3. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!         (spline routines in LIBSTELL Miscel folder)
!
         yp1 = -1.e30_dp;  ypn = -1.e30_dp
         CALL spline (snodes_vmec, y_vmec, nsv, yp1, ypn, y2_vmec)

!
!      4. Interpolate onto snodes mesh
!
         DO js = 1, nsi
            CALL splint (snodes_vmec, y_vmec, y2_vmec, nsv,             &
                         snodes(js), y_spline(js))
         END DO

         IF (istat.eq.0 .and. mp.gt.0) THEN
            IF (expm .ne. 0) y_spline = y_spline*(fac2**expm)
         END IF

         ymn_spline(modes,:) = y_spline(:)


!IDENTITY TEST THAT INTERPOLATION WORKS
         IF (nsi.EQ.nsv .AND. .NOT.lscale) THEN
            IF (ANY(ABS(ymn_vmec(modes,2:) - ymn_spline(modes,2:)).GT.1.E-10_dp)) THEN
               IF (lverbose) PRINT *,'MODE: ',modes,'y_vmec .ne. y_spline'
            ENDIF
         END IF

      END DO

      istat = 0

      END SUBROUTINE Spline_Fourier_Modes

      SUBROUTINE O_Spline_Fourier_Modes(ymn_vmec, ymn_spline, nsv, nsi,    &
                                      istat)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)        :: nsv, nsi
      INTEGER, INTENT(INOUT)     :: istat
      REAL(dp), DIMENSION(mnmax, nsv), TARGET,                           &
                   INTENT(IN)  :: ymn_vmec
      REAL(dp), DIMENSION(mnmax, nsi), TARGET,                           &
                   INTENT(OUT) :: ymn_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER          :: one=1
      INTEGER                         :: js, modes, ntype, mp
      REAL(dp), DIMENSION(nsv)     :: snodes_vmec, y2_vmec, y_save, fac1, y_vmec
      REAL(dp), DIMENSION(nsi)     :: snodes, fac2, y_spline
      REAL(dp)                     :: hsv, yp1, ypn, expm=0
      LOGICAL                         :: lfull, lscale
!-----------------------------------------------
!
!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
!

!
!     1. Set up "knots" on initial (vmec, svmec ~ phi) mesh 
!        and factor for taking out (putting back) [sqrt(s)]**m factor for odd m
!

      IF (nsv .le. 1) THEN
         istat = 1
         RETURN
      END IF

      lfull = (nsv .eq. ns_vmec)
      lscale = (iscale .eq. 0)

      hsv = one/(ns_vmec-1)
      snodes_vmec(1) = 0;  fac1(1) = 1;
      DO js = 2, ns_vmec
         IF (lfull) THEN
            snodes_vmec(js) = hsv*(js-1)
         ELSE
            snodes_vmec(js) = hsv*(js-1.5_dp)
         END IF
         fac1(js) = one/SQRT(snodes_vmec(js))     !regularization factor: sqrt(FLUX) for odd modes
      END DO

      IF (.not.lfull) THEN
         snodes_vmec(nsv) = 1
         fac1(nsv) = 1;
      END IF

!
!     2. Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!        if s_siesta ~ sfinal, this is then the original vmec mesh
!
      IF (ns_i .LE. 1) THEN
         istat = 2
         RETURN
      END IF

      ohs_i = ns_i-1
      hs_i = one/ohs_i
      fac2(1) = 0;  snodes(1) = 0
      DO js = 2, ns_i
         IF (lfull) THEN
            fac2(js) = hs_i*(js-1)
         ELSE
            fac2(js) = hs_i*(js-1.5_dp)
         END IF
         IF (lscale) THEN
            snodes(js) = fac2(js)*fac2(js)                                 !SIESTA s==fac2 ~ sqrt(s-vmec) mesh
         ELSE
            snodes(js) = fac2(js)                                          !SIESTA s==fac2 ~ s-vmec mesh
         END IF
         fac2(js) = SQRT(snodes(js))
      END DO

      IF (.not.lfull) THEN
         fac2(nsi) = 1
         snodes(nsi) = 1
      END IF

      DO modes = 1, mnmax
         y_vmec = ymn_vmec(modes,:)
         mp = xm_vmec(modes)

         IF (istat.EQ.0 .and. mp.GT.0) THEN
            IF (MOD(mp,2) .eq. 1) THEN 
               expm = 1
            ELSE
               expm = 2
            END IF
            IF (expm .ne. 0) y_vmec = y_vmec*(fac1**expm)
            IF (mp .le. 2) y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
         END IF

!
!      3. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!         (spline routines in LIBSTELL Miscel folder)
!
         yp1 = -1.e30_dp;  ypn = -1.e30_dp
         CALL spline (snodes_vmec, y_vmec, nsv, yp1, ypn, y2_vmec)

!
!      4. Interpolate onto snodes mesh
!
         DO js = 1, nsi
            CALL splint (snodes_vmec, y_vmec, y2_vmec, nsv,              &
                         snodes(js), y_spline(js))
         END DO

         IF (istat.eq.0 .and. mp.gt.0) THEN
            IF (expm .ne. 0) y_spline = y_spline*(fac2**expm)
         END IF

         ymn_spline(modes,:) = y_spline(:)


!IDENTITY TEST THAT INTERPOLATION WORKS
         IF (nsi == nsv .and. .not.lscale) THEN
            IF (ANY(ABS(ymn_vmec(modes,2:) - ymn_spline(modes,2:)).GT.1.E-10_dp)) THEN
               PRINT *,'MODE: ',modes,'y_vmec .ne. y_spline'
            ENDIF
         END IF
!
!        PRINT OUT FOR CHECKING
!
!         IF (xn_vmec(modes) .eq. 0) THEN
!             WRITE (36, *) mp, INT(xn_vmec(modes))
!             DO js = 1, nsv/10
!               WRITE (36, '(1p2e14.4)') snodes_vmec(js), ymn_vmec(modes,js)
!            END DO
!             WRITE (36, *)
!            DO js = 1, nsi/3
!               WRITE (36, '(1p2e14.4)') snodes(js), ymn_spline(modes,js)
!            END DO
!             WRITE (36, *)
!         END IF

      END DO

      istat = 0

      END SUBROUTINE O_Spline_Fourier_Modes
      

      SUBROUTINE Spline_OneD_Array (y_vmec, y_spline, istat)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT)     :: istat
      REAL(dp), DIMENSION(ns_vmec), INTENT(IN)  :: y_vmec
      REAL(dp), DIMENSION(ns_i)   , INTENT(OUT) :: y_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER          :: one = 1
      INTEGER                         :: js, modes, ntype, mp
      REAL(dp), DIMENSION(ns_vmec) :: snodes_vmec, y2_vmec
      REAL(dp), DIMENSION(ns_i)    :: snodes, fac2
      REAL(dp)                     :: hs_vmec, yp1, ypn
!-----------------------------------------------
!
!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
!

!
!     1. Set up "knots" on initial (vmec, svmec ~ phi) mesh 
!
      IF (ns_vmec .le. 1) THEN
         istat = 1
         RETURN
      END IF

      hs_vmec = one/(ns_vmec-1)
      DO js = 1, ns_vmec
         snodes_vmec(js) = hs_vmec*(js-1)
      END DO

!
!     2. Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!        if s_siesta ~ sfinal, this is then the original vmec mesh
!
      IF (ns_i .LE. 1) THEN
         istat = 2
         RETURN
      END IF

      ohs_i = ns_i-1
      hs_i = one/ohs_i
      DO js = 1, ns_i
         fac2(js) = hs_i*(js-1)
         IF (iScale .EQ. 0) THEN
            snodes(js) = fac2(js)*fac2(js)        !polar mesh, s-siesta~s-vmec**2
         ELSE
            snodes(js) = fac2(js)                 !vmec mesh   s-siesta~s-vmec
         END IF
      END DO



!     4. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!
      yp1 = -1.e30_dp;  ypn = -1.e30_dp
      CALL spline (snodes_vmec, y_vmec, ns_vmec, yp1, ypn, y2_vmec)

!
!     5. Interpolate onto snodes mesh
!
      DO js = 1, ns_i
         CALL splint (snodes_vmec, y_vmec, y2_vmec, ns_vmec,             &
                      snodes(js), y_spline(js))
      END DO


      END SUBROUTINE Spline_OneD_Array

      
      SUBROUTINE Add_Ghost_Points(ymn_vmec, ymn_ghost)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(in)            :: ymn_vmec(mnmax,ns_vmec)
      REAL(dp), INTENT(out)           :: ymn_ghost(mnmax,ns_vmec+1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER             :: p5 = 0.5_dp
      INTEGER                            :: modes, js, mp
      REAL(dp), DIMENSION(ns_vmec)    :: y_vmec
!-----------------------------------------------

      ymn_ghost(:,2:ns_vmec) = ymn_vmec(:,2:ns_vmec)
!
!     ADDS GHOST POINTS AT js=1 (s=0) AND js=ns+1 (s=1)
!
      ymn_ghost(:,ns_vmec+1) = 2*ymn_vmec(:,ns_vmec) - ymn_vmec(:,ns_vmec-1)

      DO modes = 1, mnmax
         mp = xm_vmec(modes)
         IF (mp .eq. 0) THEN
            ymn_ghost(modes,1) = 2*ymn_vmec(modes,2)-ymn_vmec(modes,3)
         ELSE
            ymn_ghost(modes,1) = 0
         END IF
      END DO

      END SUBROUTINE


      SUBROUTINE LoadRZL_VMEC (istat)
      USE fourier, ONLY: toijsp, init_fourier, f_cos, f_sin,             &
                         f_none, f_sum, f_du, f_dv
      USE island_params, ns=>ns_i
    USE island_params, mpol=>mpol_i, ntor=>ntor_i
!      USE shared_functions, ONLY: DUMPOUT3
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT)     :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER    ::  m, n, u, mn, js
      REAL(dp), DIMENSION(ns_vmec)  :: fac1, fac2
      REAL(dp) :: ohs, fac=1
      REAL(dp), PARAMETER :: twopi = 6.283185308
!-----------------------------------------------
!
!     1. LOAD FOURIER FACTORS USING VMEC mpol,ntor AND ISLAND nu,nv
!        IDEA: USE VMEC mpol, ntor VALUES AND CALL fixarray, FLIP VMEC n -> -n in REPACK
!     2. COMPUTE METRIC ELEMENTS, JACOBIAN, LAMBDA ON ISLAND MESH
      
      CALL init_fourier

!
!     COMPUTE R, Z (Metric elements, jacobian), LAMBDA (For initializing B)
!
!     FIRST, REPACK SPLINED ARRAYS FROM VMEC-ORDERING TO SIESTA-ORDERING
!     AND SWAP THE N->-N MODES TO BE CONSISTENT WITH mu+nv SIESTA ARGUMENT
!


!      DO mn=1, mnmax
!        DO js=1,ns
!          WRITE(777,*) rmnc_vmec(mn,js), rmnc_spline(mn,js)
!          CALL FLUSH(777)
!        END DO 
!      END DO

      CALL repack(rmnc_i, zmns_i, lmns_i, jcurrumnc, jcurrvmnc,         &
                  rmnc_spline, zmns_spline, lmns_spline,                &
                  currumnc_spline, currvmnc_spline, isym)
      IF (lasym)                                                        &
      CALL repack(rmns_i, zmnc_i, lmnc_i, jcurrumns, jcurrvmns,         &                                                     
                  rmns_spline, zmnc_spline, lmnc_spline,                &
                  currumns_spline, currvmns_spline, isym)

!
!     COMPUTE AND STORE R, Z AND THEIR ANGULAR DERIVATIVES, AND LAMBDA (NEED FOR VECTOR POT)
!
      ALLOCATE (r1_i(nu_i, nv_i, ns_i), z1_i(nu_i, nv_i, ns_i),         &
                ru_i(nu_i, nv_i, ns_i), zu_i(nu_i, nv_i, ns_i),         &
                rv_i(nu_i, nv_i, ns_i), zv_i(nu_i, nv_i, ns_i),         &
                stat = istat)

!      IF (istat .NE. 0) CALL ErrorAbort('Allocation failed in LoadRZL_VMEC')

      CALL toijsp (rmnc_i, r1_i, f_none, f_cos)
      CALL toijsp (rmnc_i, ru_i, f_du,   f_cos)
      CALL toijsp (rmnc_i, rv_i, f_dv,   f_cos)

      CALL toijsp (zmns_i, z1_i, f_none, f_sin)
      CALL toijsp (zmns_i, zu_i, f_du,   f_sin)
      CALL toijsp (zmns_i, zv_i, f_dv,   f_sin)

      IF (.NOT.lasym) GOTO 1000

      m = IOR(f_sum, f_du);  n = IOR(f_sum, f_dv)
      CALL toijsp (rmns_i, r1_i, f_sum, f_sin)
      CALL toijsp (rmns_i, ru_i, m,     f_sin)
      CALL toijsp (rmns_i, rv_i, n,     f_sin)

      CALL toijsp (zmnc_i, z1_i, f_sum, f_cos)
      CALL toijsp (zmnc_i, zu_i, m,     f_cos)
      CALL toijsp (zmnc_i, zv_i, n,     f_cos)

 1000 CONTINUE

      rmax = MAXVAL(r1_i);      rmin = MINVAL(r1_i)
      zmax = MAXVAL(z1_i);      zmin = MINVAL(z1_i)
      IF (.NOT.lasym) zmin = MIN(zmin, -zmax)                           !ONLY TOP HALF

#if defined(NOSKS)
!
!     COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN
!
      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)

!
!     COMPUTE FULL-MESH METRIC ELEMENTS
      CALL full_mesh_metrics
#endif
!
!     CLEAN-UP EXTRA ARRAYS
!
!      DEALLOCATE (r1_i, z1_i, ru_i, zu_i, rv_i, zv_i, stat=istat)
      

      END SUBROUTINE LoadRZL_VMEC


      SUBROUTINE repack (rmn, zmn, lmn, currumn, currvmn,               &
                         rmn_spl, zmn_spl, lmn_spl,                     &
                         currumn_spl, currvmn_spl, ipar)
      USE fourier, ONLY: orthonorm
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rmn, zmn, lmn, currumn, currvmn
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) ::   rmn_spl, zmn_spl,      &
               lmn_spl, currumn_spl, currvmn_spl
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER              :: modes, js, m, n, ntype, n1, s1, s2, istat
!-----------------------------------------------
!
!     The splined arrays (rmnc_spl, etc) are all in VMEC ordering
!
!         rmnc_spl(1:mnmax,ns_i)
!
!     Now pack them so they can be used by ISLAND Fourier routines
!
!         rmn(0:mpol, -ntor:ntor, ns_i)
!
!     NOTE: mpol_i == mpol_vmec, ntor_i == ntor_vmec here. Also, must
!           FLIP the sign of 'n' terms since VMEC uses mu-nv and SIESTA
!           uses mu+nv for args of cos and sin
!
!     KEEP lmn ON HALF MESH AND ADD EXTRAPOLATION POINTS AT 1 (rho=0)
!          AND NS+1 (=1)!
!
      CALL ASSERT(SIZE(rmn_spl,1).EQ.mnmax,'rmn_spl wrong size1 in REPACK')
      CALL ASSERT(SIZE(rmn_spl,2).EQ.ns_i,'rmn_spl wrong size2 in REPACK')
      CALL ASSERT(SIZE(lmn_spl,2).EQ.ns_i+1,'lmn_spl wrong size2 in REPACK')

      IF (ipar .EQ. isym) THEN
         s1 = 1
      ELSE
         s2 = 1
      END IF

      ALLOCATE(rmn(0:mpol_i,-ntor_i:ntor_i,ns_i),                       & 
               zmn(0:mpol_i,-ntor_i:ntor_i,ns_i),                       &
               lmn(0:mpol_i,-ntor_i:ntor_i,ns_i+1),                     &
               currumn(0:mpol_i,-ntor_i:ntor_i,ns_i),                   &
               currvmn(0:mpol_i,-ntor_i:ntor_i,ns_i), stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation error in REPACK')
      rmn = 0;  zmn = 0;  lmn = 0;  currumn = 0;  currvmn = 0

      IF (iam.EQ.0 .AND. lverbose .AND. ipar.EQ.isym) THEN
         m = INT(MAXVAL(xm_vmec))
         n = INT(MAXVAL(ABS(xn_vmec)))/nfp_vmec
         IF (m.GT.mpol_i .OR. n.GT.ntor_i)                              &
            PRINT 100,'It is recommended to increase the number ' //    &
                 'of modes to at least VMEC values mpol=',m,' ntor=',n
 100  FORMAT(1x,a,i4,a,i4)
      ENDIF
!
!     FLIP VMEC n -> SIESTA -n SO TRIG ARG IS mu+nv IN SIESTA
!     SET iflipj -> -1 FOR POSITIVE JACOBIAN SYSTEM
!
      DO modes = 1,mnmax
         m = xm_vmec(modes)
         n = xn_vmec(modes)/nfp_vmec
         n1 = ABS(n)
         IF (m.GT.mpol_i .OR. n1.GT.ntor_i) CYCLE
!
!     LOAD n>=0 ONLY FOR M=0
!
         IF (m .EQ. 0) THEN
            IF (ipar .EQ. isym)  s2 = -SIGN(1,n)
            IF (ipar .EQ. iasym) s1 = -SIGN(1,n)
            rmn(m,n1,:) = rmn(m,n1,:)                                   &
                        + s1*rmn_spl(modes,:)
            zmn(m,n1,:) = zmn(m,n1,:)                                   &
                        + s2*zmn_spl(modes,:)
            lmn(m,n1,:) = lmn(m,n1,:)                                   &   
                        + s2*lmn_spl(modes,:)*iflipj
            currumn(m,n1,:) = currumn(m,n1,:)                           &
                        + s1*currumn_spl(modes,:)
            currvmn(m,n1,:) = currvmn(m,n1,:)                           &
                        + s1*currvmn_spl(modes,:)
                       
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
      DO js = 1, ns_i+1
         lmn(:,:,js) = lmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         IF (js .GT. ns_i) EXIT
      
         rmn(:,:,js) = rmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         zmn(:,:,js) = zmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         currumn(:,:,js) = 2*hs_i*(js-1)*currumn(:,:,js)                &
                         / orthonorm(0:mpol_i,-ntor_i:ntor_i)
         currvmn(:,:,js) = 2*hs_i*(js-1)*currvmn(:,:,js)                &
                         / orthonorm(0:mpol_i,-ntor_i:ntor_i)
      END DO

      DEALLOCATE (rmn_spl, zmn_spl, lmn_spl, currumn_spl,               &
                  currvmn_spl, stat=istat)

      END SUBROUTINE repack


      SUBROUTINE half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
      USE island_params, ns=>ns_i, nuv=>nuv_i, nu=>nu_i, nv=>nv_i
      IMPLICIT NONE  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(nuv,ns), INTENT(in) ::                        &
         r1_i, ru_i, rv_i, z1_i, zu_i, zv_i
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 1.d0/2.d0, zero = 0
      INTEGER             :: js, lk, js1, istat, u, v, imin(2), imax(2)
      REAL(dp)            :: mintest, maxtest
      REAL(dp), DIMENSION(nuv) ::                                       &
                   r12, ru12, rv12, zu12, zv12, rs12, zs12
!-----------------------------------------------
!
!     ALLOCATE METRIC ELEMENT ARRAYS
!
      ALLOCATE(sqrtg(nuv,ns),                                           &
               gss(nuv,ns), gsu(nuv,ns), gsv(nuv,ns),                   &
               guu(nuv,ns), guv(nuv,ns), gvv(nuv,ns), stat=istat)
!      IF (istat .ne. 0) CALL ErrorAbort('ALLOCATION ERROR1 IN HALF_MESH_METRICS')

!     COMPUTE ALL ON THE HALF MESH
      sqrtg(:,1) = 0

      DO js = 2, ns
         js1 = js-1
         r12 = p5*(r1_i(:,js) + r1_i(:,js1))
         rs12 = (r1_i(:,js) - r1_i(:,js1))*ohs_i
         ru12= p5*(ru_i(:,js) + ru_i(:,js1))
         rv12= p5*(rv_i(:,js) + rv_i(:,js1))
         zs12 = (z1_i(:,js) - z1_i(:,js1))*ohs_i
         zu12= p5*(zu_i(:,js) + zu_i(:,js1))
         zv12= p5*(zv_i(:,js) + zv_i(:,js1))
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
         IF (iam .EQ. 0) PRINT 100,' MIN-G: ', mintest,' AT JS: ',      &
                       imin(2),'  MAX-G: ', maxtest,' AT JS: ', imax(2)
         CALL ASSERT(mintest*maxtest.GT.zero,' Jacobian changed sign in half_mesh_metrics')
 100  FORMAT(2(a,1p,e12.4,a,i5))
      END IF

      END SUBROUTINE half_mesh_metrics


      SUBROUTINE full_mesh_metrics
      USE island_params, ns=>ns_i, nuv=>nuv_i, nzeta=>nv_i, ntheta=>nu_i
      IMPLICIT NONE  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER:: istat, js
      REAL(dp)                                :: maxtest, eps, r1, s1, t1
      REAL(dp)                                :: det(nuv)
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: temps, tempu, tempv
!-----------------------------------------------
!
!     This subroutine gets the lower metric elements on the full mesh
!     Preserves positive definiteness of metric tensor
!
      ALLOCATE(gssf(nuv,ns), guuf(nuv,ns), gvvf(nuv,ns),                 &
               gsuf(nuv,ns), gsvf(nuv,ns), guvf(nuv,ns),                 &
               hss(nuv,ns),  huu(nuv,ns),  hvv(nuv,ns),                  & 
               hsu(nuv,ns),  hsv(nuv,ns),  huv(nuv,ns), stat=istat)

!      IF (istat .NE. 0) CALL ErrorAbort('ALLOCATION ERROR IN FULL_MESH_METRICS')
      !CALL to_full_mesh(gss,gssf,nuv)
      !CALL to_full_mesh(guu,guuf,nuv)
      !CALL to_full_mesh(gvv,gvvf,nuv)
      !CALL to_full_mesh(gsu,gsuf,nuv)
      !CALL to_full_mesh(gsv,gsvf,nuv)
      !CALL to_full_mesh(guv,guvf,nuv)

      CALL to_full_mesh(gss,gssf)
      CALL to_full_mesh(guu,guuf)
      CALL to_full_mesh(gvv,gvvf)
      CALL to_full_mesh(gsu,gsuf)
      CALL to_full_mesh(gsv,gsvf)
      CALL to_full_mesh(guv,guvf)
      guuf(:,1) = 0;  gsuf(:,1) = 0;  guvf(:,1) = 0

!
!     Compute upper metric elements (inverse of lower matrix) on full mesh
!     NOTE: DET == 1/det{|Gij| == 1/sqrtg(full)**2
!
      hss(:,1) = 0;  hsu(:,1) = 0;  hsv(:,1) = 0
      huu(:,1) = 0;  huv(:,1) = 0;  hvv(:,1) = 0

      DO js = 2, ns
         det(:) =                                                        &
                gssf(:,js)*(guuf(:,js)*gvvf(:,js)-guvf(:,js)*guvf(:,js)) &
              + gsuf(:,js)*(guvf(:,js)*gsvf(:,js)-gsuf(:,js)*gvvf(:,js)) &
              + gsvf(:,js)*(gsuf(:,js)*guvf(:,js)-gsvf(:,js)*guuf(:,js))
         IF (ANY(det .LE. zero)) THEN
            IF (iam .EQ. 0) PRINT *,'Determinant |gijf| <= 0 at js = ',js
            IF (js .EQ. ns) THEN
               det(:) = sqrtg(:,ns)**2
               gssf(:,ns) = gss(:,ns)
               guuf(:,ns) = guu(:,ns)
               gvvf(:,ns) = gvv(:,ns)
               gsuf(:,ns) = gsu(:,ns)
               guvf(:,ns) = guv(:,ns)
               gsvf(:,ns) = gsv(:,ns)
            ELSE
               det(:) = ((sqrtg(:,js)+sqrtg(:,js+1))/2)**2
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

!
!     CONFIRM ACCURACY OF INVERSE
!
      eps = 0.01_dp*SQRT(EPSILON(eps))

      ALLOCATE(temps(ntheta,nzeta,ns),                                   &
               tempu(ntheta,nzeta,ns),                                   &
               tempv(ntheta,nzeta,ns), stat=istat)

!     gijf * hj1   (first column of hij matrix)
      CALL tolowerf(hss, hsu, hsv, temps, tempu, tempv, 1, ns)
      r1 = MAXVAL(ABS(temps(:,:,2:)-1))
      s1 = MAXVAL(ABS(tempu(:,:,2:)))
      t1 = MAXVAL(ABS(tempv(:,:,2:)))
      maxtest = MAX(r1,s1,t1)
!      IF (maxtest.GT.eps .AND. iam.EQ.0) CALL ErrorAbort('ERROR1 IN FULL_MESH_METRICS')

!     gijf * hj2   (second column of hij matrix)
      CALL tolowerf(hsu, huu, huv, temps, tempu, tempv, 1, ns)
      r1 = MAXVAL(ABS(temps(:,:,2:)))
      s1 = MAXVAL(ABS(tempu(:,:,2:)-1))
      t1 = MAXVAL(ABS(tempv(:,:,2:)))
      maxtest = MAX(r1,s1,t1)
!      IF (maxtest.GT.eps .AND. iam.EQ.0) CALL ErrorAbort('ERROR2 IN FULL_MESH_METRICS')

!     gijf * hj3   (third column of hij matrix)
      CALL tolowerf(hsv, huv, hvv, temps, tempu, tempv, 1, ns)
      r1 = MAXVAL(ABS(temps(:,:,2:)))
      s1 = MAXVAL(ABS(tempu(:,:,2:)))
      t1 = MAXVAL(ABS(tempv(:,:,2:)-1))
      maxtest = MAX(r1,s1,t1)
!      IF (maxtest.GT.eps .AND. iam.EQ.0) CALL ErrorAbort('ERROR3 IN FULL_MESH_METRICS')
      
      DEALLOCATE(temps, tempu, tempv)

      END SUBROUTINE full_mesh_metrics

      SUBROUTINE toupper(xsubsij, xsubuij, xsubvij,                     &
                         xsupsij, xsupuij, xsupvij, nsmin, nsmax)
        USE stel_kinds
        USE island_params, nuv=>nuv_i
!
!     This subroutine converts from lower to upper components on the FULL mesh
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER :: istat, nsmin, nsmax
        REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(IN)  ::            &
        xsubsij, xsubuij, xsubvij 
        REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(OUT) ::            &
        xsupsij, xsupuij, xsupvij        
!-----------------------------------------------
        INTEGER :: js
!-----------------------------------------------

        DO js=nsmin, nsmax
        xsupsij(:,js) = hss(:,js)*xsubsij(:,js)                         &
                      + hsu(:,js)*xsubuij(:,js)                         &
                      + hsv(:,js)*xsubvij(:,js)
        xsupuij(:,js) = hsu(:,js)*xsubsij(:,js)                         &
                      + huu(:,js)*xsubuij(:,js)                         &
                      + huv(:,js)*xsubvij(:,js)
        xsupvij(:,js) = hsv(:,js)*xsubsij(:,js)                         &
                      + huv(:,js)*xsubuij(:,js)                         &
                      + hvv(:,js)*xsubvij(:,js)       
        END DO

      END SUBROUTINE toupper


      SUBROUTINE tolowerh(xsupsij, xsupuij, xsupvij,                      &
                          xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: nuv=>nuv_i
      IMPLICIT NONE        
!
!     Moves contravariant (upper) to covariant (lower) magnetic field
!     components on the half radial mesh.
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!----------------------------------------------
      INTEGER, INTENT(IN) :: nsmin, nsmax  
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(IN)  ::               &
                   xsupsij, xsupuij, xsupvij 
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(OUT) ::               &
                   xsubsij, xsubuij, xsubvij
!-----------------------------------------------
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
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE tolowerf(xsupsij, xsupuij, xsupvij,                     &
                          xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: ns=>ns_i, nuv=>nuv_i
      IMPLICIT NONE        
!
!     Computes covariant (lower) from contravariant (upper) magnetic
!     field components on the full mesh. 
!
!     xsubsij = gssf*xsupsij + gsuf*xsupuij + gsvf*xsupvij
!     xsubuij = gsuf*xsupsij + guuf*xsupuij + guvf*xsupvij
!     xsubvij = gsvf*xsupsij + guvf*xsupuij + gvvf*xsupvij       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: nsmin, nsmax
      REAL(dp), DIMENSION(nuv,ns), INTENT(IN)  ::                        &
                xsupsij, xsupuij, xsupvij 
      REAL(dp), DIMENSION(nuv,ns), INTENT(OUT) ::                        &
                xsubsij, xsubuij, xsubvij
!-----------------------------------------------
      INTEGER :: js, js1
!-----------------------------------------------
      DO js=nsmin,nsmax
      js1 = js-nsmin+1
      xsubsij(:,js1) = gssf(:,js)*xsupsij(:,js1)                         &
                     + gsuf(:,js)*xsupuij(:,js1)                         &
                     + gsvf(:,js)*xsupvij(:,js1) 

      xsubuij(:,js1) = gsuf(:,js)*xsupsij(:,js1)                         &
                     + guuf(:,js)*xsupuij(:,js1)                         &
                     + guvf(:,js)*xsupvij(:,js1) 

      xsubvij(:,js1) = gsvf(:,js)*xsupsij(:,js1)                         &
                     + guvf(:,js)*xsupuij(:,js1)                         &
                     + gvvf(:,js)*xsupvij(:,js1)       
      END DO

      END SUBROUTINE tolowerf


      SUBROUTINE cleanup_metric_elements
      USE fourier, ONLY: dealloc_fourier
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER         :: istat
!-----------------------------------------------
    
!     Note: sqrtg is deallocated in init_bcovar and stored in jacob variable

      DEALLOCATE(gss, gsu, gsv, guu, guv, gvv,                           &
                 hss, hsu, hsv, huu, huv, hvv,                           &
                 phipf_i, chipf_i, presf_i, stat=istat)
      
      DEALLOCATE (r1_i, z1_i, ru_i, zu_i, rv_i, zv_i)

      CALL dealloc_fourier
      DEALLOCATE(rmnc_i, zmns_i)     ! Used by output subroutines.        
      IF (lasym) DEALLOCATE(rmns_i, zmnc_i)


      END SUBROUTINE cleanup_metric_elements

       
       SUBROUTINE dealloc_full_lower_metrics
         USE stel_kinds
         IMPLICIT NONE
         INTEGER :: istat

         DEALLOCATE(gssf, guuf, gvvf, gsuf, gsvf, guvf, stat = istat)
!         IF (istat .NE. 0) CALL ErrorAbort('Problem dealloc. in LOWER_METRIC_FULL')

       END SUBROUTINE dealloc_full_lower_metrics

      END MODULE metrics
