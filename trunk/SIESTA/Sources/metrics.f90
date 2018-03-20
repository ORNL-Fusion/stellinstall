!>
!!  \brief Module for reading in VMEC input and computing metric elements
!!   on the SIESTA (sqrt-flux) radial grid.
!<
      MODULE metrics
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants
      USE island_params
      USE shared_data, ONLY: lasym, lverbose, unit_out
      USE read_wout_mod, ns_vmec=>ns, mpol_vmec=>mpol, ntor_vmec=>ntor, &
          rmnc_vmec=>rmnc, zmns_vmec=>zmns, lmns_vmec=>lmns,            &
          xm_vmec=>xm, xn_vmec=>xn, chipf_vmec=>chipf,                  &  ! MRC 4/1/2016
          rmns_vmec=>rmns, zmnc_vmec=>zmnc, lmnc_vmec=>lmnc,            &  ! MRC 12/1/2016
          phipf_vmec=>phipf, presf_vmec=>presf, nfp_vmec=>nfp,          &
          wb_vmec=>wb, wp_vmec=>wp, gamma_vmec=>gamma,                  &
          volume_vmec=>volume, raxis_vmec=>raxis, lasym_vmec=>lasym,    &
          assert_wout=>assert, iasym_vmec=>iasym
      USE descriptor_mod, ONLY: iam
      USE timer_mod
      USE utilities, ONLY: to_full_mesh
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
      REAL(dp), DIMENSION(:,:), ALLOCATABLE ::                          &
                   sqrtg,                                               & !sqrt(g): Jacobian on half grid
                   gss, gsu, gsv, guu, guv, gvv,                        & !elements of lower metric tensor (half mesh)
                   hss, hsu, hsv, huu, huv, hvv                            !elements of upper metric tensor (full mesh)

      REAL(dp) :: rmax, rmin, zmax, zmin

      REAL(dp), DIMENSION(:,:), ALLOCATABLE:: gssf, guuf, gvvf,         & !Lower Metric elements (full)
                   gsuf, gsvf, guvf   

      INTEGER :: iflipj=1                                                 !=-1 to flip sign of vmec jac
!
!     LOCAL (PRIVATE) HELPER ROUTINES (NOT ACCESSIBLE FROM OUTSIDE THIS MODULE)
!
      PRIVATE spline_fourier_modes, loadrzl_vmec, add_ghost_points,     &
	          getrz, getrzlsplines
      INTEGER, PARAMETER, PRIVATE :: iScale=0                             !=0, s ~ SQRT(PHI) - default; = 1, s ~ PHI (VMEC)
      INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE   ::             &
                r1_i, z1_i, ru_i, zu_i, rv_i, zv_i

      REAL(dp), PRIVATE :: skston, skstoff

      CONTAINS

      SUBROUTINE init_metric_elements(ns_in, nsext_in, mpol_in, ntor_in,       &
                                      wout_file)
      USE timer_mod, ONLY: init_timers
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*)            :: wout_file
      INTEGER, INTENT(IN)      :: ns_in, nsext_in, mpol_in, ntor_in
      INTEGER                  :: istat, imesh, js
      REAL(dp)                 :: t1, t2
      INTEGER                  :: m, n, mn
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

      IF (iScale.NE.0 .AND. lverbose)                                   &
         PRINT *,'RESET ISCALE=0 AFTER TESTING!'

!
!     READ-IN DATA FROM VMEC WOUT FILE (LIBSTELL ROUTINE)
!
#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL read_wout_file(wout_file, istat)
#if defined(SKS) 
      CALL second0(skstoff)
      read_wout_file_time=read_wout_file_time+(skstoff-skston)
#endif
      CALL ASSERT(istat.EQ.0,'Read-wout error in INIT_METRIC_ELEMENTS')

!  Asymmetric terms may be overriden in the namelist input file.
!      lasym = lasym .OR. lasym_vmec

!      DO js = 2, ns_in
!         WRITE (30, '(i4, 1pE10.2)') js, SQRT(SUM(currvmnc(:,js)**2))
!      END DO

      nfp_i = nfp_vmec
      wb_i  = (twopi*twopi)*wb_vmec
      wp_i  = (twopi*twopi)*wp_vmec
      volume_i = volume_vmec

      CALL ASSERT(nfp_i.NE.0, ALLOCATED(raxis_vmec),                    &
         'ERROR INITIALIZING FROM WOUT_FILE IN METRICS!')

      rmajor_i = raxis_vmec(0,1)

      CALL ASSERT(wb_i.gt.zero,'wb_vmec = 0:')
      gnorm_i = wb_i + wp_i/(gamma-1)
      gnorm_i = 1._dp/gnorm_i

      ns_i = ns_in + nsext_in
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

!     UPDATED WOUT TO WRITE OUT CHIPF (SPH) 
!
!      iotaf_vmec = iotaf_vmec*phipf_vmec                 MRC 4/1/2016
!
!     SPLINE R, Z. L FOURIER COMPONENTS in s FROM ORIGINAL VMEC MESH (s ~ phi, ns_vmec points) 
!     TO A "POLAR" MESH [s ~ sqrt(phi), ns_i POINTS] WITH BETTER AXIS RESOLUTION
!
      CALL GetRZLSplines(ns_in, rmnc_vmec, zmns_vmec, lmns_vmec,               &
	                     rmnc_spline, zmns_spline, lmns_spline, isym)
      IF (lasym) THEN
 	      CALL GetRZLSplines(ns_in, rmns_vmec, zmnc_vmec, lmnc_vmec,           &
		                     rmns_spline, zmnc_spline, lmnc_spline, iasym)
      END IF

!
!     Spline 1-D arrays: careful -> convert phipf VMEC and multiply
!     by ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island
!
      ALLOCATE(phipf_i(ns_in), chipf_i(ns_in), presf_i(ns_in), stat=istat)

      CALL second0(skston)
!      CALL Spline_OneD_Array (iotaf_vmec, chipf_i, istat)                       MRC 4/1/2016
      CALL Spline_OneD_Array (chipf_vmec, chipf_i, istat)                      ! MRC 4/1/2016
      CALL Spline_OneD_Array (phipf_vmec, phipf_i, istat)
      presf_vmec = mu0 * presf_vmec
      CALL Spline_OneD_Array (presf_vmec, presf_i, istat)
      CALL second0(skstoff)
#if defined(SKS) 
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
         DO js = 1, ns_in
            phipf_i(js) = 2 * hs_i*(js-1) * phipf_i(js)
            chipf_i(js) = 2 * hs_i*(js-1) * chipf_i(js)
         END DO
      END IF
121   FORMAT(a,i4, 4(a,1p,e12.3))

      t1 = twopi*hs_i*(SUM(phipf_i(2:ns_in)) + SUM(phipf_i(1:ns_in-1)))/2
      IF (iam .EQ. 0) THEN
         WRITE (unit_out, 100) volume_i,                                &
            1.E-6_dp*(wb_i+wp_i/(gamma-1))/mu0,                         &
            1.E-6_dp*wp_i/mu0,                                          &
            t1,wp_i/wb_i,chipf_i(2)/phipf_i(2),                         &
            chipf_i(ns_in)/phipf_i(ns_in), gamma
      ENDIF

 100  FORMAT(/,' INITIAL PARAMETERS (FROM VMEC)',/,                     &
               ' PHYSICS QUANTITIES',/,                                 &
               ' PLASMA VOLUME (M^3): ',1pE12.4,                        &
               ' TOTAL MHD ENERGY (MJ):   ', 1pE16.8,/,                 &
               ' THERMAL ENERGY (MJ): ', 1pE12.4,                       &
               ' EDGE TOROIDAL FLUX (Wb): ', 1pE12.4,/,                 &
               ' DIMENSIONLESS QUANTITIES',/,                           &
               ' <BETA>: ',1pE12.4, ' IOTA(0): ',1pE12.4,               &
               ' IOTA(1) ',1pE12.4, ' GAMMA: ',1pE12.4,/, 21('-'),/)

!
!     CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON SQRT(FLUX) - "POLAR" - MESH
!     AND COMPUTE METRIC ELEMENTS AND JACOBIAN
!
      CALL second0(skston)
      CALL LoadRZL_VMEC(istat)
      CALL second0(skstoff)
#if defined(SKS) 
      LoadRZL_VMEC_time=LoadRZL_VMEC_time+(skstoff-skston)
#endif
      CALL ASSERT(istat.EQ.0,'LoadRZL error in INIT_METRIC_ELEMENTS')

!      IF (lwout_opened) CALL read_wout_deallocate  !  Don't deallocate the wout file here.

      END SUBROUTINE init_metric_elements


      SUBROUTINE GetRZLSplines(ns, rmn_vmec,   zmn_vmec,   lmn_vmec,           &
	                           rmn_spline, zmn_spline, lmn_spline, ipar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar, ns
      REAL(dp), DIMENSION(mnmax, ns_vmec), INTENT(IN) ::                &
                rmn_vmec, zmn_vmec, lmn_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rmn_spline, zmn_spline,  &
	                                           lmn_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat, ntype
      REAL(dp), ALLOCATABLE     :: lmn_half(:,:)
!-----------------------------------------------
      ALLOCATE(rmn_spline(mnmax,ns), zmn_spline(mnmax,ns),              &
               lmn_spline(mnmax,ns+1), stat=istat)
      CALL ASSERT(istat.EQ.0, 'Allocation error in GetRZLSplines')

      IF (ipar.EQ.iasym .AND. .NOT.lasym_vmec) THEN
         rmn_spline = 0; zmn_spline = 0; lmn_spline = 0
         RETURN
      END IF
	     
      DO ntype = 1, 3
         IF (ntype .EQ. 1) THEN
            istat = 0
            CALL second0(skston)
            CALL Spline_Fourier_Modes(rmn_vmec, rmn_spline, ns_vmec, ns, istat)
            CALL second0(skstoff)
#if defined(SKS) 
            Spline_Fourier_Modes_time=Spline_Fourier_Modes_time+(skstoff-skston)
#endif
            CALL ASSERT(istat.EQ.0,'Error splining rmnc')
         ELSE IF (ntype .EQ. 2) THEN
            istat = 0
            CALL second0(skston)
            CALL Spline_Fourier_Modes(zmn_vmec, zmn_spline, ns_vmec, ns, istat)
            CALL second0(skstoff)
#if defined(SKS) 
            Spline_Fourier_Modes_time=Spline_Fourier_Modes_time+(skstoff-skston)
#endif
            CALL ASSERT(istat.EQ.0,'Error splining zmns')
         ELSE 
            ALLOCATE(lmn_half(mnmax,ns_vmec+1), stat=istat)
            CALL second0(skston)
            CALL Add_Ghost_Points(lmn_vmec, lmn_half)
            CALL second0(skstoff)
#if defined(SKS) 
            Add_Ghost_Points_time=Add_Ghost_Points_time+(skstoff-skston)
#endif
            istat = 0
            CALL second0(skston)
            CALL Spline_Fourier_Modes(lmn_half, lmn_spline,           &
                                      ns_vmec+1, ns+1, istat)
            CALL second0(skstoff)
#if defined(SKS) 
            Spline_Fourier_Modes_time=Spline_Fourier_Modes_time+(skstoff-skston)
#endif
            CALL ASSERT(istat.EQ.0,'Error splining lmns')
            DEALLOCATE(lmn_half)
         END IF
      END DO
     

	  END SUBROUTINE GetRZLSplines


      
      SUBROUTINE Spline_Fourier_Modes(ymn_vmec, ymn_spline, nsv, nsi,   &
                                      istat)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)        :: nsv, nsi
      INTEGER, INTENT(INOUT)     :: istat
      REAL(dp), DIMENSION(mnmax, nsv), TARGET,                          &
                   INTENT(IN)  :: ymn_vmec
      REAL(dp), DIMENSION(mnmax, nsi), TARGET,                          &
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
      IF (ns_vmec .LE. 1) THEN
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
         CALL splint (snodes_vmec, y_vmec, y2_vmec, ns_vmec,            &
                      snodes(js), y_spline(js))
      END DO


      END SUBROUTINE Spline_OneD_Array

      
      SUBROUTINE Add_Ghost_Points(ymn_vmec, ymn_ghost)
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
      USE fourier, ONLY: init_fourier
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT)     :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns_vmec)  :: fac1, fac2
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
      CALL repack(rmnc_i, zmns_i, lmns_i,                               &
	              rmnc_spline, zmns_spline, lmns_spline, isym)
	  IF (lasym) CALL repack(rmns_i, zmnc_i, lmnc_i,                    &                                                     
	              rmns_spline, zmnc_spline, lmnc_spline, iasym)
!
!     COMPUTE AND STORE R, Z AND THEIR ANGULAR DERIVATIVES, AND LAMBDA (NEED FOR VECTOR POT)
!
      ALLOCATE (r1_i(nu_i, nv_i, ns_i), z1_i(nu_i, nv_i, ns_i),         &
                ru_i(nu_i, nv_i, ns_i), zu_i(nu_i, nv_i, ns_i),         &
                rv_i(nu_i, nv_i, ns_i), zv_i(nu_i, nv_i, ns_i),         &
                stat = istat)

      CALL ASSERT(istat.EQ.0,'Allocation failed in LoadRZL_VMEC')

      CALL GetRZ(rmnc_i, zmns_i, isym)
	  IF (lasym) CALL GetRZ(rmns_i, zmnc_i, iasym)

      rmax = MAXVAL(r1_i);      rmin = MINVAL(r1_i)
      zmax = MAXVAL(z1_i);      zmin = MINVAL(z1_i)
      IF (.NOT.lasym) zmin = MIN(zmin, -zmax)                           !ONLY TOP HALF

!
!     COMPUTE HALF-MESH LOWER METRIC ELEMENTS AND JACOBIAN
      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)

!
!     COMPUTE FULL-MESH LOWER/UPPER METRIC ELEMENTS
      CALL full_mesh_metrics

!
!     CLEAN-UP EXTRA ARRAYS
!
      DEALLOCATE (r1_i, z1_i, ru_i, zu_i, rv_i, zv_i, stat=istat)

      END SUBROUTINE LoadRZL_VMEC

      
	  SUBROUTINE GetRZ(rmn, zmn, ipar)
      USE fourier, ONLY: toijsp, f_cos, f_sin,                          &
                         f_none, f_sum, f_du, f_dv
	  INTEGER, INTENT(IN)  :: ipar
	  REAL(dp), DIMENSION(:,:,:), INTENT(IN) :: rmn, zmn
!-----------------------------------------------
      INTEGER  :: four, fouz, fcomb, mfour, nfour
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         four = f_cos; fouz = f_sin; fcomb = f_none
      ELSE
         four = f_sin; fouz = f_cos; fcomb = f_sum
      END IF

      mfour = IOR(f_du, fcomb);  nfour = IOR(f_dv, fcomb)

      CALL toijsp (rmn, r1_i, fcomb,  four)
      CALL toijsp (rmn, ru_i, mfour,  four)
      CALL toijsp (rmn, rv_i, nfour,  four)

      CALL toijsp (zmn, z1_i, fcomb,  fouz)
      CALL toijsp (zmn, zu_i, mfour,  fouz)
      CALL toijsp (zmn, zv_i, nfour,  fouz)

	  END SUBROUTINE GetRZ

      
	  SUBROUTINE repack (rmn, zmn, lmn, rmn_spl, zmn_spl, lmn_spl, ipar)
      USE fourier, ONLY: orthonorm
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rmn, zmn, lmn
	  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rmn_spl, zmn_spl, lmn_spl
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
               stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation error in REPACK')
      rmn = 0;  zmn = 0;  lmn = 0

      IF (iam.EQ.0 .AND. lverbose .AND. ipar.EQ.isym) THEN
         m = INT(MAXVAL(xm_vmec))
         n = INT(MAXVAL(ABS(xn_vmec)))/nfp_vmec
         IF (m.GT.mpol_i .OR. n.GT.ntor_i)                              & 
            PRINT 100,'It is recommended to increase the number of ' // &
                    'modes to at least VMEC values mpol=',m,' ntor=',n
 100  FORMAT(a,i4,a,i4)
      ENDIF
!
!     FLIP VMEC n -> SIESTA -n SO TRIG ARG IS mu+nv IN SIESTA
!     SET iflipj -> -1 FOR POSITIVE JACOBIAN SYSTEM
!
      DO modes = 1,mnmax
         m = xm_vmec(modes)
         n = xn_vmec(modes)/nfp_vmec
         IF (m.GT.mpol_i .OR. ABS(n).GT.ntor_i) CYCLE
!
!     LOAD n>=0 ONLY FOR M=0
!
         IF (m .EQ. 0) THEN
            n1 = ABS(n)
			IF (ipar .EQ. isym)  s2 = -SIGN(1,n)
			IF (ipar .EQ. iasym) s1 = -SIGN(1,n)
            rmn(m,n1,:) = rmn(m,n1,:)                                   &
                        + s1*rmn_spl(modes,:)
            zmn(m,n1,:) = zmn(m,n1,:)                                   &
                        + s2*zmn_spl(modes,:)
            lmn(m,n1,:) = lmn(m,n1,:)                                   &   
                        + s2*lmn_spl(modes,:)*iflipj
         ELSE
            rmn(m,-n*iflipj,:) = rmn_spl(modes,:)
            zmn(m,-n*iflipj,:) = zmn_spl(modes,:)*iflipj
            lmn(m,-n*iflipj,:) = lmn_spl(modes,:)
         ENDIF
      END DO
!
! RS (10/03/06)  Divide out "orthonorm" factor used in Fourier FFT routines
! 
      DO js = 1, ns_i
         rmn(:,:,js) = rmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         zmn(:,:,js) = zmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         lmn(:,:,js) = lmn(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
      END DO

      DEALLOCATE (rmn_spl, zmn_spl, lmn_spl, stat=istat)

      END SUBROUTINE repack


      SUBROUTINE half_mesh_metrics (r1, ru, rv, z1, zu, zv)
      USE island_params, ns=>ns_i, nuv=>nuv_i
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(nuv,ns), INTENT(IN) ::                        &
         r1, ru, rv, z1, zu, zv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp, zero = 0
      INTEGER             :: js, lk, js1, istat
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
      CALL ASSERT(istat.EQ.0,'ALLOCATION ERROR1 IN HALF_MESH_METRICS')

!     COMPUTE ALL ON THE HALF MESH
      sqrtg(:,1) = 0

      DO js = 2, ns
         js1 = js-1
         r12 = p5*(r1(:,js) + r1(:,js1))
         rs12 = (r1(:,js) - r1(:,js1))*ohs_i
         ru12= p5*(ru(:,js) + ru(:,js1))
         rv12= p5*(rv(:,js) + rv(:,js1))
         zs12 = (z1(:,js) - z1(:,js1))*ohs_i
         zu12= p5*(zu(:,js) + zu(:,js1))
         zv12= p5*(zv(:,js) + zv(:,js1))
         guu(:,js) = ru12*ru12 + zu12*zu12
         guv(:,js) = ru12*rv12 + zu12*zv12
         gvv(:,js) = r12*r12 + rv12*rv12 + zv12*zv12
         gsu(:,js) = rs12*ru12 + zs12*zu12
         gsv(:,js) = rs12*rv12 + zs12*zv12
         gss(:,js) = rs12*rs12 + zs12*zs12
         sqrtg(:,js) = r12*(ru12*zs12 - rs12*zu12)
      END DO

      IF (iScale .EQ. 1) THEN
         sqrtg(:,1) = sqrtg(:,2)
      END IF

!NEED THESE FOR COMPUTING CURRENT AT ORIGIN
      gss(:,1) = gss(:,2);  gsu(:,1) = 0;  gsv(:,1) = gsv(:,2)
      guu(:,1) = 0;         guv(:,1) = 0;  gvv(:,1) = gvv(:,2)

      mintest = MINVAL(sqrtg(:,2:))
      maxtest = MAXVAL(sqrtg(:,2:))

      CALL ASSERT(mintest*maxtest.GT.zero,'Jacobian changed sign in half_mesh_metrics!')

      END SUBROUTINE half_mesh_metrics


      SUBROUTINE full_mesh_metrics
      USE island_params, ns=>ns_i, nuv=>nuv_i, nzeta=>nv_i, ntheta=>nu_i
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
      ALLOCATE(gssf(nuv,ns), guuf(nuv,ns), gvvf(nuv,ns),                &
               gsuf(nuv,ns), gsvf(nuv,ns), guvf(nuv,ns),                &
               hss(nuv,ns),  huu(nuv,ns),  hvv(nuv,ns),                 & 
               hsu(nuv,ns),  hsv(nuv,ns),  huv(nuv,ns), stat=istat)

      CALL ASSERT(istat.EQ.0, 'ALLOCATION ERROR IN FULL_MESH_METRICS')

      CALL to_full_mesh(gss,gssf)
      CALL to_full_mesh(guu,guuf)
      CALL to_full_mesh(gvv,gvvf)
      CALL to_full_mesh(gsu,gsuf)
      CALL to_full_mesh(gsv,gsvf)
      CALL to_full_mesh(guv,guvf)
!      guuf(:,1) = 0;  gsuf(:,1) = 0;  guvf(:,1) = 0                      !leave at gij(2) for inverse 

!
!     Compute upper metric elements (inverse of lower matrix) on full mesh
!     NOTE: DET == 1/det{|Gij| == 1/sqrtg(full)**2
!
      DO js = 1, ns
         det = gssf(:,js)*(guuf(:,js)*gvvf(:,js)-guvf(:,js)*guvf(:,js)) &
             + gsuf(:,js)*(guvf(:,js)*gsvf(:,js)-gsuf(:,js)*gvvf(:,js)) &
             + gsvf(:,js)*(gsuf(:,js)*guvf(:,js)-gsvf(:,js)*guuf(:,js))
         CALL ASSERT(ALL(det.GT.zero),                                  &
                     'Determinant |gijf| <= 0 in full_mesh_metrics')
         det = one/det
         hss(:,js) = det*(guuf(:,js)*gvvf(:,js) - guvf(:,js)*guvf(:,js))
         hsu(:,js) = det*(guvf(:,js)*gsvf(:,js) - gsuf(:,js)*gvvf(:,js))
         hsv(:,js) = det*(gsuf(:,js)*guvf(:,js) - gsvf(:,js)*guuf(:,js))
         huu(:,js) = det*(gssf(:,js)*gvvf(:,js) - gsvf(:,js)*gsvf(:,js))
         huv(:,js) = det*(gsvf(:,js)*gsuf(:,js) - gssf(:,js)*guvf(:,js))
         hvv(:,js) = det*(gssf(:,js)*guuf(:,js) - gsuf(:,js)*gsuf(:,js))
      END DO
!
!     CHECK ACCURACY OF INVERSE (toupper -> tolowerf yields identity)
!
!#undef _TESTM
#define _TESTM
#if defined(_TESTM)
      eps = 0.01_dp*SQRT(EPSILON(eps))

      ALLOCATE(temps(ntheta,nzeta,ns),                                  &
               tempu(ntheta,nzeta,ns),                                  &
               tempv(ntheta,nzeta,ns), stat=istat)

!     gijf * hj1   (first column of hij matrix)
      CALL tolowerf(hss, hsu, hsv, temps, tempu, tempv, 1, ns)
      r1 = MAXVAL(ABS(temps-1))
      s1 = MAXVAL(ABS(tempu))
      t1 = MAXVAL(ABS(tempv))
      maxtest = MAX(r1,s1,t1)
      CALL ASSERT(maxtest.LE.eps,'ERROR1 IN FULL_MESH_METRICS')

!     gijf * hj2   (second column of hij matrix)
      CALL tolowerf(hsu, huu, huv, temps, tempu, tempv, 1, ns)
      r1 = MAXVAL(ABS(temps))
      s1 = MAXVAL(ABS(tempu-1))
      t1 = MAXVAL(ABS(tempv))
      maxtest = MAX(r1,s1,t1)
      CALL ASSERT(maxtest.LE.eps,'ERROR2 IN FULL_MESH_METRICS')

!     gijf * hj3   (third column of hij matrix)
      CALL tolowerf(hsv, huv, hvv, temps, tempu, tempv, 1, ns)
      r1 = MAXVAL(ABS(temps))
      s1 = MAXVAL(ABS(tempu))
      t1 = MAXVAL(ABS(tempv-1))
      maxtest = MAX(r1,s1,t1)
      CALL ASSERT(maxtest.LE.eps,'ERROR3 IN FULL_MESH_METRICS')
      
      DEALLOCATE(temps, tempu, tempv)
#endif
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


      SUBROUTINE tolowerh(xsupsij, xsupuij, xsupvij,                    &
                          xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: nuv=>nuv_i
!
!     Moves contravariant (upper) to covariant (lower) magnetic field
!     components on the half radial mesh.
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!----------------------------------------------
      INTEGER, INTENT(IN) :: nsmin, nsmax  
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(IN)  ::              &
                   xsupsij, xsupuij, xsupvij 
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(OUT) ::              &
                   xsubsij, xsubuij, xsubvij
!-----------------------------------------------
      xsubsij = gss(:,nsmin:nsmax)*xsupsij                              &
              + gsu(:,nsmin:nsmax)*xsupuij                              &
              + gsv(:,nsmin:nsmax)*xsupvij

      xsubuij = gsu(:,nsmin:nsmax)*xsupsij                              &
              + guu(:,nsmin:nsmax)*xsupuij                              &
              + guv(:,nsmin:nsmax)*xsupvij
      
      xsubvij = gsv(:,nsmin:nsmax)*xsupsij                              &
              + guv(:,nsmin:nsmax)*xsupuij                              &
              + gvv(:,nsmin:nsmax)*xsupvij

      END SUBROUTINE tolowerh
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE tolowerf(xsupsij, xsupuij, xsupvij,                    &
                          xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: ns=>ns_i, nuv=>nuv_i
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
      REAL(dp), DIMENSION(nuv,ns), INTENT(IN)  ::                       &
                xsupsij, xsupuij, xsupvij 
      REAL(dp), DIMENSION(nuv,ns), INTENT(OUT) ::                       &
                xsubsij, xsubuij, xsubvij
!-----------------------------------------------
      INTEGER :: js, js1
!-----------------------------------------------
      DO js=nsmin,nsmax
      js1 = js-nsmin+1
      xsubsij(:,js1) = gssf(:,js)*xsupsij(:,js1)                        &
                     + gsuf(:,js)*xsupuij(:,js1)                        &
                     + gsvf(:,js)*xsupvij(:,js1) 

      xsubuij(:,js1) = gsuf(:,js)*xsupsij(:,js1)                        &
                     + guuf(:,js)*xsupuij(:,js1)                        &
                     + guvf(:,js)*xsupvij(:,js1) 

      xsubvij(:,js1) = gsvf(:,js)*xsupsij(:,js1)                        &
                     + guvf(:,js)*xsupuij(:,js1)                        &
                     + gvvf(:,js)*xsupvij(:,js1)       
      END DO

      END SUBROUTINE tolowerf


      SUBROUTINE cleanup_metric_elements
      USE fourier, ONLY: dealloc_fourier
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER         :: istat
!-----------------------------------------------
	
!     Note: sqrtg is deallocated in init_bcovar and stored in jacob variable

      DEALLOCATE(gss, gsu, gsv, guu, guv, gvv,                          &
                 hss, hsu, hsv, huu, huv, hvv,                          &
                 phipf_i, chipf_i, presf_i, stat=istat)

      CALL dealloc_fourier
      DEALLOCATE(rmnc_i, zmns_i)     ! Used by output subroutines.        
      IF (lasym) DEALLOCATE(rmns_i, zmnc_i)
	
      END SUBROUTINE cleanup_metric_elements

       
      SUBROUTINE dealloc_full_lower_metrics
      USE stel_kinds
      INTEGER :: istat
       
      DEALLOCATE(gssf, guuf, gvvf, gsuf, gsvf, guvf, stat = istat)
      CALL ASSERT(istat.EQ.0, 'Problem dealloc. in DEALLOC_METRIC')
        
      END SUBROUTINE dealloc_full_lower_metrics

      END MODULE metrics
