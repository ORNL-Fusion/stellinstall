MODULE grid_extension
  USE vmec_info, ONLY: jcurrumnc, jcurrvmnc
  USE stel_kinds
  USE stel_constants
  USE shared_data, ONLY: lasym, r1_i, z1_i, ru_i, zu_i, rv_i, zv_i, nsp
  USE island_params, ns=>ns_i, mpol=>mpol_i, ntor=>ntor_i, &
    ntheta=>nu_i, nzeta=>nv_i, mnmax=>mnmax_i, nsv=>ns_v, ohs=>ohs_i
  USE island_params, ONLY: phipf_i, chipf_i, presf_i, nsh
  USE fourier, ONLY: f_none, f_cos, f_sin, f_sum, f_du, f_dv,           &
    cosmui, vtoijsp, orthonorm
  USE v3_utilities, ONLY: assert
  USE siesta_namelist, ONLY: l_vessel
  USE descriptor_mod, ONLY: iam

  IMPLICIT NONE


  ! Vessel information
  INTEGER, PARAMETER :: LIN=1, QUAD=2, CUB=3
  INTEGER :: interp_type
  INTEGER :: nse

  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ri, zi

  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1_vmec, z1_vmec
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ru_vmec, zu_vmec
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rv_vmec, zv_vmec
  
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1s_vmec, r1ss_vmec
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: z1s_vmec, z1ss_vmec
  
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1su_vmec, z1su_vmec
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1sv_vmec, z1sv_vmec

  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rmnc_v, zmns_v
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rmns_v, zmnc_v

  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: r_ext, z_ext
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: ru_ext, zu_ext
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rv_ext, zv_ext

  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1_ves, z1_ves
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ru_ves, zu_ves
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rv_ves, zv_ves

  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: r_full, z_full
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: ru_full, zu_full
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: rv_full, zv_full

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE  :: sqrtg_full

  REAL(dp), PARAMETER :: p5 = 0.5_dp
  REAL(dp)            :: s_edge, ming, maxg

!  LOGICAL :: l_vessel

CONTAINS

  SUBROUTINE grid_extender(wout_file, itype, istat)
    USE fourier, ONLY: toijsp, f_cos, f_sin,                            &
      f_none, f_sum, f_du, f_dv, tomnsp_par
    USE v3_utilities, ONLY: assert
    USE vmec_info, ONLY: rmnc_i, rmns_i, zmns_i, zmnc_i,                &
                         jcurrumnc, jcurrvmnc, jcurrumns, jcurrvmns
    USE descriptor_mod, ONLY: SIESTA_COMM
    USE siesta_namelist, ONLY: nsin_ext

    IMPLICIT NONE

    ! Dummy variables
    CHARACTER(*), INTENT(IN) :: wout_file,  itype
    INTEGER, INTENT(OUT) :: istat

    CHARACTER(LEN=10)   :: chfit
    CHARACTER(LEN=100)  :: device

    INTEGER :: j, m, n, u, v, i, k, js, js1, nsp1
    REAL(dp) :: rho, vol0, vol1, vols, rho1, rho2, rho3, vp1, s_ext, delv
    REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: r12, rs12, zs12, ru12, zu12
    REAL(dp), DIMENSION(:), ALLOCATABLE  :: vp 
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE  :: tmp, tmp2
    REAL(dp), DIMENSION(:), ALLOCATABLE  :: tmp1 

    IF (.NOT.l_vessel) RETURN

    CALL init_extender

    i = INDEX(wout_file, '.', BACK=.TRUE.)
    j = INDEX(wout_file,'_')
    IF (j .LT. i) THEN
      device = wout_file(j+1:i-1)
    ELSE
      device = wout_file(1:i)
    END IF

    !chfit = 'quad'
    !chfit = TRIM(itype)
    IF (INDEX(itype,'lin').EQ.1 .OR. INDEX(itype,'LIN').EQ.1) THEN
      interp_type = LIN; chfit = "linear"
    ELSE IF (INDEX(itype,'quad').EQ.1 .OR. INDEX(itype,'QUAD').EQ.1) THEN
      interp_type = QUAD; chfit = "quadratic"
    ELSE IF (INDEX(itype,'cub').EQ.1 .OR. INDEX(itype,'CUB').EQ.1) THEN
      interp_type = CUB; chfit = "cubic"
    ELSE
      CALL ASSERT(.FALSE., 'Unknown interpolation <TYPE> specified')
    END IF

    ALLOCATE(ri(ntheta,nzeta),zi(ntheta,nzeta))
    ALLOCATE(r1_ves(ntheta,nzeta),z1_ves(ntheta,nzeta), &
      ru_ves(ntheta,nzeta), zu_ves(ntheta,nzeta), &
      rv_ves(ntheta,nzeta), zv_ves(ntheta,nzeta))

    !Set internal normalization
    rmnc_v = rmnc_v/orthonorm(0:mpol_v,-ntor_v:ntor_v)
    zmns_v = zmns_v/orthonorm(0:mpol_v,-ntor_v:ntor_v)

    CALL vtoijsp (rmnc_v, r1_ves, f_none, f_cos)
    CALL vtoijsp (zmns_v, z1_ves, f_none, f_sin)

    CALL vtoijsp (rmnc_v, ru_ves, f_du, f_cos)
    CALL vtoijsp (zmns_v, zu_ves, f_du, f_sin)

    CALL vtoijsp (rmnc_v, rv_ves, f_dv, f_cos)
    CALL vtoijsp (zmns_v, zv_ves, f_dv, f_sin)

    !Compute extended "s" coordinate for s>1, assuming Vp(s) = Vp(s=1)*s
    !where Vp(s=1) is obtained from extrapolating VMEC values. Results is
    !      V(s) = Vp*(s^2-1)/2 + V0 => s_edge = SQRT(2(V1-V0)/Vp + 1)
    !
    vol0 = Get_Volume(zero)    !rho = 0: s=1 vmec
    vol1 = Get_Volume(one)     !rho = 1: vessel

    CALL GetVP1(vp1)
    s_edge = SQRT(1+2*(vol1-vol0)/vp1)

    !NOTE: at js=nse, s(js) == (js-1)*hs = s_edge
    nse = NINT(s_edge*ohs) + 1
!    print *,"SKS> ns: ", ns, " nse: ",nse, "s_edge: ",s_edge, "vp1: ", vp1

    nsv = nse-ns
    s_edge = (nsv-1)*hs_i + 1

    !rounding discrepancy in volume due to integer nse 
    delv = vol1 - (vol0 + p5*vp1*(s_edge**2-1))
    delv = delv/(s_edge-1)**2

!    print '(a,3i4,1p,3E12.4)'," --- 1]: ", nse, nsv, ns,               &
!      s_edge, vp1/(2*vol0), delv/vol1

    ALLOCATE(r_ext (ntheta,nzeta,2:nsv), z_ext (ntheta,nzeta,2:nsv),    &
             ru_ext(ntheta,nzeta,2:nsv), zu_ext(ntheta,nzeta,2:nsv),    &
             rv_ext(ntheta,nzeta,2:nsv), zv_ext(ntheta,nzeta,2:nsv),    &
             stat=istat)

    DO j = 2, nsv

      s_ext = 1 + (j-1)*hs_i    !at j=nsv, s_ext = s_edge

      !V(s) = Vp*(s^2-1)/2 + V0
      vols = vol0 + p5*vp1*(s_ext**2-1)   !matches vp exactly at s=1
      vols = vols + delv*(s_ext-1)**2     !matches vol1 exactly at s=s_edge
      vols = MIN(vols, vol1)
      CALL remap_radial(rho, vols)

      rho1 = 1 - rho
      IF (interp_type.EQ.LIN) THEN
        r_ext(:,:,j) = rho1*r1_vmec + rho*r1_ves
        z_ext(:,:,j) = rho1*z1_vmec + rho*z1_ves
        ru_ext(:,:,j) = rho1*ru_vmec + rho*ru_ves 
        zu_ext(:,:,j) = rho1*zu_vmec + rho*zu_ves 
        rv_ext(:,:,j) = rho1*rv_vmec + rho*rv_ves 
        zv_ext(:,:,j) = rho1*zv_vmec + rho*zv_ves 
      ELSE
        rho2 = rho*rho
        IF (interp_type .EQ. QUAD) THEN
          !QUADRATIC (R' match at rho = 0)
          r_ext(:,:,j) = p5*r1s_vmec*(1 - rho2 - rho1*rho1) + &
            r1_vmec*(1-rho2) + r1_ves*rho2
          z_ext(:,:,j) = p5*z1s_vmec*(1 - rho2 - rho1*rho1) + &
            z1_vmec*(1-rho2) + z1_ves*rho2
          ru_ext(:,:,j) = p5*r1su_vmec*(1 - rho2 - rho1*rho1) + &
            ru_vmec*(1-rho2) + ru_ves*rho2
          zu_ext(:,:,j) = p5*z1su_vmec*(1 - rho2 - rho1*rho1) + &
            zu_vmec*(1-rho2) + zu_ves*rho2
          rv_ext(:,:,j) = p5*r1sv_vmec*(1 - rho2 - rho1*rho1) + &
            rv_vmec*(1-rho2) + rv_ves*rho2
          zv_ext(:,:,j) = p5*z1sv_vmec*(1 - rho2 - rho1*rho1) + &
            zv_vmec*(1-rho2) + zv_ves*rho2
        ELSE IF (interp_type .EQ. CUB) THEN
          rho3 = rho*rho2
          CALL ASSERT(.FALSE., 'Unknown interpolation <TYPE> specified')
        END IF
      END IF
    END DO

    nse = nse-1 !Not to double count s=1 surface
    ALLOCATE(r_full(ntheta,nzeta,nse), z_full(ntheta,nzeta,nse),        &
      ru_full(ntheta,nzeta,nse), zu_full(ntheta,nzeta,nse), &
      rv_full(ntheta,nzeta,nse), zv_full(ntheta,nzeta,nse))
 
    r_full (:,:,1:ns) = r1_i; r_full (:,:,ns+1:nse) = r_ext
    z_full (:,:,1:ns) = z1_i; z_full (:,:,ns+1:nse) = z_ext
    ru_full(:,:,1:ns) = ru_i; ru_full(:,:,ns+1:nse) = ru_ext
    zu_full(:,:,1:ns) = zu_i; zu_full(:,:,ns+1:nse) = zu_ext
    rv_full(:,:,1:ns) = rv_i; rv_full(:,:,ns+1:nse) = rv_ext
    zv_full(:,:,1:ns) = zv_i; zv_full(:,:,ns+1:nse) = zv_ext

!    CALL PLOTRZ (device, chfit)
    nse = MIN(ns + nsin_ext, nse)
    WRITE (*,1000) nse
1000 FORMAT('\Using ',i4,' Surfaces\')

    !REASSIGN THE NEW VALUE OF NS
    nsp = ns
    ns = nse
    nsh = ns-1

    IF (nsp .NE. ns) THEN
       DEALLOCATE(r1_i);
       DEALLOCATE(z1_i);
       DEALLOCATE(ru_i);
       DEALLOCATE(zu_i);
       DEALLOCATE(rv_i);
       DEALLOCATE(zv_i);
       ALLOCATE(r1_i(ntheta,nzeta,ns), z1_i(ntheta,nzeta,ns),           &
                ru_i(ntheta,nzeta,ns), zu_i(ntheta,nzeta,ns),           &
                rv_i(ntheta,nzeta,ns), zv_i(ntheta,nzeta,ns))
       r1_i=r_full(:,:,1:ns)
       z1_i=z_full(:,:,1:ns)
       ru_i=ru_full(:,:,1:ns)
       zu_i=zu_full(:,:,1:ns)
       rv_i=rv_full(:,:,1:ns)
       zv_i=zv_full(:,:,1:ns)
    END IF

    DEALLOCATE(r_full,z_full)
    DEALLOCATE(ru_full,zu_full)
    DEALLOCATE(rv_full,zv_full)
    
    IF (nsp .NE. ns) THEN
       nsp1 = nsp+1
       ALLOCATE(tmp(0:mpol,-ntor:ntor,ns), tmp2(0:mpol,-ntor:ntor, ns))
       tmp(:,:,1:nsp) = rmnc_i; tmp2(:,:,1:nsp) = zmns_i
       DEALLOCATE(rmnc_i, zmns_i)
       ALLOCATE(rmnc_i(0:mpol,-ntor:ntor,ns),                           &
                zmns_i(0:mpol,-ntor:ntor,ns), stat=istat)
       rmnc_i(:,:,1:nsp) = tmp(:,:,1:nsp); zmns_i(:,:,1:nsp) = tmp2(:,:,1:nsp)
       CALL tomnsp_par(r1_i(:,:,nsp1:ns), rmnc_i(:,:,nsp1:ns), f_cos)
       CALL tomnsp_par(z1_i(:,:,nsp1:ns), zmns_i(:,:,nsp1:ns), f_sin)
       DO js = nsp1, ns
          rmnc_i(:,:,js) = rmnc_i(:,:,js)/orthonorm(0:mpol,-ntor:ntor)
          zmns_i(:,:,js) = zmns_i(:,:,js)/orthonorm(0:mpol,-ntor:ntor)
       END DO

! PUT PROFILES ONTO EXTENDED MESH
       tmp=0
       tmp(:,:,1:nsp)=jcurrumnc
       DEALLOCATE(jcurrumnc)
       ALLOCATE(jcurrumnc(0:mpol,-ntor:ntor,ns))
       jcurrumnc=tmp

       tmp=0
       tmp(:,:,1:nsp)=jcurrvmnc
       DEALLOCATE(jcurrvmnc)
       ALLOCATE(jcurrvmnc(0:mpol,-ntor:ntor,ns))
       jcurrvmnc=tmp

! DO THE SAME FOR THE phipf_i, chipf_i AND presf_i
! phipf, chipf will be computed in PCHELMS in the vacuum region
       ALLOCATE(tmp1(ns))
       tmp1=0
       tmp1(1:nsp)=phipf_i
       DEALLOCATE(phipf_i)
       ALLOCATE(phipf_i(ns))
       phipf_i=tmp1

       tmp1=0
       tmp1(1:nsp)=chipf_i
       DEALLOCATE(chipf_i)
       ALLOCATE(chipf_i(ns))
       chipf_i=tmp1

       tmp1(1:nsp)=presf_i
       tmp1(nsp+1:) = presf_i(nsp)
       DEALLOCATE(presf_i)
       ALLOCATE(presf_i(ns))
       presf_i=tmp1
       DEALLOCATE(tmp1)

       IF (lasym) THEN
          tmp(:,:,1:nsp) = rmns_i; tmp2(:,:,1:nsp) = zmnc_i
          DEALLOCATE(rmns_i, zmnc_i)
          ALLOCATE(rmns_i(0:mpol,-ntor:ntor,ns),                        &
                   zmnc_i(0:mpol,-ntor:ntor,ns), stat=istat)
          rmns_i(:,:,1:nsp) = tmp(:,:,1:nsp); zmnc_i(:,:,1:nsp) = tmp2(:,:,1:nsp)
          CALL tomnsp_par(r1_i(:,:,nsp1:ns), rmns_i(:,:,nsp1:ns), f_sin)
          CALL tomnsp_par(z1_i(:,:,nsp1:ns), zmnc_i(:,:,nsp1:ns), f_cos)
          DO js = nsp1, ns
             rmns_i(:,:,js) = rmns_i(:,:,js)/orthonorm(0:mpol,-ntor:ntor)
             zmnc_i(:,:,js) = zmnc_i(:,:,js)/orthonorm(0:mpol,-ntor:ntor)
          END DO

          tmp=0
          tmp(:,:,1:nsp)=jcurrumns
          DEALLOCATE(jcurrumns)
          ALLOCATE(jcurrumns(0:mpol,-ntor:ntor,ns))
          jcurrumns=tmp

          tmp=0
          tmp(:,:,1:nsp)=jcurrvmns
          DEALLOCATE(jcurrvmns)
          ALLOCATE(jcurrvmns(0:mpol,-ntor:ntor,ns))
          jcurrvmns=tmp
       END IF

       DEALLOCATE(tmp, tmp2)

    END IF

#if defined(NOSKS)
    DO u=1, ntheta
      DO v=1, nzeta
        DO js=1, ns
          WRITE(1600,*) r1_i(u,v,js), z1_i(u,v,js)
          CALL FLUSH(1600)
        END DO
      END DO
    END DO

    ! u = theta, v = zeta
    ALLOCATE (r12 (ntheta,nzeta), &
      rs12(ntheta,nzeta), &
      zs12(ntheta,nzeta), &
      ru12(ntheta,nzeta), &
      zu12(ntheta,nzeta), &
      sqrtg_full(ntheta,nzeta,nse), &
      stat = istat)
    CALL ASSERT(istat.eq.0,' Allocation error in grid_extender')
    sqrtg_full(:,:,1) = 0
    DO js = 2, nse
      js1 = js-1
      r12  = p5*(r1_i(:,:,js) + r1_i(:,:,js1))
      rs12 = (r1_i(:,:,js) - r1_i(:,:,js1))*ohs
      zs12 = (z1_i(:,:,js) - z1_I(:,:,js1))*ohs
      ru12 = p5*(ru_i(:,:,js) + ru_i(:,:,js1))
      zu12 = p5*(zu_i(:,:,js) + zu_i(:,:,js1))
      sqrtg_full(:,:,js) = r12*(ru12*zs12 - rs12*zu12)
    END DO

    DEALLOCATE (r12, rs12, zs12, ru12, zu12)

    ming = MINVAL(sqrtg_full(:,:,ns:nse) ) 
    maxg = MAXVAL(sqrtg_full(:,:,ns:nse))
    PRINT *,' IAM: ', iam, 'MING: ', ming,' MAXG: ',maxg
    CALL ASSERT(ming*maxg.GT.zero,' Jacobian changed sign in grid_extender')

    ALLOCATE (vp(nse))
    CALL SurfAvgLocal(vp,sqrtg_full,1,nse)
    vp(1) = 0
    DO js=2, nse
      ming = MINVAL(sqrtg_full(:,:,js) ) 
      maxg = MAXVAL(sqrtg_full(:,:,js))
      WRITE(5000,'(i4,1p,3E14.4)') js, -vp(js), ming, maxg
    END DO
    CALL FLUSH(5000)
    DEALLOCATE (vp,r_ext, z_ext, ru_ext, zu_ext)
    DEALLOCATE (ri,zi)
#endif

    DEALLOCATE(r1s_vmec, r1su_vmec, r1sv_vmec, z1s_vmec, z1su_vmec, z1sv_vmec)

  END SUBROUTINE grid_extender

  SUBROUTINE init_extender
    REAL(dp) :: fac

    ALLOCATE(r1_vmec(ntheta,nzeta),z1_vmec(ntheta,nzeta))
    ALLOCATE(ru_vmec(ntheta,nzeta),zu_vmec(ntheta,nzeta))
    ALLOCATE(rv_vmec(ntheta,nzeta),zv_vmec(ntheta,nzeta))
    r1_vmec = r1_i(:,:,ns); z1_vmec = z1_i(:,:,ns)
    ru_vmec = ru_i(:,:,ns); zu_vmec(:,:) = zu_i(:,:,ns)
    rv_vmec = rv_i(:,:,ns); zv_vmec(:,:) = zv_i(:,:,ns)

    !First derivatives on s=1 surface
    ALLOCATE(r1s_vmec(ntheta,nzeta),z1s_vmec(ntheta,nzeta))
    ALLOCATE(r1su_vmec(ntheta,nzeta),z1su_vmec(ntheta,nzeta))
    ALLOCATE(r1sv_vmec(ntheta,nzeta),z1sv_vmec(ntheta,nzeta))

    fac =.5  !COULD SCAN THIS TO MAKE MAX(sqrt(g) < 0) IF NECESSARY
!!    fac =.17  !COULD SCAN THIS TO MAKE MAX(sqrt(g) < 0) IF NECESSARY
    r1s_vmec=(r1_vmec - r1_i(:,:,ns-1))*ohs * fac
    z1s_vmec=(z1_vmec - z1_i(:,:,ns-1))*ohs * fac

    r1su_vmec = (ru_vmec - ru_i(:,:,ns-1))*ohs * fac
    z1su_vmec = (zu_vmec - zu_i(:,:,ns-1))*ohs * fac

    r1sv_vmec = (rv_vmec - rv_i(:,:,ns-1))*ohs * fac
    z1sv_vmec = (zv_vmec - zv_i(:,:,ns-1))*ohs * fac

  END SUBROUTINE init_extender


  SUBROUTINE GetVP1(vp1)
    REAL(dp), INTENT(OUT) :: vp1
    INTEGER :: js, js1
    REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: r12, rs12, zs12, ru12, zu12
    REAL(dp), DIMENSION(:), ALLOCATABLE    :: vp 

    ALLOCATE(vp(ns-1:ns), sqrtg_full(ntheta,nzeta,ns-1:ns))
    ALLOCATE(r12(ntheta,nzeta), ru12(ntheta,nzeta), rs12(ntheta,nzeta), &
      zu12(ntheta,nzeta), zs12(ntheta,nzeta))

    !Guess assuming V(s) = V0*s**2:    vp1 = 2*vol0
    !NOTE: jacobian not computed yet

    DO js = ns-1, ns
      js1 = js-1
      r12 = p5*(r1_i(:,:,js) + r1_i(:,:,js1))
      rs12 = (r1_i(:,:,js) - r1_i(:,:,js1))*ohs
      zs12 = (z1_i(:,:,js) - z1_i(:,:,js1))*ohs
      ru12= p5*(ru_i(:,:,js) + ru_i(:,:,js1))
      zu12= p5*(zu_i(:,:,js) + zu_i(:,:,js1))
      sqrtg_full(:,:,js) = r12*(ru12*zs12 - rs12*zu12)
    END DO

    CALL SurfAvgLocal(vp,sqrtg_full,ns-1,ns)   
    vp = twopi**2*ABS(vp)
    vp1 = p5*(3*vp(ns) - vp(ns-1))

    DEALLOCATE(vp, sqrtg_full, r12, ru12, rs12, zu12, zs12)

  END SUBROUTINE GetVP1

#if 0
  SUBROUTINE PlotRZ (device, chfit)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: device, chfit

    CHARACTER*200 :: plabel
    INTEGER :: PGOPEN, s, s0, s1, soff, u, v, vr, k, npts, nc, nc0
    REAL :: rmin, rmax, zmin, zmax, rho, xc
    REAL, ALLOCATABLE :: X1(:), Y1(:)
    REAL, ALLOCATABLE :: X2(:,:), Y2(:,:)

    rmin = MINVAL(r_full); rmax = MAXVAL(r_full)
    zmin = MINVAL(z_full); zmax = MAXVAL(z_full)
    IF (.NOT.lasym) THEN
      zmax = MAX(zmax, -zmin)
      zmin = MIN(-zmax, zmin)
    ENDIF
    !Expand viewport a bit
    rho = 1.1*MAX((rmax-rmin), (zmax-zmin))/2
    xc = (rmax+rmin)/2
    rmin = xc - rho;  rmax = xc + rho
    xc = (zmax+zmin)/2
    zmin = xc - rho;  zmax = xc + rho

    npts = ntheta+1
    IF (.NOT.lasym) npts = 2*ntheta

    nc = 0
    DO s = 1, nse, ns/10
      nc = nc+1
    END DO
    nc0 = nc

    ALLOCATE (x1(npts), y1(npts), x2(nc,npts), y2(nc,npts))

    IF (PGOPEN('?') .LT. 1) STOP

    ! Reverse background/foreground colors
    CALL PGSCR(1,0.0,0.0,0.0)
    CALL PGSCR(0,1.0,1.0,1.0)

    PHILOOP: DO v = 1, NZETA

      ! Start new page
      CALL PGPAGE
      CALL PGPAP(8.0,1.0)   !Sets Aspect ratio of axis to 1 on 8" page

      ! Erase all graphics from the current page and fill with background color
      CALL PGERAS
      !
      ! Set up window and viewport.
      !
      CALL PGSCH(1.0)
      CALL PGVSTD
      CALL PGSWIN(rmin,rmax,zmin,zmax)
      !
      ! Fill the viewport in a different color.
      !
      CALL PGSCI(0)
      CALL PGRECT(rmin,rmax,zmin,zmax)

      ! Draw thin grid lines
      ! SET COLOR INDEX (SCI) AND LINE WIDTH (SLW)
      CALL PGSCI(1)
      CALL PGSLW(1)
      CALL PGBOX('G', 0.0, 0, 'G', 0.0, 0)

      ! Draw thick box, axis line
      CALL PGSCI(1)
      CALL PGSLW(3)
      CALL PGSCF(2)
      CALL PGBOX('BCNST', 0.0, 0, 'BCNSTV', 0.0, 0)

      K = NINT(REAL(360*(V-1))/NZETA)
      WRITE (plabel,'(3a,i3,a)') '\fr Device: ',TRIM(device),  &
        ' \fi \gz = ',K,'\u\(0015)'
      !	         TRIM(chfit),' fit: \fi \(0647) = ',K,'\u\(0015)'
      CALL PGLAB('\fiR', '\fiZ', TRIM(plabel))
      WRITE (plabel,'(2a)') '\fr Extended Coordinate: ', &
        TRIM(chfit)
      CALL PGMTXT('T',1.0,0.23,0.0,plabel)

      nc = 0
      soff = nse-1 - (nc0-1)*(ns/10)
      s1 = 0
      SURF_LOOP: DO s = 1, nse, ns/10

        CALL PGSCI(4)
        IF (s .GT. ns) CALL PGSCI(3)
        CALL PGSLW(3)

        !Plot the outer-most surface thicker, different color
        s0 = s
        nc = nc+1
        IF (soff .GT. 0) THEN
          soff = soff-1
          s1 = s1 + 1
        END IF
        s0 = s0 + s1
        IF (nc .EQ. nc0) THEN
          s0 = nse
          CALL PGSCI(6)
          CALL PGSLW(4)
        ELSE IF (nc .GT. nc0) THEN
          EXIT
        END IF


        x1(1:ntheta) = r_full(:,v,s0)
        y1(1:ntheta) = z_full(:,v,s0)
        IF (lasym) THEN
          x1(npts) = x1(1)
          y1(npts) = y1(1)
        ELSE
          vr = nzeta - (v-1)          !reflected vr = 2pi-v
          vr = MOD(vr, nzeta)+1
          x1(ntheta+1:npts) =  r_full(ntheta:1:-1,vr,s0)
          y1(ntheta+1:npts) = -z_full(ntheta:1:-1,vr,s0)
        END IF

        x2(nc,:) = x1;  y2(nc,:) = y1

        CALL PGLINE(npts,x1,y1)    

      END DO SURF_LOOP

      !Plot Angle contours
      DO u = 1, npts
        CALL PGSCI(2)
        CALL PGSLW(2)
        CALL PGLINE(NC,x2(:,u),y2(:,u))
      END DO

    END DO PHILOOP

    CALL PGEND

    DEALLOCATE(x1, y1, x2, y2)

  END SUBROUTINE PlotRZ
#endif

  SUBROUTINE remap_radial(rho, vrho)
    IMPLICIT NONE
    INTEGER                 :: it
    REAL(dp), INTENT(INOUT) :: rho
    REAL(dp), INTENT(IN)    :: vrho
    REAL(dp)                :: rootl, rootr, vol, delv
    REAL(dp), PARAMETER     :: tol=1.E-6_dp
    !
    !     ROOT FINDER: FINDS THE NEW rho VALUE CORRESPONDING TO THE 
    !     INITIAL VOLUME VRHO
    !
    rho = p5
    rootr = 1
    rootl = 0

    DO it = 1, 100
      vol = get_volume(rho)
      delv = vol-vrho
      IF (ABS(delv) .LT. tol*vrho) EXIT
      IF (delv .GT. 0) THEN
        rootr = rho
        rho = (rootl+rho)/2
      ELSE
        rootl = rho
        rho = rootl+(rootr-rho)/2
      END IF
    END DO      
  END SUBROUTINE remap_radial

  REAL(dp) FUNCTION get_volume (rho)
    IMPLICIT NONE
    !     
    !     USES FORMULA VOLUME(rho) = INTEGRAL[.5*(R^2 * ZU)] 
    !     WHERE INTEGRATION IS OVER U,V FOR A FIXED rho=X
    !     FOR SIMPLICITY, WE'LL DO THE INTEGRAL PLANE BY PLANE (WHICH IS NOT QUITE RIGHT)
    !
    REAL(dp) :: rho
    INTEGER  :: i, j
    REAL(dp) :: zu, r1, dnorm, odelu

    get_volume = 0
    odelu = ntheta
    odelu = odelu/twopi
    dnorm = p5*twopi/(nzeta*odelu)
    IF (.NOT. lasym) odelu = 2*odelu
    CALL LoadSurfacePoints(rho)
    DO j = 1, nzeta
      DO i = 1, ntheta-1
        zu = (zi(i+1,j)-zi(i,j))*odelu
        r1 = p5*(ri(i+1,j)+ri(i,j))
        get_volume = get_volume + r1**2 * zu
      END DO
      IF (lasym) THEN
        zu = (zi(1,j)-zi(ntheta,j))*odelu
        r1 = p5*(ri(1,j)+ri(ntheta,j))
        get_volume = get_volume + r1**2*zu
      END IF
    END DO
    get_volume = dnorm*get_volume

  END FUNCTION get_volume

  SUBROUTINE LoadSurfacePoints(rho)
    USE island_params, ns=>ns_i, mpol=>mpol_i, ntor=>ntor_i, &
      ntheta=>nu_i, nzeta=>nv_i, mnmax=>mnmax_i, nsv=>ns_v
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: rho
    REAL(dp)             :: rho1, rho2, rho3

    !LOADS THE THETA AND ZETA VALUES OF RI, ZI ARRAYS 
    !FOR THE INPUT VALUE OF rho ("radial" surface)
    rho1 = 1 - rho
    IF (interp_type .EQ. LIN) THEN
      ri = rho1*r1_vmec + rho*r1_ves
      zi = rho1*z1_vmec + rho*z1_ves
    ELSE
      rho2 = rho*rho
      IF (interp_type .EQ. QUAD) THEN
        !QUADRATIC (R' match at rho = 0)
        ri = p5*r1s_vmec*(1 - rho2 - rho1*rho1) &
          + r1_vmec*(1-rho2) + r1_ves*rho2
        zi = p5*z1s_vmec*(1 - rho2 - rho1*rho1) &
          + z1_vmec*(1-rho2) + z1_ves*rho2
      ELSE IF (interp_type .EQ. CUB) THEN
        !CUBIC     (R', R'' match at rho=0)
        rho3 = rho*rho2
        CALL ASSERT(.FALSE., 'Unknown interpolation <TYPE> specified')
      END IF
    END IF
  END SUBROUTINE LoadSurfacePoints

  SUBROUTINE read_vessel_file(wout_file, istat)
    USE island_params, ns=>ns_i, mpol=>mpol_i, ntor=>ntor_i, &
      ntheta=>nu_i, nzeta=>nv_i, mnmax=>mnmax_i
    USE siesta_namelist, ONLY: vessel_file
    IMPLICIT NONE

    ! Dummy variables
    CHARACTER(*), INTENT(IN) :: wout_file
    INTEGER, INTENT(OUT) :: istat

    ! Local variables
    CHARACTER(LEN=400) :: cprefix
    INTEGER :: j, l1, m, n, mn, mnmax_v
    LOGICAL :: lex
    REAL(dp) :: rv1, zv1

    mpol_v = 0; ntor_v =0

    l1=LEN(TRIM(wout_file))
    WRITE(cprefix,*) wout_file
    cprefix=TRIM(ADJUSTL(cprefix))

    INQUIRE(file=vessel_file, exist=lex)
    IF (.not.lex) THEN
!        PRINT *, 'Message from line ', __LINE__,' of ', __FILE__
      STOP "Vessel file not found."
    END IF

    IF (lex) THEN
      istat = 0
      mnmax_v = 0
      OPEN(UNIT=36, FILE=vessel_file,STATUS="old",FORM="formatted", iostat=j)
      IF (j .NE. 0) RETURN
      DO WHILE(.TRUE.)
        READ(36, *, END=101) n, m, rv1, zv1
        mpol_v=MAX(m,mpol_v)
        ntor_v=MAX(ABS(n),ntor_v)
        mnmax_v = mnmax_v + 1
      END DO
      101   CONTINUE

      REWIND(UNIT=36,iostat=istat)
      IF(istat .NE. 0) RETURN
      ALLOCATE(rmnc_v(0:mpol_v,-ntor_v:ntor_v), zmns_v(0:mpol_v,-ntor_v:ntor_v))
      rmnc_v=0; zmns_v=0 
      IF (lasym) THEN
        ALLOCATE(rmns_v(0:mpol_v,-ntor_v:ntor_v), zmnc_v(0:mpol_v,-ntor_v:ntor_v))
        rmns_v=0; zmnc_v=0 
      END IF

      DO mn = 1, mnmax_v
        READ(36, *, END=102) n, m, rv1, zv1
        CALL ASSERT(m.GE.0.AND.m.LE.mpol_v,'m out of bounds: read_vessel_file')
        CALL ASSERT(n.GE.-ntor_v.AND.n.LE.ntor_v, 'n out of bounds: read_vessel_file')
        rmnc_v(m,-n) = rv1;  zmns_v(m,-n) = zv1
      END DO
      102   CONTINUE 
      CLOSE(UNIT=36)
      l_vessel = .TRUE.

    ELSE IF (.NOT.lex) THEN
      l_vessel = .FALSE.
      istat = 1
    END IF

  END SUBROUTINE read_vessel_file

  SUBROUTINE SurfAvgLocal(average, q3d, nsmin, nsmax)
    USE island_params, ONLY: cosmui 
    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: nsmin, nsmax
    REAL(dp), INTENT(OUT) :: average(nsmin:nsmax)
    REAL(dp), INTENT(IN)  :: q3d(ntheta,nzeta,nsmin:nsmax)

    INTEGER :: l, js
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: wint

    ALLOCATE(wint(ntheta,nzeta,nsmin:nsmax))

    DO l = 1, ntheta
      wint(l,:,:) = cosmui(l,0)
    END DO         
    DO js = nsmin, nsmax
      average(js) = SUM(q3d(:,:,js)*wint(:,:,js))
    END DO

    DEALLOCATE(wint)
  END SUBROUTINE SurfAvgLocal

END MODULE grid_extension
