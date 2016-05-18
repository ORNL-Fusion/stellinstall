      SUBROUTINE spline_fields(arg1)
      USE stel_kinds
      USE stel_constants
      USE bspline
      USE spline3d_fit_coefs
      USE mpi_params
#if defined(MPI_OPT)
      USE mpi_inc
#endif
      IMPLICIT NONE

      CHARACTER(LEN=60), INTENT(IN):: arg1
      INTEGER:: i, j, k, istat, iphi0, itht0, is0, i1, j1,
     &   k1, isize
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
      REAL(rprec):: maxabs, locabs, diff, c1
      REAL(rprec), ALLOCATABLE, DIMENSION(:):: s, u, v
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: r, z, phi,
     &  br, bz, bphi
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: st, ut, vt
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: r_test,
     &  z_test, phi_test, br_test, bz_test, bphi_test

      OPEN(unit=4, file=TRIM(arg1), status="old", iostat=istat)            ! File produced by SIESTA containing data
	IF (istat .ne. 0) STOP 'Cannot open specified file!'
      READ(4,*) iphi0, is, itht0                                           ! IPHI0 = Number of toroidal nodes in file; IS = number o
      iphi = iphi0 + 1; itht = itht0 + 1; is0 = is -1                      ! Include U = 2*pi, V = 2*pi when splining. Thus, ITHT =
!
      CALL alloc_splines                                                   ! Allocate quantities for splines which will be passed vi
      ALLOCATE(s(is), u(itht), v(iphi), STAT=istat)                        ! Allocate local quantities
      IF (istat .NE. 0) STOP 'Allocation error'
      ALLOCATE(r(is,itht,iphi), z(is,itht,iphi),
     &  phi(is,itht,iphi), br(is,itht,iphi),
     &  bz(is,itht,iphi), bphi(is,itht,iphi), STAT = istat)
      IF (istat .NE. 0) STOP 'Allocation error'
      s = zero; v = zero; u = zero
      r = zero; z = zero; phi = zero
      br = zero; bz = zero; bphi = zero
!
      DO k = 1, iphi                                                       ! Construct meshes. Assumed uniform.
        v(k) = twopi*REAL(k-1,dp)/REAL(iphi0,dp)                           ! Extend last points in u, v to include 2pi:
      ENDDO
      DO i = 1, is
        s(i) = REAL(i-1,dp)/REAL(is0,dp)
      ENDDO
      DO j = 1, itht
        u(j) = twopi*REAL(j-1,dp)/REAL(itht0,dp)
      ENDDO
      if(myid .eq. 0) then
      WRITE(*,*)
      WRITE(*,*)' 3D fields being stored in file: ',
     &   TRIM(ADJUSTL(arg1))
      WRITE(*,'(a, i6, a8, i6, a8, i6)') '  Dimensions:  Ns = ', is,
     &   '   Nu = ', itht0, '   Nv = ', iphi0
      WRITE(*,*)' Reading ',TRIM(arg1),'..........'
      endif      !if(myid .eq. 0)
!OLD ORDER      DO k = 1, iphi0      !EARLIER THAN 4-16-2013 - SPH
!OLD ORDER        DO i = 1, is
!OLD ORDER          DO j = 1, itht0
      DO j = 1, itht0
        DO k = 1, iphi0
          DO i = 1, is
            READ(4,*) r(i,j,k), z(i,j,k), phi(i,j,k),
     &        br(i,j,k), bz(i,j,k), bphi(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(4)

      if(myid .eq. 0) then
      WRITE(*,'("  .........completed")')
      WRITE(*,*)
      endif   !if(myid .eq. 0) then

      CALL enforce_periodicity(is, itht, iphi, r, 0)                       ! Enforce periodicity of fields in U and V
      CALL enforce_periodicity(is, itht, iphi, z, 0)
      CALL enforce_periodicity(is, itht, iphi, phi, 1)                     ! Use 1 cause PHI is not periodic in V
      CALL enforce_periodicity(is, itht, iphi, br, 0)
      CALL enforce_periodicity(is, itht, iphi, bz, 0)
      CALL enforce_periodicity(is, itht, iphi, bphi, 0)

      if(myid .eq. 0) then
      WRITE(*,*) ' Creating spline knots'
      endif   !if(myid .eq. 0) then

      CALL dbsnak(is, s, ks, s_knots)
      CALL dbsnak(iphi, v, kv, v_knots)
      CALL dbsnak(itht, u, ku, u_knots)

      if(myid .eq. 0) then
      WRITE(*,*) ' Computing spline coefficients'
      endif   !if(myid .eq. 0) then

      call dbs3in(is, s, itht, u, iphi, v, r, is, itht, ks, ku, kv,
     >  s_knots, u_knots, v_knots, r_coef)
      call dbs3in(is, s, itht, u, iphi, v, z, is, itht, ks, ku, kv,
     >  s_knots, u_knots, v_knots, z_coef)
      call dbs3in(is, s, itht, u, iphi, v, phi, is, itht, ks, ku, kv,
     >  s_knots, u_knots, v_knots, phi_coef)
      call dbs3in(is, s, itht, u, iphi, v, br, is, itht, ks, ku, kv,
     >  s_knots, u_knots, v_knots, br_coef)
      call dbs3in(is, s, itht, u, iphi, v, bz, is, itht, ks, ku, kv,
     >  s_knots, u_knots, v_knots, bz_coef)
      call dbs3in(is, s, itht, u, iphi, v, bphi, is, itht, ks, ku, kv,
     >  s_knots, u_knots, v_knots, bphi_coef)
!
!     Test quality of spline at mid point of ALL 3D cells
!       by measuring the maximum relative difference between averaged
!       and splined value at center of 3D cell for ALL splined fields
!
      IF (l_test_splines) THEN

      if(myid .eq. 0) then
        WRITE(*,*) ' Testing quality of splines'
      endif    !  if(myid .eq. 0) then
        isize = is0*itht0*iphi0
        ALLOCATE(st(is0, itht0, iphi0), ut(is0, itht0, iphi0),
     &    vt(is0, itht0, iphi0),STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        st = zero; ut = zero; vt = zero

        DO i = 1, is0
           st(i,:,:) = (REAL(i,dp)-p5)/REAL(is,dp)
        ENDDO
        DO j = 1, itht0
           ut(:,j,:) = twopi*(REAL(j,dp)-p5)/REAL(itht,dp)
        ENDDO
        DO k = 1, iphi0
           vt(:,:,k) = twopi*(REAL(k,dp)-p5)/REAL(iphi,dp)
        ENDDO

        ALLOCATE(br_test(is0, itht0, iphi0), STAT=istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        br_test = zero
        CALL dbs3vd(0, 0, 0, isize, st, ut, vt, ks, ku, kv,
     &    s_knots, u_knots, v_knots, is, itht, iphi,
     &    br_coef, br_test, is_uniformx, is_uniformy,
     &    is_uniformz)

        DO i = 1, is0
          i1 = i + 1
          DO j = 1, itht0
            j1 = j + 1
            DO k = 1, iphi0
              k1 = k + 1
              c1 = br(i,j,k) + br(i1,j,k)
     &          + br(i,j1,k) + br(i,j,k1)
     &          + br(i,j1,k1) + br(i1,j,k1)
     &          + br(i1,j1,k) + br(i1,j1,k1)
              c1 = c1/8.d0                                                 ! Average value at center of 3D cell
              diff = ABS(c1 - br_test(i,j,k))
              IF (i==1 .AND. j == 1 .AND. k == 1) THEN
                maxabs  = diff
                locabs = c1
              ELSE
                maxabs = MAX(maxabs, diff)                              ! Maximum (absolute) discrepancy
                IF (maxabs == diff) locabs = c1                         ! Store value at which max. discrepancy happens
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      if(myid .eq. 0) then
        WRITE(*, 12) maxabs, locabs
      endif    !  if(myid .eq. 0) then
12      FORMAT('   - Max. absolute difference (Br):', f12.4,
     &    ' at value:', f12.4)

        DEALLOCATE(br_test)

        ALLOCATE(bz_test(is0, itht0, iphi0), STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        bz_test = zero
        CALL dbs3vd(0, 0, 0, isize, st, ut, vt, ks, ku, kv,
     &    s_knots, u_knots, v_knots, is, itht, iphi,
     &    bz_coef, bz_test, is_uniformx, is_uniformy,
     &    is_uniformz)

        DO i = 1, is0
          i1 = i + 1
          DO j = 1, itht0
            j1 = j + 1
            DO k = 1, iphi0
              k1 = k + 1
              c1 = bz(i,j,k) +  bz(i1,j,k)
     &          + bz(i,j1,k) + bz(i,j,k1)
     &          + bz(i,j1,k1) + bz(i1,j,k1)
     &          + bz(i1,j1,k) + bz(i1,j1,k1)
              c1 = c1/8.d0                                                 ! Average value at center of 3D cell
              diff = ABS(c1 - bz_test(i,j,k))
              IF (i==1 .AND. j == 1 .AND. k == 1) THEN
                maxabs  = diff
                locabs = c1
              ELSE
                maxabs = MAX(maxabs, diff)                              ! Maximum (absolute) discrepancy
                IF (maxabs == diff) locabs = c1                         ! Store value at which max. discrepancy happens
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      if(myid .eq. 0) then
        WRITE(*, 13) maxabs, locabs
      endif    !  if(myid .eq. 0) then
13      FORMAT('   - Max. absolute difference (Bz):', f12.4,
     &    ' at value:', f12.4)
        DEALLOCATE(bz_test)

        ALLOCATE(bphi_test(is0, itht0, iphi0), STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        bphi_test = zero
        CALL dbs3vd(0, 0, 0, isize, st, ut, vt, ks, ku, kv,
     &    s_knots, u_knots, v_knots, is, itht, iphi,
     &    bphi_coef, bphi_test, is_uniformx, is_uniformy,
     &    is_uniformz)

        DO i = 1, is0
          i1 = i + 1
          DO j = 1, itht0
            j1 = j + 1
            DO k = 1, iphi0
              k1 = k + 1
              c1 = bphi(i,j,k) + bphi(i1,j,k)
     &          + bphi(i,j1,k) + bphi(i,j,k1)
     &          + bphi(i,j1,k1) + bphi(i1,j,k1)
     &          + bphi(i1,j1,k) + bphi(i1,j1,k1)
              c1 = c1/8.d0
              diff = ABS(c1 - bphi_test(i,j,k))
              IF (i==1 .AND. j == 1 .AND. k == 1) THEN
                maxabs  = diff
                locabs = c1
              ELSE
                maxabs = MAX(maxabs, diff)                              ! Maximum (absolute) discrepancy
                IF (maxabs == diff) locabs = c1                         ! Store value at which max. discrepancy happens
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      if(myid .eq. 0) then
        WRITE(*, 14) maxabs, locabs
      endif    !  if(myid .eq. 0) then
14      FORMAT('   - Max. absolute difference (Bphi):', f12.4,
     &    ' at value:', f12.4)
        DEALLOCATE(bphi_test)

        ALLOCATE(r_test(is0, itht0, iphi0), STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        r_test = zero
        CALL dbs3vd(0, 0, 0, isize, st, ut, vt, ks, ku, kv,
     &    s_knots, u_knots, v_knots, is, itht, iphi,
     &    r_coef, r_test, is_uniformx, is_uniformy,
     &    is_uniformz)

        DO i = 1, is0
          i1 = i + 1
          DO j = 1, itht0
            j1 = j + 1
            DO k = 1, iphi0
              k1 = k + 1
              c1 = r(i,j,k) + r(i1,j,k)
     &          + r(i,j1,k) + r(i,j,k1)
     &          + r(i,j1,k1) + r(i1,j,k1)
     &          + r(i1,j1,k) + r(i1,j1,k1)
              c1 = c1/8._dp
              diff = ABS(c1 - r_test(i,j,k))
              IF (i==1 .AND. j == 1 .AND. k == 1) THEN
                maxabs  = diff
                locabs = c1
              ELSE
                maxabs = MAX(maxabs, diff)                              ! Maximum (absolute) discrepancy
                IF (maxabs == diff) locabs = c1                         ! Store value at which max. discrepancy happens
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      if(myid .eq. 0) then
        WRITE(*, 15) maxabs, locabs
      endif    !  if(myid .eq. 0) then
15      FORMAT('   - Max. absolute difference (R):', f12.4,
     &    ' at value:', f12.4)
        DEALLOCATE(r_test)

        ALLOCATE(z_test(is0, itht0, iphi0), STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        z_test = zero
        CALL dbs3vd(0, 0, 0, isize, st, ut, vt, ks, ku, kv,
     &    s_knots, u_knots, v_knots, is, itht, iphi,
     &    z_coef, z_test, is_uniformx, is_uniformy,
     &    is_uniformz)

        DO i = 1, is0
          i1 = i + 1
          DO j = 1, itht0
            j1 = j + 1
            DO k = 1, iphi0
              k1 = k + 1
              c1 = z(i,j,k) + z(i1,j,k)
     &          + z(i,j1,k) + z(i,j,k1)
     &          + z(i,j1,k1) + z(i1,j,k1)
     &          + z(i1,j1,k) + z(i1,j1,k1)
              c1 = c1/8._dp
              diff = ABS(c1 - z_test(i,j,k))
              IF (i==1 .AND. j == 1 .AND. k == 1) THEN
                maxabs  = diff
                locabs = c1
              ELSE
                maxabs = MAX(maxabs, diff)                              ! Maximum (absolute) discrepancy
                IF (maxabs == diff) locabs = c1                         ! Store value at which max. discrepancy happens
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      if(myid .eq. 0) then
        WRITE(*, 16) maxabs, locabs
      endif    !  if(myid .eq. 0) then
16      FORMAT('   - Max. absolute difference (Z):', f12.4,
     &    ' at value:', f12.4)
        DEALLOCATE(z_test)

        ALLOCATE(phi_test(is0, itht0, iphi0), STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        phi_test = zero
        CALL dbs3vd(0, 0, 0, isize, st, ut, vt, ks, ku, kv,
     &    s_knots, u_knots, v_knots, is, itht, iphi,
     &    phi_coef, phi_test, is_uniformx, is_uniformy,
     &    is_uniformz)

        DO i = 1, is0
          i1 = i + 1
          DO j = 1, itht0
            j1 = j + 1
            DO k = 1, iphi0
              k1 = k + 1
              c1 = phi(i,j,k) + phi(i1,j,k)
     &          + phi(i,j1,k) + phi(i,j,k1)
     &          + phi(i,j1,k1) + phi(i1,j,k1)
     &          + phi(i1,j1,k) + phi(i1,j1,k1)
              c1 = c1/8._dp
              diff = ABS(c1 - phi_test(i,j,k))
              IF (i==1 .AND. j == 1 .AND. k == 1) THEN
                maxabs  = diff
                locabs = c1
              ELSE
                maxabs = MAX(maxabs, diff)                              ! Maximum (absolute) discrepancy
                IF (maxabs == diff) locabs = c1                         ! Store value at which max. discrepancy happens
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      if(myid .eq. 0) then
        WRITE(*, 17) maxabs, locabs
      endif    !  if(myid .eq. 0) then
17      FORMAT('   - Max. absolute difference (PHI):', f12.4,
     &    ' at value:', f12.4)
        DEALLOCATE(phi_test)

        DEALLOCATE(st, ut, vt)

      ENDIF
      DEALLOCATE(s, u, v, r, z, phi, br, bz, bphi)

      END SUBROUTINE spline_fields
