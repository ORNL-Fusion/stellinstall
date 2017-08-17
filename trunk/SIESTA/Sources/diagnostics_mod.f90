      MODULE diagnostics_mod
      USE fourier, ONLY: tomnsp, toijsp, tomnsp_par, toijsp_par,         &
                   f_sum, f_cos, f_sin, f_none, f_du, f_dv
      USE descriptor_mod, ONLY: iam, nprocs, SIESTA_COMM
      USE timer_mod
      USE quantities
      IMPLICIT NONE
      INTEGER     :: unit35=35, unit36=36
      REAL(dp) :: divb_rms, avbeta, bdotj_rms, divj_rms, bdotj2_rms
      REAL(dp) :: toroidal_flux, toroidal_flux0=0, bgradp_rms, wbgradp
      REAL(dp) :: max_bgradp, min_bgradp
      REAL(dp) :: tnorm_bgradp
!
!     The next three allocatable variables are kept in memory and used to
!     produce output whenever L_DUMP_DIAGNO = .TRUE. 
!      
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjmnch, bdotjijh
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjmnsh
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: divjmnsh, divbmnsf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: divjmnch, divbmncf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bgradpf, jpmnch_loc
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: fgradpmnch
      REAL(dp), DIMENSION(:), POINTER     :: xc, gc
      LOGICAL :: lextend = .FALSE.                   !Set in call to EVOLVE_BGRADP
      LOGICAL :: lcurr_init = .FALSE.
!
!     CALL THIS AFTER THE B-FIELD COMPONENTS bi = Jac*B-sup-i HAVE BEEN UPDATED
!     IN UPDATE_STATE (CALL FROM UPDATE_STATE IS SAFEST WAY!) 

      CONTAINS

        SUBROUTINE divb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER                  :: m, n, istat
        REAL(rprec)              :: ton, toff, tnorm
!-----------------------------------------------
!
!     COMPUTE divB ON FULL MESH POINTS
!
        CALL second0(ton)

        IF (.NOT. ALLOCATED(divbmnsf))  THEN
          ALLOCATE(divbmnsf(0:mpol,-ntor:ntor,ns),  stat=istat)          ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation failed in divb'
          divbmnsf = 0
        ENDIF

        DO m = 0, mpol
          DO n = -ntor,ntor
            divbmnsf(m,n,2:nsh) =                                        &
                   -m*(jbsupumnch(m,n,3:ns)+jbsupumnch(m,n,2:nsh))/2     &
               -n*nfp*(jbsupvmnch(m,n,3:ns)+jbsupvmnch(m,n,2:nsh))/2     &
               + ohs *(jbsupsmnsh(m,n,3:ns)-jbsupsmnsh(m,n,2:nsh))
          END DO
        END DO

        tnorm = hs_i*SUM(jbsupumnch**2 + jbsupvmnch**2 + jbsupsmnsh**2)
        divbmnsf = divbmnsf/SQRT(tnorm)
        divb_rms = SQRT(SUM(divbmnsf**2))

!        IF (iam .EQ. 0) WRITE(4000,*)'DIVB: tnorm=',SQRT(tnorm),' divb_rms: ',divb_rms
        CALL second0(toff)
        time_divb = time_divb + (toff-ton)

        END SUBROUTINE divb

        SUBROUTINE divb_par
        USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
        USE mpi_inc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER  :: istat, m, n, nsmin, nsmax, n1, n2
        REAL(dp) :: tnorm, s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE divB ON FULL MESH POINTS
!
        CALL second0(ton)

        nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
        IF (.NOT. ALLOCATED(divbmnsf)) THEN
        ALLOCATE(divbmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   ! Otherwise, allocate first time called
        IF (istat .NE. 0) STOP 'Allocation failed in divb'
        divbmnsf = 0
        END IF

        n1 = MAX(2, nsmin)
        n2 = MIN(ns-1, nsmax)
        DO m = 0, mpol
           DO n = -ntor, ntor
           divbmnsf(m,n,n1:n2) =                                        &
                -m*(jbsupumnch(m,n,n1+1:n2+1)+jbsupumnch(m,n,n1:n2))/2  &
            -n*nfp*(jbsupvmnch(m,n,n1+1:n2+1)+jbsupvmnch(m,n,n1:n2))/2  &
            + ohs *(jbsupsmnsh(m,n,n1+1:n2+1)-jbsupsmnsh(m,n,n1:n2))
           END DO
        END DO

        tnorm = hs_i*SUM(jbsupumnch(:,:,nsmin:nsmax)**2                 &
              +          jbsupvmnch(:,:,nsmin:nsmax)**2                 &
              +          jbsupsmnsh(:,:,nsmin:nsmax)**2)
       
        divb_rms = SUM(divbmnsf**2)
#if defined(MPI_OPT)        
        s1(1) = divb_rms; s1(2) = tnorm
        CALL MPI_ALLREDUCE(s1,r1,2,MPI_REAL8, MPI_SUM,                  &
                           SIESTA_COMM,MPI_ERR)
        divb_rms=r1(1);   tnorm = r1(2)
#endif
        IF (tnorm .NE. zero) THEN
           divb_rms = SQRT(divb_rms/tnorm)
           divbmnsf = divbmnsf/SQRT(tnorm)
        END IF

        CALL second0(toff)
        time_divb = time_divb + (toff-ton)

        END SUBROUTINE divb_par

         
        SUBROUTINE WRITE_PROFILES(fsq_total)
        USE safe_open_mod
        USE metrics, ONLY: lmns_i, lmnc_i
        USE fourier, ONLY: tomnsp_pest
        USE siesta_init, ONLY: init_state
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(rprec) :: fsq_total
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER                  :: istat, m, js
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: pijh, pmnch
        REAL(dp)                 :: tnorm, average(ns)
!-----------------------------------------------
!
!       IT IS ASSUMED COMING IN HERE THAT THE PROFILES HAVE BEEN GATHERED
!       (IN UPDATE_STATE)
!
        IF (iam .EQ. 0 .and. lverbose) PRINT *,'--- GENERALIZE WRITE_PROFILES FOR lasym=T ---'

        ALLOCATE(pijh(ntheta,nzeta,ns), stat=istat)
        ALLOCATE(pmnch(0:mpol,-ntor:ntor,ns), stat=istat)      
        IF (istat .ne. 0) STOP 'Allocation failed in WRITE_PROFILES' 

        CALL safe_open(unit35, istat, "siesta_profiles.txt", 'replace', &
                       'formatted')
        CALL safe_open(unit36, istat, "siesta_profiles_pest.txt",       &
                       'replace','formatted')

        IF (istat .NE. 0) STOP 'Error writing siesta profiles data file'

        pijh = 0
        CALL toijsp(jpmnch, pijh, f_none, f_cos)
        IF (lasym) THEN
           CALL toijsp(jpmnsh, pijh, f_sum, f_sin)
        END IF
        pijh = pijh/jacobh                                                 ! Need first to remove jacobian...
        pmnch = 0
        CALL tomnsp(pijh, pmnch, f_cos)
        CALL siesta_profiles (unit35,pmnch, fsq_total, 'pmnch(jrad);  fsq_total: ')
        IF (lasym) THEN
           CALL tomnsp(pijh, pmnch, f_sin)
           CALL siesta_profiles (unit35,pmnch, fsq_total, 'pmnsh(jrad);  fsq_total: ')
           CALL tomnsp_pest(pijh, pmnch, lmns_i, f_cos, lmnc=lmnc_i)
        ELSE
           CALL tomnsp_pest(pijh, pmnch, lmns_i, f_cos)
        END IF
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'pmnch(jrad);  fsq_total: ')
        IF (lasym) THEN
           CALL tomnsp_pest(pijh, pmnch, lmns_i, f_sin, lmnc=lmnc_i)
           CALL siesta_profiles (unit36, pmnch, fsq_total, 'pmnsh(jrad);  fsq_total: ')
        END IF

!USE pijh, pmnch to compute components of bsupX on half-mesh
!BSUPS
        CALL toijsp(jbsupsmnsh, pijh, f_none, f_sin)
        pijh = pijh/jacobh
        CALL tomnsp(pijh, pmnch, f_sin)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'B^s_mn;  fsq_total: ')
        CALL tomnsp_pest(pijh, pmnch, lmns_i, f_sin)
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'B^s_mn;  fsq_total: ')

!BSUPU
        CALL toijsp(jbsupumnch, pijh, f_none, f_cos)
        pijh = pijh/jacobh
        CALL tomnsp(pijh, pmnch, f_cos)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'B^u_mn;  fsq_total: ')
        CALL tomnsp_pest(pijh, pmnch, lmns_i, f_cos)
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'B^u_mn;  fsq_total: ')

!BSUPV
        CALL toijsp(jbsupvmnch, pijh, f_none, f_cos)
        pijh = pijh/jacobh
        CALL tomnsp(pijh, pmnch, f_cos)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'B^v_mn;  fsq_total: ')
        CALL tomnsp_pest(pijh, pmnch, lmns_i, f_cos)
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'B^v_mn;  fsq_total: ')

!KSUPS CRCook
        CALL toijsp(ksupsmnsf, pijh, f_none, f_sin)
        pijh = pijh/jacobh
        CALL tomnsp(pijh, pmnch, f_sin)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'J^s_mn;  fsq_total: ')
        CALL tomnsp_pest(pijh, pmnch, lmns_i, f_sin)
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'J^s_mn;  fsq_total: ')

!KSUPU CRCook
        CALL toijsp(ksupumncf, pijh, f_none, f_cos)
        pijh = pijh/jacobh
        CALL tomnsp(pijh, pmnch, f_cos)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'J^u_mn;  fsq_total: ')
        CALL tomnsp_pest(pijh, pmnch, lmns_i, f_cos)
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'J^u_mn;  fsq_total: ')

!KSUPV CRCook
        CALL toijsp(ksupvmncf, pijh, f_none, f_cos)
        pijh = pijh/jacobh
        CALL tomnsp(pijh, pmnch, f_cos)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'J^v_mn;  fsq_total: ')
        CALL tomnsp_pest(pijh, pmnch, lmns_i, f_cos)
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'J^v_mn;  fsq_total: ')

!K|| = (J*B)/B**2 (SPH:040915)
        lcurr_init=.TRUE.                             
        CALL init_state(lcurr_init)
        CALL bdotj(bdotjijh)
        pijh = 0
        WHERE (bsq .NE. zero) pijh = bdotjijh/bsq
        DEALLOCATE(bdotjijh)
        CALL tomnsp(pijh, pmnch, f_cos)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'K||_mn;  fsq_total: ')
        CALL tomnsp_pest(pijh, pmnch, lmns_i, f_cos)
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'K||_mn;  fsq_total: ')

!K|| (non-axisymmetric part only: SPH090816)
        CALL tomnsp(pijh, pmnch, f_cos)
        pmnch(:,0,:) = 0
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'K||_mn - nonaxisymmetric;  fsq_total: ')
        CALL tomnsp_pest(pijh, pmnch, lmns_i, f_cos)
        pmnch(:,0,:) = 0
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'K||_mn - nonaxisymmetric;  fsq_total: ')

!<K||-nonaxi> : RWILCOX metric in PEST coordinates
        CALL toijsp(pmnch, pijh, f_none, f_cos)
        CALL SurfAverage(average, ABS(pijh), 1, ns)
        WRITE (unit36, *) 
        WRITE (unit36, *) '<ABS(K||) - nonaxisymmetric>'
        DO js = 2,ns
           WRITE (unit36, '(i4,1p,e10.2)') js, average(js)
        END DO
        

!Displacements X jac (SPH:052114)
        CALL siesta_profiles (unit35, jvsupsmncf, fsq_total, 'vel^s_mn(jrad);  fsq_total: ')
        CALL siesta_profiles (unit35, jvsupumnsf, fsq_total, 'vel^u_mn(jrad);  fsq_total: ')
        CALL siesta_profiles (unit35, jvsupvmnsf, fsq_total, 'vel^v_mn(jrad);  fsq_total: ')

        DEALLOCATE (pijh, pmnch)

        CLOSE (unit35)
        CLOSE (unit36)

      END SUBROUTINE WRITE_PROFILES

      SUBROUTINE SIESTA_PROFILES (iunit, value, fsq_total, label)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)       :: iunit
      REAL(rprec), INTENT(in)   :: value(0:mpol,-ntor:ntor,ns)
      REAL(rprec), INTENT(in)   :: fsq_total
      CHARACTER*(*), INTENT(in) :: label
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, m, js
      INTEGER :: mnmax
      CHARACTER (len=19) :: p_format
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: m_array
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: n_array
!-----------------------------------------------

      WRITE (iunit, *)
      WRITE (iunit, *) TRIM(label), fsq_total

      ALLOCATE(m_array(0:mpol,-ntor:ntor))
      ALLOCATE(n_array(0:mpol,-ntor:ntor))

      DO n = -ntor, ntor
         DO m = 0, mpol
            m_array(m,n) = m
            n_array(m,n) = n
         END DO
      END DO

      mnmax = (mpol + 1)*(2*ntor + 1)

      WRITE (p_format,1000) mnmax
      WRITE (iunit,p_format) 'MPOL', m_array
      WRITE (p_format,1001) mnmax
      WRITE (iunit,p_format) 'NTOR', n_array

      WRITE (p_format,1002) mnmax
      DO js = 2, ns
         WRITE(iunit,p_format) js, value(:,:,js)
      END DO

      DEALLOCATE(m_array, n_array)

1000  FORMAT('(a,',i3,'(2x,i12))')
1001  FORMAT('(a,',i3,'(2x,i12))')
1002  FORMAT('(i4,',i3,'(2x,es12.5))')

      END SUBROUTINE SIESTA_PROFILES


      SUBROUTINE bgradp
!     
!     WRITTEN 08-16-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes B*GRAD P (normalized to B*P) on the FULL mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                  :: istat
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      ALLOCATE(bgradpf(ntheta,nzeta,ns), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in BGRADP'

      CALL get_bgradp

      DEALLOCATE (bgradpf)

      END SUBROUTINE bgradp


      SUBROUTINE get_bgradp
!     
!     PURPOSE: Computes B*GRAD P on the FULL mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: istat, n1
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                       &
          dpduijh, dpdvijh, pijh,                                    &
          bsupsijh0, bsupuijh0, bsupvijh0
      REAL(dp)    :: ton, toff
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      CALL second0(ton)

      bgradpf = 0 

!     Compute pressure and angle derivatives
!
      ALLOCATE(pijh(ntheta,nzeta,ns),                                  &
               dpduijh(ntheta,nzeta,ns),                               &
               dpdvijh(ntheta,nzeta,ns),                               &
               bsupsijh0(ntheta,nzeta,ns),                             &
               bsupuijh0(ntheta,nzeta,ns),                             &
               bsupvijh0(ntheta,nzeta,ns), stat=istat)                         
      IF (istat .ne. 0) STOP 'Allocation failed in GET_BGRADP'

      CALL GetPressure (pijh, dpduijh, dpdvijh)
! 
!     Compute contravariant components of magnetic field
!           
      IF (.NOT.ALLOCATED(jbsupsmnsh)) STOP 'jbsupXmn UNALLOCATED IN GET_BGRADP'
      CALL toijsp(jbsupsmnsh, bsupsijh0, f_none, f_sin)
      CALL toijsp(jbsupumnch, bsupuijh0, f_none, f_cos)
      CALL toijsp(jbsupvmnch, bsupvijh0, f_none, f_cos)

	  IF (.NOT.lasym) GOTO 100
      IF (.NOT.ALLOCATED(jbsupsmnch)) STOP 'jbsupXmn UNALLOCATED IN GET_BGRADP'
      CALL toijsp(jbsupsmnch, bsupsijh0, f_sum, f_cos)
      CALL toijsp(jbsupumnsh, bsupuijh0, f_sum, f_sin)
      CALL toijsp(jbsupvmnsh, bsupvijh0, f_sum, f_sin)
100   CONTINUE

      bsupsijh0 = bsupsijh0/jacobh
      bsupuijh0 = bsupuijh0/jacobh
      bsupvijh0 = bsupvijh0/jacobh

      bgradpf(:,:,2:nsh) =                                              &   ! Get all values except first and
               (bsupsijh0(:,:,3:ns) + bsupsijh0(:,:,2:nsh))           &   ! last surface from half mesh values.
               *    (pijh(:,:,3:ns) -      pijh(:,:,2:nsh))*ohs/2     &
             + (bsupuijh0(:,:,3:ns) + bsupuijh0(:,:,2:nsh))           &
              *(  dpduijh(:,:,3:ns) +   dpduijh(:,:,2:nsh))/4         &
             + (bsupvijh0(:,:,3:ns) + bsupvijh0(:,:,2:nsh))           &
              *(  dpdvijh(:,:,3:ns) +   dpdvijh(:,:,2:nsh))/4                   
    
      tnorm_bgradp = SUM((bsupsijh0**2 + bsupuijh0**2 + bsupvijh0**2)*pijh**2)
      IF (tnorm_bgradp .eq. zero) tnorm_bgradp = EPSILON(tnorm_bgradp)

      wbgradp = SUM(bgradpf(:,:,3:nsh)**2)
      bgradp_rms = SQRT(wbgradp/tnorm_bgradp)

      tnorm_bgradp = SQRT(tnorm_bgradp/(ns*ntheta*nzeta))

      bgradpf(:,:,2:nsh) = bgradpf(:,:,2:nsh)/tnorm_bgradp
!IGNORE FIRST POINT
      n1 = MAX(3, ns/20)
      max_bgradp = MAXVAL(bgradpf(:,:,n1:nsh))
      min_bgradp = MINVAL(bgradpf(:,:,n1:nsh))

      DEALLOCATE (pijh, dpdvijh, dpduijh,                            &
                  bsupsijh0, bsupuijh0, bsupvijh0)

      CALL second0(toff)
      time_bgradp = time_bgradp+(toff-ton)

      END SUBROUTINE get_bgradp

      SUBROUTINE bgradp_par
!     
!     WRITTEN 03-21-13 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes B*GRAD P (normalized to B*P) on the FULL mesh
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE mpi_inc
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: istat, nsmin, nsmax
      REAL(dp)    :: s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
      
      ALLOCATE(bgradpf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in BGRADP'
      CALL get_bgradp_par

      DEALLOCATE (bgradpf, stat=istat)

      END SUBROUTINE bgradp_par


      SUBROUTINE get_bgradp_par
      USE nscalingtools, ONLY: startglobrow, endglobrow, mpi_err
      USE mpi_inc
!     
!     PURPOSE: Computes B*GRAD P on the FULL mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER  :: jsh=3
      INTEGER     :: istat, nsmin, nsmax, ns_span, n1, n2
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                        &
          dpduijh, dpdvijh, pijh,                                        &
          bsupsijh0, bsupuijh0, bsupvijh0
      REAL(dp)    :: s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      CALL second0(ton)

      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
      ns_span = nsmax-nsmin+1

      bgradpf = 0 

!     Compute pressure and angle derivatives
!
      ALLOCATE(pijh(ntheta,nzeta,nsmin:nsmax),                           &
               dpduijh(ntheta,nzeta,nsmin:nsmax),                        &
               dpdvijh(ntheta,nzeta,nsmin:nsmax),                        &
               bsupsijh0(ntheta,nzeta,nsmin:nsmax),                      &
               bsupuijh0(ntheta,nzeta,nsmin:nsmax),                      &
               bsupvijh0(ntheta,nzeta,nsmin:nsmax), stat=istat)                         
      IF (istat .NE. 0) STOP 'Allocation failed in GET_BGRADP_PAR'

      IF (.NOT.ALLOCATED(jbsupsmnsh)) STOP 'jbsupXmn UNALLOCATED IN GET_BGRADP_PAR'

      CALL GetPressure_par (pijh, dpduijh, dpdvijh, nsmin, nsmax)
! 
!     Compute contravariant components of magnetic field
!           
      CALL toijsp_par(jbsupsmnsh(:,:,nsmin:nsmax), bsupsijh0, f_none, f_sin,1, ns_span)
      CALL toijsp_par(jbsupumnch(:,:,nsmin:nsmax), bsupuijh0, f_none, f_cos,1, ns_span)
      CALL toijsp_par(jbsupvmnch(:,:,nsmin:nsmax), bsupvijh0, f_none, f_cos,1, ns_span)

      IF (.NOT.lasym) GOTO 100
      CALL toijsp_par(jbsupsmnch(:,:,nsmin:nsmax), bsupsijh0, f_sum, f_cos,1, ns_span)
      CALL toijsp_par(jbsupumnsh(:,:,nsmin:nsmax), bsupuijh0, f_sum, f_sin,1, ns_span)
      CALL toijsp_par(jbsupvmnsh(:,:,nsmin:nsmax), bsupvijh0, f_sum, f_sin,1, ns_span)

 100  CONTINUE

      bsupsijh0 = bsupsijh0/jacobh(:,:,nsmin:nsmax)
      bsupuijh0 = bsupuijh0/jacobh(:,:,nsmin:nsmax)
      bsupvijh0 = bsupvijh0/jacobh(:,:,nsmin:nsmax)

      nsmax=MIN(endglobrow,ns)
      n1 = MAX(2, nsmin); n2 = MIN(nsh, nsmax)
      bgradpf(:,:,n1:n2) =                                               &   ! Get all values except first and
               (bsupsijh0(:,:,n1+1:n2+1) + bsupsijh0(:,:,n1:n2))         &   ! last surface from half mesh values.
             *ohs/2*(pijh(:,:,n1+1:n2+1) -      pijh(:,:,n1:n2))         &
             + (bsupuijh0(:,:,n1+1:n2+1) + bsupuijh0(:,:,n1:n2))         &
              *(  dpduijh(:,:,n1+1:n2+1) +   dpduijh(:,:,n1:n2))/4       &
             + (bsupvijh0(:,:,n1+1:n2+1) + bsupvijh0(:,:,n1:n2))         &
              *(  dpdvijh(:,:,n1+1:n2+1) +   dpdvijh(:,:,n1:n2))/4                   
    
      tnorm_bgradp = SUM((bsupsijh0(:,:,n1:nsmax)**2                     &
                   +      bsupuijh0(:,:,n1:nsmax)**2                     &
                   +      bsupvijh0(:,:,n1:nsmax)**2)                    &
                   *           pijh(:,:,n1:nsmax)**2)

!SPH10-31-16: ignore inner s<.1 for this diagnostic
      n1 = MAX(ns/20, nsmin)
      wbgradp = SUM(bgradpf(:,:, n1:n2)**2)
      s1(1)=tnorm_bgradp;  s1(2)=wbgradp
#if defined(SKS)
      CALL MPI_ALLREDUCE(s1,r1,2,MPI_REAL8, MPI_SUM,                     &
                         SIESTA_COMM,MPI_ERR)
#else
      r1 = s1
#endif
      tnorm_bgradp=r1(1); wbgradp=r1(2)

      IF (tnorm_bgradp .GT. zero) THEN
         bgradp_rms = SQRT(wbgradp/tnorm_bgradp)
         tnorm_bgradp = SQRT(tnorm_bgradp/(ns*nuv))

!         IF (iam .EQ. 0) WRITE (3200, *)' tnorm_bgradp: ',tnorm_bgradp,  &
!                         ' bgradp_rms: ', bgradp_rms

!IGNORE FIRST POINT
         max_bgradp = MAXVAL(bgradpf(:,:,n1:n2))/tnorm_bgradp
         min_bgradp = MINVAL(bgradpf(:,:,n1:n2))/tnorm_bgradp
         s1(1)=max_bgradp; s1(2) = -min_bgradp
#if defined(SKS)
         CALL MPI_ALLREDUCE(s1,r1,2,MPI_REAL8, MPI_MAX,                  &
                            SIESTA_COMM,MPI_ERR)
#else
         r1 = s1
#endif
         max_bgradp=r1(1)
         min_bgradp=-r1(2)
         IF (.NOT.lasym) THEN                                              !reflect onto u=pi-2*pi: bgradp is sin parity
            max_bgradp = MAX(max_bgradp, -min_bgradp)
            min_bgradp = MIN(min_bgradp, -max_bgradp)
         END IF
      END IF

      DEALLOCATE (pijh, dpdvijh, dpduijh,                                &
                  bsupsijh0, bsupuijh0, bsupvijh0, stat=istat)

      CALL second0(toff)
      time_bgradp = time_bgradp+(toff-ton)

      END SUBROUTINE get_bgradp_par
      

      SUBROUTINE GetPressure (pijh, dpduijh, dpdvijh)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(OUT)  ::              &
           dpduijh, dpdvijh, pijh
!-----------------------------------------------

      CALL GetPressure_par(pijh, dpduijh, dpdvijh, 1, ns)

      END SUBROUTINE GetPressure

      SUBROUTINE GetPressure_par (pijh, dpduijh, dpdvijh, nsmin, nsmax)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: nsmin, nsmax
      REAL(dp), DIMENSION(ntheta,nzeta, nsmin:nsmax),                    &
               INTENT(OUT)  :: dpduijh, dpdvijh, pijh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: istat, ns_span, m1, m2
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: pmnch, pmnsh
!-----------------------------------------------
      ns_span = nsmax-nsmin+1

      CALL toijsp_par(jpmnch(:,:,nsmin:nsmax), pijh, f_none, f_cos,      &
                      1, ns_span)
      IF (lasym) THEN
      CALL toijsp_par(jpmnsh(:,:,nsmin:nsmax), pijh, f_sum,  f_sin,      &
                      1, ns_span)
      END IF
      pijh = pijh/jacobh(:,:,nsmin:nsmax)                                  ! Remove jacobian factor
      
      ALLOCATE(pmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)      
      IF (istat .NE. 0) STOP 'Allocation error in GetPressure_par'
      CALL tomnsp_par(pijh, pmnch, f_cos,1, ns_span)

      CALL toijsp_par(pmnch, pijh, f_none, f_cos,1, ns_span)
      CALL toijsp_par(pmnch, dpduijh, f_du, f_cos,1, ns_span)              ! Get theta and phi pressure derivatives
      CALL toijsp_par(pmnch, dpdvijh, f_dv, f_cos,1, ns_span)

      IF (lasym) THEN
      ALLOCATE(pmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)      
      IF (istat .NE. 0) STOP 'Allocation error in GetPressure_par'
      CALL tomnsp_par(pijh, pmnsh, f_sin,1, ns_span)

      m1 = IOR(f_sum,f_du);  m2 = IOR(f_sum,f_dv)
      CALL toijsp_par(pmnsh, pijh, f_sum, f_sin,1, ns_span)
      CALL toijsp_par(pmnsh, dpduijh, m1, f_sin,1, ns_span)                ! Get theta and phi pressure derivatives
      CALL toijsp_par(pmnsh, dpdvijh, m2, f_sin,1, ns_span)
      END IF

     
      DEALLOCATE (pmnch)
      IF (lasym) DEALLOCATE (pmnsh)

      END SUBROUTINE GetPressure_par

      SUBROUTINE BETA

        avbeta = wp/wb
 
      END SUBROUTINE BETA
 
      
      SUBROUTINE TFLUX
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, n0=0
!-----------------------------------------------

!Averages over all toroidal cross sections (which should be the same)
      toroidal_flux = SUM(jbsupvmnch(m0,n0,2:ns))

      toroidal_flux = signjac*twopi*toroidal_flux*hs_i/b_factor  
      
      IF (toroidal_flux0 .EQ. zero) toroidal_flux0 = toroidal_flux

      END SUBROUTINE TFLUX


      SUBROUTINE bdotj (bdotjijhA)
!
!     WRITTEN 08-20-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes parallel current normalized to the total current (J*B/|JxB|) on the half mesh
!-----------------------------------------------
!   D u m m y   A r g u m e n t
!-----------------------------------------------
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjijhA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER      :: istat
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupsijh,             &
          bsupuijh, bsupvijh,                                            &
          bsubsijh, bsubuijh, bsubvijh, ksubsijh, ksubuijh,              &
          ksubvijh, ksupsijh, ksupuijh, ksupvijh, work1, work2,          &
          work3, work4, work5, work6
        REAL(dp)  :: tnorm, tnorm2, ton, toff
!-----------------------------------------------
        CALL second0(ton)

        IF (.NOT. lcurr_init) STOP 'MUST CALL init_state(.TRUE.) BEFORE bdotj'
! 
!     Compute contravariant components of magnetic field
!       
        ALLOCATE(bsupsijh(ntheta,nzeta,ns),                              &
                 bsupuijh(ntheta,nzeta,ns),                              &
                 bsupvijh(ntheta,nzeta,ns), stat=istat)
        bsupsijh=0; bsupuijh=0; bsupvijh=0
        IF (istat .NE. 0) STOP 'Allocation #1 failed in BDOTJ' 

        CALL toijsp(jbsupsmnsh, bsupsijh, f_none, f_sin)
        CALL toijsp(jbsupumnch, bsupuijh, f_none, f_cos)
        CALL toijsp(jbsupvmnch, bsupvijh, f_none, f_cos)
        IF (lasym) THEN
        CALL toijsp(jbsupsmnch, bsupsijh, f_sum, f_sin)
        CALL toijsp(jbsupumnsh, bsupuijh, f_sum, f_cos)
        CALL toijsp(jbsupvmnsh, bsupvijh, f_sum, f_cos)
        END IF
        bsupsijh = bsupsijh/jacobh
        bsupuijh = bsupuijh/jacobh
        bsupvijh = bsupvijh/jacobh
!
!   Compute covariant field. Used for normalization.
!        
        ALLOCATE(bsubsijh(ntheta,nzeta,ns),                              &
                 bsubuijh(ntheta,nzeta,ns),                              &
                 bsubvijh(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #2 failed in BDOTJ' 
        bsubsijh=0; bsubvijh=0; bsubuijh=0
      
        CALL tolowerh(bsupsijh, bsupuijh, bsupvijh,                      &
                      bsubsijh, bsubuijh, bsubvijh)

!
!    In initial state, calculate currents (remember they are multiplied by a jacobian!)
!
!
!    Move contravariant currents (*jacob) to the half mesh
!       
        ALLOCATE(ksupsijh(ntheta,nzeta,ns),                              &
                 ksupuijh(ntheta,nzeta,ns),                              &
                 ksupvijh(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #3 failed in BDOTJ'

        CALL to_half_mesh(ksupsijf0, ksupsijh, nuv)
        CALL to_half_mesh(ksupuijf0, ksupuijh, nuv)
        CALL to_half_mesh(ksupvijf0, ksupvijh, nuv)
!
!   Remove jacobian 
!        
        ksupsijh = ksupsijh/jacobh      
        ksupuijh = ksupuijh/jacobh  
        ksupvijh = ksupvijh/jacobh
!          
!   Get covariant currents
!  
        ALLOCATE(ksubsijh(ntheta,nzeta,ns),                              &
                 ksubuijh(ntheta,nzeta,ns),                              &
                 ksubvijh(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #4 failed in BDOTJ'
        ksubsijh=0; ksubuijh=0; ksubvijh=0

        CALL tolowerh(ksupsijh, ksupuijh, ksupvijh,                      &
                      ksubsijh, ksubuijh, ksubvijh)
!
!   Compute B*J in the half mesh.
!        
        IF (ALLOCATED(bdotjmnch) .AND. SIZE(bdotjmnch,3).NE. ns)         &
            DEALLOCATE(bdotjmnch)
        IF (.NOT. ALLOCATED(bdotjmnch)) THEN                             ! Should already be allocated in QUANTITIES       
          ALLOCATE(bdotjmnch(0:mpol,-ntor:ntor,ns), stat=istat)          ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation #5 failed in BDOTJ'
        ENDIF
        IF (lasym) THEN
        IF (ALLOCATED(bdotjmnsh) .AND. SIZE(bdotjmnsh,3).NE. ns)         &
            DEALLOCATE(bdotjmnsh)
        IF (.NOT. ALLOCATED(bdotjmnsh)) THEN                             ! Should already be allocated in QUANTITIES       
          ALLOCATE(bdotjmnsh(0:mpol,-ntor:ntor,ns), stat=istat)          ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation #5 failed in BDOTJ'
        ENDIF
        END IF
        ALLOCATE(bdotjijhA(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #6 failed in BDOTJ'
        bdotjijhA=0
        
        bdotjijhA(:,:,3:nsh-1) =                                         &  ! RS: If I include 2 and nsh, it gets huge
          bsupsijh(:,:,3:nsh-1)*ksubsijh(:,:,3:nsh-1) +                  &  ! something funny happens at the boundaries.
          bsupuijh(:,:,3:nsh-1)*ksubuijh(:,:,3:nsh-1) +                  &
          bsupvijh(:,:,3:nsh-1)*ksubvijh(:,:,3:nsh-1)                       ! B*J = (B^s*J_s + B^u*J_u + B^v*J_v)
! 
!   Compute tnorm = |JXB|^2 for normalization; thus we print out J_parallel/|J_perp|.
!
        ALLOCATE(work1(ntheta,nzeta,ns), work2(ntheta,nzeta,ns),         &
                 work3(ntheta,nzeta,ns), work4(ntheta,nzeta,ns),         &
                 work5(ntheta,nzeta,ns), work6(ntheta,nzeta,ns),         &
                 stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #7 failed in BDOTJ'
        work1=0; work2=0; work3=0
        work4=0; work5=0; work6=0
        
        work1 = (bsubuijh*ksubvijh - bsubvijh*ksubuijh)/jacobh              ! (JxB)^s  = (1/sqrt(g))*(K_uB_b-K_vB_u)
        work2 = (bsubvijh*ksubsijh - bsubsijh*ksubvijh)/jacobh              ! (JxB)^u  = .....
        work3 = (bsubsijh*ksubuijh - bsubuijh*ksubsijh)/jacobh              ! (JxB)^v  = .....
        CALL tolowerh(work1, work2, work3, work4, work5, work6)
        
        work1 = work1*work4 + work2*work5 + work3*work6                     ! |JxB|**2
        tnorm = SUM(work1(:,:,3:nsh-1))
        IF (iam.EQ.0 .AND. ANY(work1 .LT. zero))                         &
           PRINT *,'ERROR1 IN BDOTJ, JXB^2 < 0!'
        tnorm = SQRT(ABS(tnorm))
        DEALLOCATE(work3, work4, work5, work6)
        
        work1 = bsupsijh*bsubsijh + bsupuijh*bsubuijh +                  &  ! |B|**2
                bsupvijh*bsubvijh
        work2 = ksupsijh*ksubsijh + ksupuijh*ksubuijh +                  &  ! |J|**2
                ksupvijh*ksubvijh
       
        tnorm2 = SUM(work1(:,:,3:nsh-1)*work2(:,:,3:nsh-1))
        work1(:,:,1) = 0;  work2(:,:,1) = 0
        IF (iam.EQ.0 .AND. ANY(work1.LT.zero))                           &
           PRINT *, 'ERROR2 IN BDOTJ: B^2 < 0!'               
        IF (iam.EQ.0 .AND. ANY(work2.LT.zero))                           &
           PRINT *, 'ERROR3 IN BDOTJ: J^2 < 0!'               
        tnorm2 = SQRT(ABS(tnorm2))                                           ! |J|*|B|
        bdotj_rms = SUM(bdotjijhA*bdotjijhA)

        DEALLOCATE (bsupsijh, bsupuijh, bsupvijh, bsubsijh,              &
                    bsubuijh, bsubvijh, ksupsijh, ksupuijh,              &
                    ksupvijh, ksubsijh, ksubuijh, ksubvijh,              &      
                    work1, work2)

        tnorm = MAX(tnorm, EPSILON(tnorm))
        tnorm2= MAX(tnorm2,EPSILON(tnorm2))
        bdotj2_rms = SQRT(ABS(bdotj_rms))/tnorm2                           ! RMS of bdotj/|J|*|B|  
        bdotj_rms  = SQRT(ABS(bdotj_rms))/tnorm                            ! RMS of bdotj/|JxB|  
        CALL tomnsp(bdotjijhA, bdotjmnch, f_cos)                  ! Keep harmonics of BDOTJ for output
        bdotjmnch = bdotjmnch/tnorm
        IF (lasym) THEN
        CALL tomnsp(bdotjijhA, bdotjmnsh, f_cos)                  ! Keep harmonics of BDOTJ for output
        bdotjmnsh = bdotjmnsh/tnorm
        END IF
        CALL second0(toff)
        time_bdotj = time_bdotj+(toff-ton)

        bdotjijh = bdotjijhA

        END SUBROUTINE bdotj

      SUBROUTINE bdotj_par
!#if defined(SKS)
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE mpi_inc
!     
!     WRITTEN 08-20-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     UPDATED 03-21-13 BY S. HIRSHMAN FOR PARALLEL EXECUTION
!     
!     PURPOSE: Computes parallel current normalized to the total current (J*B/|JxB|) on the half mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: istat, nsmin, nsmax, ns_span, n1, n2
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupsijh,               &
          bsupuijh, bsupvijh,                                            &
          bsubsijh, bsubuijh, bsubvijh, ksubsijh, ksubuijh,              &
          ksubvijh, ksupsijh, ksupuijh, ksupvijh, work1, work2,          &
          work3, work4, work5, work6
      REAL(dp)  :: tnorm, tnorm2, s1(3), r1(3), ton, toff
!-----------------------------------------------
      CALL second0(ton)

      IF (.NOT. lcurr_init) STOP 'MUST CALL init_state_par(.TRUE.) BEFORE bdotj_par'
      nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow)
      ns_span=nsmax-nsmin+1
! 
!     Compute contravariant components of magnetic field
!       
        ALLOCATE(bsupsijh(ntheta,nzeta,nsmin:nsmax),                     &
                 bsupuijh(ntheta,nzeta,nsmin:nsmax),                     &
                 bsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #1 failed in BDOTJ_PAR' 

        CALL toijsp_par(jbsupsmnsh(:,:,nsmin:nsmax), bsupsijh, f_none, f_sin,1, ns_span)
        CALL toijsp_par(jbsupumnch(:,:,nsmin:nsmax), bsupuijh, f_none, f_cos,1, ns_span)
        CALL toijsp_par(jbsupvmnch(:,:,nsmin:nsmax), bsupvijh, f_none, f_cos,1, ns_span)
        IF (lasym) THEN
        CALL toijsp_par(jbsupsmnch(:,:,nsmin:nsmax), bsupsijh, f_sum, f_sin,1, ns_span)
        CALL toijsp_par(jbsupumnsh(:,:,nsmin:nsmax), bsupuijh, f_sum, f_cos,1, ns_span)
        CALL toijsp_par(jbsupvmnsh(:,:,nsmin:nsmax), bsupvijh, f_sum, f_cos,1, ns_span)
        END IF

        bsupsijh = bsupsijh/jacobh(:,:,nsmin:nsmax)
        bsupuijh = bsupuijh/jacobh(:,:,nsmin:nsmax)
        bsupvijh = bsupvijh/jacobh(:,:,nsmin:nsmax)
!
!   Compute covariant field. Used for normalization.
!        
        ALLOCATE(bsubsijh(ntheta,nzeta,nsmin:nsmax),                     &
                 bsubuijh(ntheta,nzeta,nsmin:nsmax),                     &
                 bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #2 failed in BDOTJ_PAR' 
      
        CALL tolowerh_par(bsupsijh, bsupuijh, bsupvijh,                  &
                          bsubsijh, bsubuijh, bsubvijh, nsmin, nsmax)
!
!    In initial state, calculate currents (remember they are multiplied by a jacobian!)
!
!
!    Move contravariant currents (*jacob) to the half mesh
!       
        nsmin=MAX(1,startglobrow-1)
        ALLOCATE(ksupsijh(ntheta,nzeta,nsmin:nsmax),                     &
                 ksupuijh(ntheta,nzeta,nsmin:nsmax),                     &
                 ksupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #3 failed in BDOTJ_PAR'

        CALL to_half_mesh_par(ksupsijf0, ksupsijh, nuv, nsmin, nsmax)
        CALL to_half_mesh_par(ksupuijf0, ksupuijh, nuv, nsmin, nsmax)
        CALL to_half_mesh_par(ksupvijf0, ksupvijh, nuv, nsmin, nsmax)

!
!   Remove jacobian 
!        
        ksupsijh = ksupsijh/jacobh(:,:,nsmin:nsmax)
        ksupuijh = ksupuijh/jacobh(:,:,nsmin:nsmax)
        ksupvijh = ksupvijh/jacobh(:,:,nsmin:nsmax)

!          
!   Get covariant currents
!  
        ALLOCATE(ksubsijh(ntheta,nzeta,nsmin:nsmax),                     &
                 ksubuijh(ntheta,nzeta,nsmin:nsmax),                     &
                 ksubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #4 failed in BDOTJ_PAR'

        CALL tolowerh_par(ksupsijh, ksupuijh, ksupvijh,                  &
                          ksubsijh, ksubuijh, ksubvijh, nsmin, nsmax)
        nsmin=MAX(1,startglobrow)

!
!   Compute B*J in the half mesh.
!      
        IF (ALLOCATED(bdotjmnch) .AND. SIZE(bdotjmnch,3).NE.ns_span)     &
            DEALLOCATE(bdotjmnch)
        IF (.NOT. ALLOCATED(bdotjmnch)) THEN                            
          ALLOCATE(bdotjmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation #5A failed in BDOTJ_PAR'
        ENDIF
        IF (lasym) THEN
        IF (ALLOCATED(bdotjmnsh) .AND. SIZE(bdotjmnsh,3).NE.ns_span)     &
            DEALLOCATE(bdotjmnsh)
        IF (.NOT. ALLOCATED(bdotjmnsh)) THEN                            
          ALLOCATE(bdotjmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation #6B failed in BDOTJ_PAR'
        ENDIF
        ENDIF

        ALLOCATE(bdotjijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #6 failed in BDOTJ_PAR'
        bdotjijh=0
        
        n1 = MAX(3,nsmin); n2 = MIN(nsh-1, nsmax)
        bdotjijh(:,:,n1:n2) =                                            &  ! RS: If include 2 and nsh, it gets huge
          bsupsijh(:,:,n1:n2)*ksubsijh(:,:,n1:n2) +                      &  ! something funny happens at the boundaries.
          bsupuijh(:,:,n1:n2)*ksubuijh(:,:,n1:n2) +                      &
          bsupvijh(:,:,n1:n2)*ksubvijh(:,:,n1:n2)                           ! B*J = (B^s*J_s + B^u*J_u + B^v*J_v)
        
! 
!   Compute tnorm = |JXB|^2 for normalization; thus we print out J_parallel/|J_perp|.
!
        ALLOCATE(work1(ntheta,nzeta,nsmin:nsmax),                        &
                 work2(ntheta,nzeta,nsmin:nsmax),                        &
                 work3(ntheta,nzeta,nsmin:nsmax),                        &
                 work4(ntheta,nzeta,nsmin:nsmax),                        &
                 work5(ntheta,nzeta,nsmin:nsmax),                        &
                 work6(ntheta,nzeta,nsmin:nsmax),                        &
                 stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #7 failed in BDOTJ_PAR'
        
        work1 = bsubuijh*ksubvijh(:,:,nsmin:nsmax)                       &
              - bsubvijh*ksubuijh(:,:,nsmin:nsmax)                        ! (JxB)^s  = (1/sqrt(g))*(K_uB_b-K_vB_u)
        work2 = bsubvijh*ksubsijh(:,:,nsmin:nsmax)                       & 
              - bsubsijh*ksubvijh(:,:,nsmin:nsmax)                        ! (JxB)^u  = .....
        work3 = bsubsijh*ksubuijh(:,:,nsmin:nsmax)                       &
              - bsubuijh*ksubsijh(:,:,nsmin:nsmax)                        ! (JxB)^v  = .....
        work1 = work1/jacobh(:,:,nsmin:nsmax)
        work2 = work2/jacobh(:,:,nsmin:nsmax)
        work3 = work3/jacobh(:,:,nsmin:nsmax)
        CALL tolowerh_par(work1, work2, work3,                           &
                          work4, work5, work6, nsmin, nsmax)
        
        work1 = work1*work4 + work2*work5 + work3*work6                     ! |JxB|**2
        tnorm = SUM(work1(:,:,n1:n2))
        IF (nsmin .EQ. 1) work1(:,:,1) = 0
        IF (iam.EQ.0 .AND. ANY(work1 .LT. zero))                         &
           PRINT *,'ERROR1 IN BDOTJ_PAR, JXB^2 < 0!'
        DEALLOCATE(work3, work4, work5, work6)
        
        work1 = bsupsijh*bsubsijh + bsupuijh*bsubuijh +                  &   ! |B|**2
                bsupvijh*bsubvijh
        work2 = ksupsijh(:,:,nsmin:nsmax)*ksubsijh(:,:,nsmin:nsmax)      &
              + ksupuijh(:,:,nsmin:nsmax)*ksubuijh(:,:,nsmin:nsmax)      &
              + ksupvijh(:,:,nsmin:nsmax)*ksubvijh(:,:,nsmin:nsmax)         ! |J|**2
       
        IF (nsmin .EQ. 1) THEN
           work1(:,:,1) = 0;  work2(:,:,1) = 0
        END IF
        IF (iam.EQ.0 .AND. ANY(work1.LT.zero))                           &
           PRINT *, 'ERROR2 IN BDOTJ_PAR: B^2 < 0!'               
        IF (iam.EQ.0 .AND. ANY(work2.LT.zero))                           &
           PRINT *, 'ERROR3 IN BDOTJ_PAR: J^2 < 0!'               
        tnorm2 = SUM(work1(:,:,n1:n2)*work2(:,:,n1:n2))
        bdotj_rms = SUM(bdotjijh(:,:,n1:n2)*bdotjijh(:,:,n1:n2))
#if defined(MPI_OPT)
        s1(1)=tnorm; s1(2)=tnorm2; s1(3)=bdotj_rms
        CALL MPI_ALLREDUCE(s1,r1,3,MPI_REAL8, MPI_SUM,                   &
                           SIESTA_COMM,MPI_ERR)
        tnorm  = SQRT(ABS(r1(1)))
        tnorm2 = SQRT(ABS(r1(2)))                                             ! |J|*|B|
        bdotj_rms = r1(3) 
#else
        tnorm = SQRT(ABS(tnorm)); tnorm2 = SQRT(ABS(tnorm2))
#endif
        DEALLOCATE (bsupsijh, bsupuijh, bsupvijh, bsubsijh,              &
                    bsubuijh, bsubvijh, ksupsijh, ksupuijh,              &
                    ksupvijh, ksubsijh, ksubuijh, ksubvijh,              &      
                    work1, work2)


        tnorm2= MAX(tnorm2,EPSILON(tnorm2))
        bdotj2_rms = SQRT(ABS(bdotj_rms))/tnorm2                            ! RMS of bdotj/|J|*|B|  
        bdotj_rms  = SQRT(ABS(bdotj_rms))/tnorm                             ! RMS of bdotj/|JxB|  

        bdotjijh = bdotjijh/tnorm
        
        CALL tomnsp_par(bdotjijh, bdotjmnch, f_cos,1, ns_span)     ! Keep harmonics of BDOTJ for output
        IF (lasym) THEN
        CALL tomnsp_par(bdotjijh, bdotjmnsh, f_cos,1, ns_span)     ! Keep harmonics of BDOTJ for output
        END IF

        DEALLOCATE (bdotjijh)

        CALL second0(toff)
        time_bdotj = time_bdotj+(toff-ton)
!#endif
        END SUBROUTINE bdotj_par
      
        SUBROUTINE divj
!     
!     WRITTEN 09-10-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes divergence of the current normalized to the total current on the half mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER               :: m, n, istat
        REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
        REAL(dp)              :: tnorm, ton, toff
!-----------------------------------------------
!
!     COMPUTE divJ ON HALF MESH POINTS
!
        CALL second0(ton)

        IF (.NOT. lcurr_init) STOP 'MUST CALL init_state(.TRUE.) BEFORE divJ'

        IF (.NOT. ALLOCATED(divjmnsh))  THEN                                  
          ALLOCATE(divjmnsh(0:mpol,-ntor:ntor,ns),  stat=istat)          ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation failed in DIVJ'
          divjmnsh = zero
        ENDIF
!
!  Compute divergence of J. It is a SINE parity quantity.
!
        DO m = 0, mpol
           DO n = -ntor,ntor
             divjmnsh(m,n,3:nsh) =                                       &
               -m*(ksupumncf(m,n,3:nsh)+ksupumncf(m,n,2:nsh-1))/2        &
               -n*nfp*(ksupvmncf(m,n,3:nsh)+ksupvmncf(m,n,2:nsh-1))/2    &
               + ohs *(ksupsmnsf(m,n,3:nsh)-ksupsmnsf(m,n,2:nsh-1))
           END DO
        END DO

 
        ALLOCATE(temp1(ntheta,nzeta,ns),                                 &
                 temp2(ntheta,nzeta,ns), stat=istat)                     
        IF (istat .NE. 0) STOP 'Allocation failed in DIVJ'
        temp1 = 0; temp2 = 0
 
        temp1 = ksupsijf0*ksubsijf + ksupuijf0*ksubuijf +                & ! Norm: |J|^2
                ksupvijf0*ksubvijf
        CALL to_half_mesh(temp1, temp2, nuv)
        tnorm = SUM(temp2/jacobh**2)
        IF (iam.EQ.0 .AND. ANY(temp2.LT.zero)) PRINT *,'ERROR1 in DIVJ: |J|^2 < 0!'
        tnorm = SQRT(ABS(tnorm))
 
        DEALLOCATE(temp1, temp2, stat=istat)
 !
 !   Compute rms of divergence of J
 !     
        divjmnsh = divjmnsh/tnorm
        divj_rms = SQRT(SUM(divjmnsh**2))

!        IF (iam .EQ. 0) WRITE(4100,*)'DIVJ: tnorm=',tnorm,' divj_rms: ',divj_rms
        CALL second0(toff)
        time_divj = time_divj+(toff-ton)

        END SUBROUTINE divj

        SUBROUTINE divj_par
!     
!     WRITTEN 03-18-13 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes divergence of the current normalized to the total current on the half mesh
!
!#if defined(SKS)
        USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
        USE mpi_inc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER                  :: m, n, istat, nsmin, nsmax, n1, n2
        REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::  temp1, temp2
        REAL(rprec)              :: tnorm, s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE divJ ON HALF MESH POINTS
!
        CALL second0(ton)

        IF (.NOT. lcurr_init) STOP 'MUST CALL init_state_par(.TRUE.) BEFORE divJ_par'

        nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

        IF (.NOT. ALLOCATED(divjmnsh))  THEN                                  
          ALLOCATE(divjmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation failed in DIVJ'
          divjmnsh = 0
        ENDIF

!
!  Compute divergence of J. It is a SINE parity quantity.
!
        n1 = MAX(3, nsmin); n2 = MIN(ns-1, nsmax)
        DO m = 0, mpol
           DO n = -ntor,ntor
             divjmnsh(m,n,n1:n2) =                                       &
               -m*(ksupumncf(m,n,n1:n2)+ksupumncf(m,n,n1-1:n2-1))/2      &
           -n*nfp*(ksupvmncf(m,n,n1:n2)+ksupvmncf(m,n,n1-1:n2-1))/2      &
           + ohs *(ksupsmnsf(m,n,n1:n2)-ksupsmnsf(m,n,n1-1:n2-1))
           END DO
        END DO
!
!    Compute covariant components of current (times jacobian). Used for normalization. 
!  
        nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+1,ns)
        ALLOCATE(temp1(ntheta,nzeta,nsmin:nsmax),                        &
                 temp2(ntheta,nzeta,nsmin:nsmax), stat=istat)                     
        IF (istat .NE. 0) STOP 'Allocation failed in DIVJ'
        temp1 = 0; temp2 = 0

        temp1 = ksupsijf0(:,:,nsmin:nsmax)*ksubsijf(:,:,nsmin:nsmax)     &
              + ksupuijf0(:,:,nsmin:nsmax)*ksubuijf(:,:,nsmin:nsmax)     & ! Norm: |J|^2
              + ksupvijf0(:,:,nsmin:nsmax)*ksubvijf(:,:,nsmin:nsmax)

        CALL to_half_mesh_par(temp1, temp2, nuv, nsmin, nsmax)

        temp2 = temp2/jacobh(:,:,nsmin:nsmax)**2

        nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
        tnorm = SUM(temp2(:,:,nsmin:nsmax))

        DEALLOCATE(temp1, temp2, stat=istat)

        divj_rms = SUM(divjmnsh(:,:,nsmin:nsmax)**2)
#if defined(MPI_OPT)
        s1(1) = divj_rms; s1(2) = tnorm
        CALL MPI_ALLREDUCE(s1,r1,2,MPI_REAL8, MPI_SUM,                   &
                           SIESTA_COMM,MPI_ERR)
        divj_rms=r1(1);   tnorm = r1(2)
#endif   
 !   Compute rms of divergence of J
 !     
        IF (tnorm .GE. zero) THEN
           divjmnsh = divjmnsh/SQRT(tnorm)
           divj_rms = SQRT(divj_rms/tnorm)
        END IF

!        IF (iam .EQ. 0) WRITE(3100,*)'DIVJ_PAR: tnorm=',SQRT(tnorm),' divj_rms: ',divj_rms

        CALL second0(toff)
        time_divj = time_divj+(toff-ton)
!#endif

        END SUBROUTINE divj_par

       SUBROUTINE dealloc_diagnostics

       IF (ALLOCATED(divbmnsf)) DEALLOCATE(divbmnsf)
       IF (ALLOCATED(divjmnsh)) DEALLOCATE(divjmnsh)
       IF (ALLOCATED(bdotjmnch)) DEALLOCATE(bdotjmnch)
       IF (lasym) THEN
       IF (ALLOCATED(bdotjmnsh)) DEALLOCATE(bdotjmnsh)
       END IF

       END SUBROUTINE dealloc_diagnostics


      END MODULE diagnostics_mod
