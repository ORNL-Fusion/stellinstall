!#if defined(SKS)
!*******************************************************************************
!>  \brief Computed harmonics of the MHD force in parallel.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = J X B - Grad(p)                                                 (1)
!>
!>  where J is the plasma current density. In curvilinear coodinates, the forces
!>  components becomes
!>
!>     F_s = (K^u*B^v - K^v*B^u) - dpds                                    (2a)
!>
!>     F_u = (K^v*B^s - K^s*B^v) - dpdu                                    (2b)
!>
!>     F_v = (K^s*B^u - K^u*B^s) - dpdv                                    (2c)
!>
!>  where K^i = sqrt(g)*J^i . Note this routine only works on a
!>  single parity at a time.
!>
!>  \param[in]  pmnf        Fourier pressure. (cos, ipar=0)
!>  \param[in]  jkxbsubsmnc Fourier magnetic pressure and radial pressure
!>                          gradient in the radial direction.
!>  \param[in]  jkxbsubumns Fourier magnetic pressure in the poloidal direction.
!>  \param[in]  jkxbsubvmns Fourier magnetic pressure in the toroidal direction.
!>  \param[out] f_smnf      F_s component of the MHD force. (cos, ipar=0)
!>  \param[out] f_umnf      F_u component of the MHD force. (sin, ipar=0)
!>  \param[out] f_vmnf      F_v component of the MHD force. (sin, ipar=0)
!>  \param[in]  ipar        =0, "even" (default l_asym=F), =1, "odd" 
!*******************************************************************************
       SUBROUTINE get_force_harmonics_par(pmnf,                          &
                                          KxB_smn, KxB_umn, KxB_vmn,     &
                                          f_smnf, f_umnf, f_vmnf, ipar) 
       USE quantities
       USE hessian, ONLY: muPar, l_Compute_Hessian
       USE descriptor_mod, ONLY: iam
       USE evolution, ONLY: l_push_s, l_push_u, l_push_v, l_push_edge, l_Gmres
       USE shared_data, ONLY: gc, nprecon, l_PrintOriginForces, col_scale, fsq_total1
       USE timer_mod, ONLY: get_force_harmonics_time
       USE descriptor_mod, ONLY: iam, nprocs
       USE nscalingtools, ONLY: startglobrow, endglobrow, leftproc,      &
           rightproc, MPI_ERR, FIRSTPARTITION_GE_2, LASTPARTITION_GE_2
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)  ::         &
           pmnf, KxB_smn, KxB_umn, KxB_vmn
       REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(OUT) ::         &
           f_smnf, f_umnf, f_vmnf
       INTEGER, INTENT(IN)                                    ::  ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER  :: m0=0, m1=1, m2=2, moff=1, proc0=0
       REAL(dp), PARAMETER :: p5 = 0.5_dp
       INTEGER             :: m, n, js, js1
       REAL(dp)            :: skston, skstoff 
       REAL(dp)            :: f1, f2, f3, r0, tmps(6), parity
       INTEGER             :: nsmin, nsmax, from
!-----------------------------------------------
!
!      NOTE: there is an m=1 part of pmnf retained at js=1, which is REALLY at js=hs/2
!      It contributes a spurious piece to f_vmnf, which is correctly zeroed
!
       CALL second0(skston)
       IF (ipar .EQ. 0) THEN
          parity = 1
       ELSE
          parity = -1
       END IF
       
       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

       DO js = nsmin,nsmax
          js1 = js-nsmin+1
 	      f_smnf(:,:,js) = KxB_smn(:,:,js1)
          DO n = -ntor, ntor       
             DO m = 0, mpol
                f_umnf(m,n,js) = KxB_umn(m,n,js1)                        &
                               + parity*m*pmnf(m,n,js1)                   ! m -> -(d/du)COS(mu+nfp*nv) for ipar=0
                f_vmnf(m,n,js) = KxB_vmn(m,n,js1)                        &
                               + parity*n*nfp*pmnf(m,n,js1)               ! n*nfp -> -(d/dv)COS(mu+nfp*n*v) for ipar=0
              ENDDO
           ENDDO
       ENDDO

       IF (ANY(f_smnf(m0,-ntor:-1,1:nsmax-nsmin+1) .NE. 0)) STOP 'f_smnf != 0'
       IF (ANY(f_umnf(m0,-ntor:-1,1:nsmax-nsmin+1) .NE. 0)) STOP 'f_umnf != 0'
       IF (ANY(f_vmnf(m0,-ntor:-1,1:nsmax-nsmin+1) .NE. 0)) STOP 'f_vmnf != 0'

!      ADD SYMMETRIC DAMPING TERMS -B v dot B TO SUPPRESS NULL SPACE
!      WORKS BEST: comment out IF statement below, allowing -mu v|| term
!      to be added to BOTH linear and nonlinear forces. However, must
!      guarantee that v->0 at end for this to work well

!CALL funct_island WITH vsupXmn...=0 to get un-augmented MHD force
       IF (muPar.NE.zero .AND. l_Compute_Hessian)                        &
!       IF (muPar.NE.zero)                                               &
           CALL AddParDamping_par(f_smnf, f_umnf, f_vmnf, KxB_vmn,       &
                                  nsmin, nsmax, ipar)

!GATHER BDY FORCES AND PRINT OUT ON PROC=0
       PRINT_O: IF (l_PrintOriginForces) THEN
             fbdy(1) = SQRT(SUM(f_smnf(m1,:,1)**2))
             fbdy(2) = SQRT(SUM(f_smnf(m1,:,2)**2))
             fbdy(3) = SQRT(SUM(f_smnf(m0,:,2)**2) + SUM(f_smnf(m2:,:,2)**2))

             fbdy(4) = SQRT(SUM(f_umnf(m1,:,1)**2))
             fbdy(5) = SQRT(SUM(f_umnf(m1,:,2)**2))
             fbdy(6) = SQRT(SUM(f_umnf(m0,:,2)**2) + SUM(f_umnf(m2:,:,2)**2))

             fbdy(7) = SQRT(SUM(f_vmnf(m0,:,1)**2))
             fbdy(8) = SQRT(SUM(f_vmnf(m0,:,2)**2))
             fbdy(9) = SQRT(SUM(f_vmnf(m1:,:,2)**2))

             fbdy(10) = SQRT(SUM(f_umnf(:,:,ns)**2))
             fbdy(11) = SQRT(SUM(f_umnf(:,:,ns-1)**2))
             fbdy(12) = SQRT(SUM(f_vmnf(:,:,ns)**2))
             fbdy(13) = SQRT(SUM(f_vmnf(:,:,ns-1)**2))

#if defined(MPI_OPT)
           NPROCS_0: IF (nprocs.GT.1) THEN
             IF (.NOT.FIRSTPARTITION_GE_2) THEN
                tmps(1)=fbdy(2); tmps(2)=fbdy(3)
                tmps(3)=fbdy(5); tmps(4)=fbdy(6)
                tmps(5)=fbdy(8); tmps(6)=fbdy(9)
                from=2
                CALL MPI_BCast(tmps,6,MPI_REAL8,from-1,MPI_COMM_WORLD,MPI_ERR)
                fbdy(2)=tmps(1); fbdy(3)=tmps(2)
                fbdy(5)=tmps(3); fbdy(6)=tmps(4)
                fbdy(8)=tmps(5); fbdy(9)=tmps(6)
             END IF

             IF (.NOT.LASTPARTITION_GE_2) THEN
                tmps(1)=fbdy(11); tmps(2)=fbdy(13)
                from=nprocs-1
                CALL MPI_BCast(tmps,2,MPI_REAL8,from-1,MPI_COMM_WORLD,MPI_ERR)
                fbdy(11)=tmps(1); fbdy(13)=tmps(2)
             END IF

             tmps(1)=fbdy(10); tmps(2)=fbdy(11)
             tmps(3)=fbdy(12); tmps(4)=fbdy(13)
             from=nprocs
             CALL MPI_BCast(tmps,4,MPI_REAL8,from-1,MPI_COMM_WORLD,MPI_ERR)
             fbdy(10)=tmps(1); fbdy(11)=tmps(2)
             fbdy(12)=tmps(3); fbdy(13)=tmps(4)
           END IF NPROCS_0
#endif
           IF (iam .EQ. 0) THEN
              WRITE(*,121) fbdy(1),fbdy(2),fbdy(3)
              WRITE(*,122) fbdy(4),fbdy(5),fbdy(6)
              WRITE(*,123) fbdy(7),fbdy(8),fbdy(9)
              WRITE(*,124) fbdy(10),fbdy(11)
              WRITE(*,125) fbdy(12),fbdy(13)
           END IF

       END IF PRINT_O

 121   FORMAT(/,' fs(1,m=1): ',1p,e12.3,' fs(2,m=1): ',1pe12.3,' fs(2,m!=1):',1pe12.3)
 122   FORMAT(  ' fu(1,m=1): ',1p,e12.3,' fu(2,m=1): ',1pe12.3,' fu(2,m!=1):',1pe12.3)
 123   FORMAT(  ' fv(1,m=0): ',1p,e12.3,' fv(2,m=0): ',1pe12.3,' fv(2,m>0): ',1pe12.3)
 124   FORMAT(  ' fu(ns)   : ',1p,e12.3,' fu(ns-1):  ',1pe12.3)
 125   FORMAT(  ' fv(ns)   : ',1p,e12.3,' fv(ns-1):  ',1pe12.3,/)

!      ORIGIN BOUNDARY CONDITIONS
       r0 = hs_i/2
       NSMIN_1: IF (nsmin .EQ. 1) THEN
       f_smnf(m0,:,1) = 0; f_smnf(m2:,:,1) = 0
!       IF (.not.l_natural) f_umnf(m1,:,1) = f_umnf(m1,:,1)-r0*f_smnf(m1,:,1) 
       IF (.NOT.l_push_s) THEN 
          f_smnf(m1,:,1) = 0
       ELSE
          f_smnf(m1,:,1) = p5*f_smnf(m1,:,1)
       END IF

       f_umnf(m0,:,1) = 0; f_umnf(m2:,:,1) = 0
       IF (.NOT.l_push_u) THEN
          f_umnf(m1,:,1) = 0
       ELSE
          f_umnf(m1,:,1) = p5*f_umnf(m1,:,1)
       END IF

       f_vmnf(m1:,:,1) = 0
       IF (.NOT.l_push_v) THEN
          f_vmnf(m0,:,1) = 0
       ELSE
          f_vmnf(m0,:,1) = p5*f_vmnf(m0,:,1)   !NEED THIS TO MAINTAIN SYMMETRY OF MATRIX!
       END IF

       END IF NSMIN_1
 
!      EDGE BOUNDARY CONDITIONS
       NSMAX_NS: IF (nsmax .EQ. ns) THEN
       f_smnf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             f_umnf(:,:,ns) = 0
             f_vmnf(:,:,ns) = 0
          ELSE
             f_umnf(:,:,ns) = p5*f_umnf(:,:,ns)
             f_vmnf(:,:,ns) = p5*f_vmnf(:,:,ns)
          END IF
       END IF NSMAX_NS

       m = 1
       IF (ipar .EQ. 1) m = 4
       f_smnf(:,:,nsmin:nsmax) = signjac*f_smnf(:,:,nsmin:nsmax)*col_scale(:,:,m,nsmin:nsmax)
       f_umnf(:,:,nsmin:nsmax) = signjac*f_umnf(:,:,nsmin:nsmax)*col_scale(:,:,m+1,nsmin:nsmax)
       f_vmnf(:,:,nsmin:nsmax) = signjac*f_vmnf(:,:,nsmin:nsmax)*col_scale(:,:,m+2,nsmin:nsmax)

       CALL second0(skstoff)
       get_force_harmonics_time=get_force_harmonics_time+(skstoff-skston)

       END SUBROUTINE get_force_harmonics_par


       SUBROUTINE addpardamping_par (f_smnf, f_umnf, f_vmnf, work,       &
                                     nsmin, nsmax, ipar)
       USE quantities
       USE hessian, ONLY: muPar, mupar_norm, l_Compute_Hessian
       USE fourier, ONLY: tomnsp_par, f_sin, f_cos, f_full
       USE shared_data, ONLY:  col_scale
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: nsmin, nsmax, ipar
       REAL(dp)            :: work(0:mpol,-ntor:ntor,nsmin:nsmax)
       REAL(dp), INTENT(INOUT), DIMENSION(0:mpol,-ntor:ntor,ns) ::       &
                 f_smnf, f_umnf, f_vmnf          
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: js, istat, ns_span, m
       REAL(dp) :: av(ns)
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: vdotB, tmp
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: xsubsij, xsubuij, xsubvij
!-----------------------------------------------
!
!      Compute parallel (to B) damping force and add it to (linearized) forces
!
       ns_span = nsmax-nsmin+1
       m = 1
       IF (ipar .EQ. 1) m = 4

       ALLOCATE (vdotB(ntheta,nzeta,nsmin:nsmax),                        &
                 tmp(ntheta,nzeta,nsmin:nsmax), stat=istat)

!SPH 05-6-16 UPDATE mupar_norm IN cv_currents and put in 1/jacobh factor
!SIMPLER (FASTER, APPROXIMATE) FORM
!
!       GOTO 200

       tmp = bsubsijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
       CALL SurfAverage(av(nsmin:nsmax), tmp, nsmin, nsmax)

       DO js = nsmin,nsmax
          av(js) = av(js)*muPar*mupar_norm(js)
          IF (ipar .EQ. 0) THEN
             work(:,:,js) = jvsupsmncf(:,:,js)
          ELSE
             work(:,:,js) = jvsupsmnsf(:,:,js)
          END IF
          f_smnf(:,:,js) = f_smnf(:,:,js)                                &
                         + work(:,:,js)*av(js)*col_scale(:,:,m,js)
       END DO   

       tmp = bsubuijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
       CALL SurfAverage(av(nsmin:nsmax), tmp, nsmin, nsmax)

       DO js = nsmin,nsmax
          av(js) = av(js)*muPar*mupar_norm(js)
          IF (ipar .EQ. 0) THEN
             work(:,:,js) = jvsupumnsf(:,:,js)
          ELSE
             work(:,:,js) = jvsupumncf(:,:,js)
          END IF
          f_umnf(:,:,js) = f_umnf(:,:,js)                                &
                         + work(:,:,js)*av(js)*col_scale(:,:,m+1,js) 
       END DO   

       tmp = bsubvijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
       CALL SurfAverage(av(nsmin:nsmax), tmp, nsmin, nsmax)

       DO js = nsmin,nsmax
          av(js) = av(js)*muPar*mupar_norm(js)
          IF (ipar .EQ. 0) THEN
             work(:,:,js) = jvsupvmnsf(:,:,js)
          ELSE
             work(:,:,js) = jvsupvmncf(:,:,js)
          END IF
          f_vmnf(:,:,js) = f_vmnf(:,:,js)                                &
                         + work(:,:,js)*av(js)*col_scale(:,:,m+2,js)
       END DO   

       GO TO 100       

 200   CONTINUE

       IF (nsmin .EQ. 1) bsubuijf(:,:,1) = 0
       vdotB(:,:,nsmin:nsmax) =                                          &
           jvsupsijf(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)          &
         + jvsupuijf(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)          &
         + jvsupvijf(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)

       DO js = nsmin,nsmax
          vdotB(:,:,js) = muPar*mupar_norm(js)*vdotB(:,:,js)/jacobf(:,:,js)
       END DO

       tmp = bsubsijf(:,:,nsmin:nsmax)*vdotB
       m = f_cos
       IF (ipar .EQ. 1) m = f_sin
       CALL tomnsp_par(tmp, work, m, f_full, 1, ns_span)
       DO js=nsmin,nsmax
          f_smnf(:,:,js) = f_smnf(:,:,js) + work(:,:,js)
       END DO

       m = f_sin
       IF (ipar .EQ. 1) m = f_cos
       tmp = bsubuijf(:,:,nsmin:nsmax)*vdotB
       CALL tomnsp_par(tmp, work, m, f_full, 1, ns_span)
       DO js=nsmin,nsmax
          f_umnf(:,:,js) = f_umnf(:,:,js) + work(:,:,js)
       END DO

       tmp = bsubvijf(:,:,nsmin:nsmax)*vdotB
       CALL tomnsp_par(tmp, work, m, f_full, 1, ns_span)
       DO js=nsmin,nsmax
          f_vmnf(:,:,js) = f_vmnf(:,:,js) + work(:,:,js)
       END DO

 100   CONTINUE
       DEALLOCATE (vdotB, tmp)

       END SUBROUTINE addpardamping_par
!#endif

       SUBROUTINE get_force_harmonics(pmnf,                              &
                                      KxB_smn, KxB_umn, KxB_vmn,         &
                                      f_smnf, f_umnf, f_vmnf, ipar) 
       USE quantities,ONLY: dp, mpol, ntor, ns
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)  ::         &
           pmnf,  KxB_smn, KxB_umn, KxB_vmn
       REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(OUT)  ::        &
           f_smnf, f_umnf, f_vmnf
       INTEGER, INTENT(IN) :: ipar
!-----------------------------------------------
       CALL get_force_harmonics_par(pmnf,                                &
                                    KxB_smn, KxB_umn, KxB_vmn,           &
                                    f_smnf, f_umnf, f_vmnf, ipar) 
       
       END SUBROUTINE get_force_harmonics

