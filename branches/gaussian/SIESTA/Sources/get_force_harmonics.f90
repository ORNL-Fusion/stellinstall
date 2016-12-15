!*******************************************************************************
!>  \file get_force_harmonics.f90
!>  \brief Computed harmonics of the MHD force.
!*******************************************************************************

#if defined(SKS)
!*******************************************************************************
!>  \brief Computed harmonics of the MHD force in parallel.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = K X B - Grad(p)                                                   (1)
!>
!>  where K is the plasma current density. In curalinear coodinates, the forces
!>  components becomes
!>
!>     F_s = J(K^u*B^v - K^v*B^u) - dpds                                    (2a)
!>
!>     F_u = J(K^v*B^s - K^s*B^v) - dpdu                                    (2b)
!>
!>     F_v = J(K^s*B^u - K^u*B^s) - dpdv                                    (2c)
!>
!>  where J is the coordinate system Jacobian. Note this routine only works on a
!>  single parity.
!>
!>  \param[in]  pmncf       Fourier pressure.
!>  \param[in]  jkxbsubsmnc Fourier magnetic pressure and radial pressure
!>                          gradient in the radial direction.
!>  \param[in]  jkxbsubumns Fourier magnetic pressure in the poloidal direction.
!>  \param[in]  jkxbsubvmns Fourier magnetic pressure in the toroidal direction.
!>  \param[out] fsubsmncf   F_s component of the MHD force.
!>  \param[out] fsubumnsf   F_u component of the MHD force.
!>  \param[out] fsubvmnsf   F_v component of the MHD force.
!*******************************************************************************
       SUBROUTINE get_force_harmonics_par(pmncf,                                 &
      &                                   jkxbsubsmnc, jkxbsubumns, jkxbsubvmns, &
      &                                   fsubsmncf, fsubumnsf, fsubvmnsf)
       USE quantities, dum1=>fsubsmncf, dum2=>fsubumnsf, dum3=>fsubvmnsf
       USE hessian, ONLY: muPar, l_Compute_Hessian
       USE descriptor_mod, ONLY: iam
       USE evolution, ONLY: l_push_s, l_push_u, l_push_v, l_push_edge
       USE shared_data, ONLY: gc, nprecon, l_PrintOriginForces, col_scale
       USE timer_mod, ONLY: get_force_harmonics_time
       USE descriptor_mod, ONLY: iam, nprocs
       USE nscalingtools, ONLY: startglobrow, endglobrow, leftproc,     &
           rightproc, MPI_ERR, FIRSTPARTITION_GE_2, LASTPARTITION_GE_2
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)  ::     &
           pmncf, jkxbsubsmnc, jkxbsubumns, jkxbsubvmns
       REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(OUT)  ::    &
           fsubsmncf, fsubumnsf, fsubvmnsf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER     :: m0=0, m1=1, m2=2, moff=1, proc0=0
       REAL(rprec), PARAMETER :: p5 = 0.5_dp
       INTEGER                :: m, n, js, js1
       REAL(rprec)            :: skston, skstoff
       REAL(rprec)            :: f1, f2, f3, r0, tmps(6)
       INTEGER                :: nsmin, nsmax, from
!-----------------------------------------------
!
!      NOTE: there is an m=1 part of pmncf retained at js=1, which is REALLY at js=hs/2
!      It contributes a spurious piece to fsubvmnsf, which is correctly zeroed
!
       CALL second0(skston)
       gc = 0                                                            ! zero all forces
       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

       DO js = nsmin,nsmax
          js1 = js-nsmin+1
 	      fsubsmncf(:,:,js) = jkxbsubsmnc(:,:,js1) ! J(K^u*B^v - K^v*B^u) - dpds Note already contains the dpds term.
          DO n = -ntor, ntor       
             DO m = 0, mpol
                fsubumnsf(m,n,js) = jkxbsubumns(m,n,js1) + m*pmncf(m,n,js1)    ! J(K^v*B^s - K^s*B^v) - dpdu  m -> -(d/du)COS(mu+nv)
                fsubvmnsf(m,n,js) = jkxbsubvmns(m,n,js1) + n*nfp*pmncf(m,n,js1)! J(K^s*B^u - K^u*B^s) - dpdv  n*nfp -> -(d/dv)COS(mu+n*NFP*v)
             ENDDO
          ENDDO
          fsubsmncf(m0,-ntor:-1,js) = 0                                  !redundant COS modes
          fsubumnsf(m0,-ntor:0 ,js) = 0                                  !redundant SIN modes
          fsubvmnsf(m0,-ntor:0 ,js) = 0                                  !"                 "
       ENDDO

!      ADD SYMMETRIC DAMPING TERMS -B v dot B TO SUPPRESS NULL SPACE
!      WORKS BEST: comment out IF statement below, allowing -mu v|| term
!      to be added to BOTH linear and nonlinear forces. However, must
!      guarantee that v->0 at end for this to work well

!CALL funct_island WITH vsupXmn...=0 to get un-augmented MHD force
       IF (l_Compute_Hessian .AND. muPar.NE.zero) THEN                   !Always (lin/nonlin) use mu|| = 0 EXCEPT for Hessian calculation
          CALL AddParDamping_par(fsubsmncf, fsubumnsf, fsubvmnsf, jkxbsubvmns, nsmin, nsmax) ! jkxbsubvmns used as a work buffer.
       END IF

!GATHER BDY FORCES AND PRINT OUT ON PROC=0
       PRINT_O: IF (l_PrintOriginForces) THEN
             fbdy(1) = SQRT(SUM(fsubsmncf(m1,:,1)**2))
             fbdy(2) = SQRT(SUM(fsubsmncf(m1,:,2)**2))
             fbdy(3) = SQRT(SUM(fsubsmncf(m0,:,2)**2) + SUM(fsubsmncf(m2:,:,2)**2))

             fbdy(4) = SQRT(SUM(fsubumnsf(m1,:,1)**2))
             fbdy(5) = SQRT(SUM(fsubumnsf(m1,:,2)**2))
             fbdy(6) = SQRT(SUM(fsubumnsf(m0,:,2)**2) + SUM(fsubumnsf(m2:,:,2)**2))

             fbdy(7) = SQRT(SUM(fsubvmnsf(m0,:,1)**2))
             fbdy(8) = SQRT(SUM(fsubvmnsf(m0,:,2)**2))
             fbdy(9) = SQRT(SUM(fsubvmnsf(m1:,:,2)**2))

             fbdy(10) = SQRT(SUM(fsubumnsf(:,:,ns)**2))
             fbdy(11) = SQRT(SUM(fsubumnsf(:,:,ns-1)**2))
             fbdy(12) = SQRT(SUM(fsubvmnsf(:,:,ns)**2))
             fbdy(13) = SQRT(SUM(fsubvmnsf(:,:,ns-1)**2))

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
       NSMIN_1: IF (nsmin .EQ. 1) THEN
       fsubsmncf(m0,:,1) = 0; fsubsmncf(m2:,:,1) = 0
!       r0 = hs_i/2
!       IF (.not.l_natural) fsubumnsf(m1,:,1) = fsubumnsf(m1,:,1)-r0*fsubsmncf(m1,:,1) 
       IF (.not.l_push_s) THEN 
          fsubsmncf(m1,:,1) = 0
       ELSE
          fsubsmncf(m1,:,1) = p5*fsubsmncf(m1,:,1)
       END IF

       fsubumnsf(m0,:,1) = 0; fsubumnsf(m2:,:,1) = 0
       IF (.not.l_push_u) THEN
          fsubumnsf(m1,:,1) = 0
       ELSE
          fsubumnsf(m1,:,1) = p5*fsubumnsf(m1,:,1)
       END IF

       fsubvmnsf(m1:,:,1) = 0
       IF (.not.l_push_v) THEN
          fsubvmnsf(m0,:,1) = 0
       ELSE
          fsubvmnsf(m0,:,1) = p5*fsubvmnsf(m0,:,1)
       END IF

       END IF NSMIN_1
 
!      EDGE BOUNDARY CONDITIONS (factor of p5 due to one-sided average at edge in var principle)
       NSMAX_NS: IF (nsmax .EQ. ns) THEN
       fsubsmncf(:,:,ns) = 0
       IF (.NOT.l_push_edge) THEN
          fsubumnsf(:,:,ns) = 0
          fsubvmnsf(:,:,ns) = 0
       ELSE
          fsubumnsf(:,:,ns) = p5*fsubumnsf(:,:,ns)
          fsubvmnsf(:,:,ns) = p5*fsubvmnsf(:,:,ns)
       END IF
       END IF NSMAX_NS


       fsubsmncf(:,:,nsmin:nsmax) = signjac*fsubsmncf(:,:,nsmin:nsmax)*col_scale(:,:,1,nsmin:nsmax)
       fsubumnsf(:,:,nsmin:nsmax) = signjac*fsubumnsf(:,:,nsmin:nsmax)*col_scale(:,:,2,nsmin:nsmax)
       fsubvmnsf(:,:,nsmin:nsmax) = signjac*fsubvmnsf(:,:,nsmin:nsmax)*col_scale(:,:,3,nsmin:nsmax)

       CALL second0(skstoff)
       get_force_harmonics_time=get_force_harmonics_time+(skstoff-skston)

       END SUBROUTINE get_force_harmonics_par

#if defined(MRC)
!*******************************************************************************
!>  \brief Computed harmonics of the MHD force in parallel for other parity.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = K X B - Grad(p)                                                   (1)
!>
!>  where K is the plasma current density. In curalinear coodinates, the forces
!>  components becomes
!>
!>     F_s = J(K^u*B^v - K^v*B^u) - dpds                                    (2a)
!>
!>     F_u = J(K^v*B^s - K^s*B^v) - dpdu                                    (2b)
!>
!>     F_v = J(K^s*B^u - K^u*B^s) - dpdv                                    (2c)
!>
!>  where J is the coordinate system Jacobian. Note this routine only works on a
!>  single parity.
!>
!>  \param[in]  pmnsf       Fourier pressure.
!>  \param[in]  jkxbsubsmns Fourier magnetic pressure and radial pressure
!>                          gradient in the radial direction.
!>  \param[in]  jkxbsubumnc Fourier magnetic pressure in the poloidal direction.
!>  \param[in]  jkxbsubvmnc Fourier magnetic pressure in the toroidal direction.
!>  \param[out] fsubsmnsf   F_s component of the MHD force.
!>  \param[out] fsubumncf   F_u component of the MHD force.
!>  \param[out] fsubvmncf   F_v component of the MHD force.
!*******************************************************************************
       SUBROUTINE get_force_harmonics_par_other(pmnsf,                           &
      &                                   jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc, &
      &                                   fsubsmnsf, fsubumncf, fsubvmncf)
       USE quantities, dum1=>fsubsmnsf, dum2=>fsubumncf, dum3=>fsubvmncf
       USE hessian, ONLY: muPar, l_Compute_Hessian
       USE descriptor_mod, ONLY: iam
       USE evolution, ONLY: l_push_s, l_push_u, l_push_v, l_push_edge
       USE shared_data, ONLY: gc, nprecon, l_PrintOriginForces, col_scale
       USE timer_mod, ONLY: get_force_harmonics_time
       USE descriptor_mod, ONLY: iam, nprocs
       USE nscalingtools, ONLY: startglobrow, endglobrow, leftproc,     &
           rightproc, MPI_ERR, FIRSTPARTITION_GE_2, LASTPARTITION_GE_2
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)  ::     &
           pmnsf, jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc
       REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(OUT)  ::    &
           fsubsmnsf, fsubumncf, fsubvmncf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER     :: m0=0, m1=1, m2=2, moff=1, proc0=0
       REAL(rprec), PARAMETER :: p5 = 0.5_dp
       INTEGER                :: m, n, js, js1
       REAL(rprec)            :: skston, skstoff
       REAL(rprec)            :: f1, f2, f3, r0, tmps(6)
       INTEGER                :: nsmin, nsmax, from
!-----------------------------------------------
!
!      NOTE: there is an m=1 part of pmncf retained at js=1, which is REALLY at js=hs/2
!      It contributes a spurious piece to fsubvmnsf, which is correctly zeroed
!
       CALL second0(skston)
       gc = 0                                                            ! zero all forces
       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

       DO js = nsmin,nsmax
          js1 = js-nsmin+1
          fsubsmnsf(:,:,js) = jkxbsubsmns(:,:,js1) ! J(K^u*B^v - K^v*B^u) - dpds Note already contains the dpds term.
          DO n = -ntor, ntor
             DO m = 0, mpol
                fsubumncf(m,n,js) = jkxbsubumnc(m,n,js1) - m*pmnsf(m,n,js1)    ! J(K^v*B^s - K^s*B^v) - dpdu  m -> (d/du)SIN(mu+nv)
                fsubvmncf(m,n,js) = jkxbsubvmnc(m,n,js1) - n*nfp*pmnsf(m,n,js1)! J(K^s*B^u - K^u*B^s) - dpdv  n*nfp -> (d/dv)SIN(mu+n*NFP*v)
             ENDDO
          ENDDO
          fsubsmnsf(m0,-ntor:0,js) = 0                                    !redundant SIN modes
          fsubumncf(m0,-ntor:-1 ,js) = 0                                  !redundant COS modes
          fsubvmncf(m0,-ntor:-1 ,js) = 0                                  !"                 "
       ENDDO

!      ADD SYMMETRIC DAMPING TERMS -B v dot B TO SUPPRESS NULL SPACE
!      WORKS BEST: comment out IF statement below, allowing -mu v|| term
!      to be added to BOTH linear and nonlinear forces. However, must
!      guarantee that v->0 at end for this to work well

!CALL funct_island WITH vsupXmn...=0 to get un-augmented MHD force
       IF (l_Compute_Hessian .AND. muPar.NE.zero) THEN                   !Always (lin/nonlin) use mu|| = 0 EXCEPT for Hessian calculation
          CALL AddParDamping_par_other(fsubsmnsf, fsubumncf, fsubvmncf,        &
                                       jkxbsubvmnc, nsmin, nsmax) ! jkxbsubvmnc used as a work buffer.
       END IF

!GATHER BDY FORCES AND PRINT OUT ON PROC=0
       PRINT_O: IF (l_PrintOriginForces) THEN
          fbdy(1) = SQRT(SUM(fsubsmnsf(m1,:,1)**2))
          fbdy(2) = SQRT(SUM(fsubsmnsf(m1,:,2)**2))
          fbdy(3) = SQRT(SUM(fsubsmnsf(m0,:,2)**2) + SUM(fsubsmnsf(m2:,:,2)**2))

          fbdy(4) = SQRT(SUM(fsubumncf(m1,:,1)**2))
          fbdy(5) = SQRT(SUM(fsubumncf(m1,:,2)**2))
          fbdy(6) = SQRT(SUM(fsubumncf(m0,:,2)**2) + SUM(fsubumncf(m2:,:,2)**2))

          fbdy(7) = SQRT(SUM(fsubvmncf(m0,:,1)**2))
          fbdy(8) = SQRT(SUM(fsubvmncf(m0,:,2)**2))
          fbdy(9) = SQRT(SUM(fsubvmncf(m1:,:,2)**2))

          fbdy(10) = SQRT(SUM(fsubumncf(:,:,ns)**2))
          fbdy(11) = SQRT(SUM(fsubumncf(:,:,ns-1)**2))
          fbdy(12) = SQRT(SUM(fsubvmncf(:,:,ns)**2))
          fbdy(13) = SQRT(SUM(fsubvmncf(:,:,ns-1)**2))

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
       NSMIN_1: IF (nsmin .EQ. 1) THEN
          fsubsmnsf(m0,:,1) = 0; fsubsmnsf(m2:,:,1) = 0
!       r0 = hs_i/2
!       IF (.not.l_natural) fsubumncf(m1,:,1) = fsubumncf(m1,:,1)-r0*fsubsmnsf(m1,:,1)
          IF (.not.l_push_s) THEN
             fsubsmnsf(m1,:,1) = 0
          ELSE
             fsubsmnsf(m1,:,1) = p5*fsubsmnsf(m1,:,1)
          END IF

          fsubumncf(m0,:,1) = 0; fsubumncf(m2:,:,1) = 0
          IF (.not.l_push_u) THEN
             fsubumncf(m1,:,1) = 0
          ELSE
             fsubumncf(m1,:,1) = p5*fsubumncf(m1,:,1)
          END IF

          fsubvmncf(m1:,:,1) = 0
          IF (.not.l_push_v) THEN
             fsubvmncf(m0,:,1) = 0
          ELSE
             fsubvmncf(m0,:,1) = p5*fsubvmncf(m0,:,1)
          END IF

       END IF NSMIN_1

!      EDGE BOUNDARY CONDITIONS (factor of p5 due to one-sided average at edge in var principle)
       NSMAX_NS: IF (nsmax .EQ. ns) THEN
          fsubsmnsf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
              fsubumncf(:,:,ns) = 0
              fsubvmncf(:,:,ns) = 0
          ELSE
              fsubumncf(:,:,ns) = p5*fsubumncf(:,:,ns)
              fsubvmncf(:,:,ns) = p5*fsubvmncf(:,:,ns)
          END IF
       END IF NSMAX_NS


       fsubsmnsf(:,:,nsmin:nsmax) = signjac*fsubsmnsf(:,:,nsmin:nsmax)*col_scale(:,:,1,nsmin:nsmax)
       fsubumncf(:,:,nsmin:nsmax) = signjac*fsubumncf(:,:,nsmin:nsmax)*col_scale(:,:,2,nsmin:nsmax)
       fsubvmncf(:,:,nsmin:nsmax) = signjac*fsubvmncf(:,:,nsmin:nsmax)*col_scale(:,:,3,nsmin:nsmax)

       CALL second0(skstoff)
       get_force_harmonics_time=get_force_harmonics_time+(skstoff-skston)

       END SUBROUTINE get_force_harmonics_par_other
#endif

!*******************************************************************************
!>  \brief Add parallel damping force in parallel.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = -B*v.B*mu_||*mu_norm                                              (1)
!>
!>  where mu is the parallel
!>
!>     F_s = B_s*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>     F_u = B_u*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>     F_v = B_v*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>  where J is the coordinate system Jacobian. Note this routine only works on a
!>  single parity.
!>
!>  \param[inout] fsubsmncf F_s component of the MHD force.
!>  \param[inout] fsubumnsf F_u component of the MHD force.
!>  \param[inout] fsubvmnsf F_v component of the MHD force.
!>  \param        work      Temporary work buffer.
!>  \param[in]    nsmin     Min s index.
!>  \param[in]    nsmax     Max s index.
!*******************************************************************************
       SUBROUTINE addpardamping_par (fsubsmncf, fsubumnsf, fsubvmnsf, work, nsmin, nsmax)
       USE quantities, dum1=>fsubsmncf, dum2=>fsubumnsf, dum3=>fsubvmnsf
       USE hessian, ONLY: muPar, mupar_norm
       USE fourier, ONLY: tomnsp_par
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: nsmin, nsmax
       REAL(dp)            :: work(0:mpol,-ntor:ntor,nsmin:nsmax)
       REAL(dp), INTENT(INOUT), DIMENSION(0:mpol,-ntor:ntor,ns) ::      &
                 fsubsmncf, fsubumnsf, fsubvmnsf          
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: js, istat, ns_span
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: vdotB, tmp
!-----------------------------------------------
!
!      Compute parallel (to B) damping force and add it to (linearized) forces
!
       ns_span = nsmax-nsmin+1
       ALLOCATE (vdotB(ntheta,nzeta,nsmin:nsmax),                       &
                 tmp(ntheta,nzeta,nsmin:nsmax), stat=istat)

       IF (istat .NE. 0) STOP 'ALLOCATION ERROR IN ADDPARDAMPING'
!
!      compute (sin-parity) v * B: note the -signJac factor is in mupar_norm
!
       IF (nsmin .EQ. 1) bsubuijf(:,:,1) = 0
       vdotB(:,:,nsmin:nsmax) =                                         &
           jvsupsijf(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)         &
         + jvsupuijf(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)         &
         + jvsupvijf(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)

       DO js = nsmin,nsmax
          vdotB(:,:,js) = muPar*mupar_norm(js)*vdotB(:,:,js)
       END DO

       tmp = bsubsijf(:,:,nsmin:nsmax)*vdotB                             !COS SYM
       CALL tomnsp_par(tmp, work, 0, 0, 1, ns_span)
       DO js=nsmin,nsmax
          fsubsmncf(:,:,js) = fsubsmncf(:,:,js) + work(:,:,js)
       END DO

       tmp = bsubuijf(:,:,nsmin:nsmax)*vdotB                             !SIN SYM
       CALL tomnsp_par(tmp, work, 1, 0, 1, ns_span)
       DO js=nsmin,nsmax
          fsubumnsf(:,:,js) = fsubumnsf(:,:,js) + work(:,:,js)
       END DO

       tmp = bsubvijf(:,:,nsmin:nsmax)*vdotB                             !SIN SYM
       CALL tomnsp_par(tmp, work, 1, 0, 1, ns_span)
       DO js=nsmin,nsmax
          fsubvmnsf(:,:,js) = fsubvmnsf(:,:,js) + work(:,:,js)
       END DO

       DEALLOCATE (vdotB, tmp)

       END SUBROUTINE addpardamping_par

#if defined(MRC)
!*******************************************************************************
!>  \brief Add parallel damping force in parallel for other parity.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = -B*v.B*mu_||*mu_norm                                              (1)
!>
!>  where mu is the parallel
!>
!>     F_s = B_s*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>     F_u = B_u*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>     F_v = B_v*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>  where J is the coordinate system Jacobian. Note this routine only works on a
!>  single parity.
!>
!>  \param[inout] fsubsmnsf F_s component of the MHD force.
!>  \param[inout] fsubumncf F_u component of the MHD force.
!>  \param[inout] fsubvmncf F_v component of the MHD force.
!>  \param        work      Temporary work buffer.
!>  \param[in]    nsmin     Min s index.
!>  \param[in]    nsmax     Max s index.
!*******************************************************************************
       SUBROUTINE AddParDamping_par_other(fsubsmnsf, fsubumncf, fsubvmncf, work, nsmin, nsmax)
       USE quantities, dum1=>fsubsmnsf, dum2=>fsubumncf, dum3=>fsubvmncf
       USE hessian, ONLY: muPar, mupar_norm
       USE fourier, ONLY: tomnsp_par, fourier_parity_odd, fourier_parity_even, &
                          fourier_full_grid
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: nsmin, nsmax
       REAL(dp)            :: work(0:mpol,-ntor:ntor,nsmin:nsmax)
       REAL(dp), INTENT(INOUT), DIMENSION(0:mpol,-ntor:ntor,ns) ::      &
                 fsubsmnsf, fsubumncf, fsubvmncf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: js, istat, ns_span
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: vdotB, tmp
!-----------------------------------------------
!
!      Compute parallel (to B) damping force and add it to (linearized) forces
!
       ns_span = nsmax-nsmin+1
       ALLOCATE (vdotB(ntheta,nzeta,nsmin:nsmax),                       &
       tmp(ntheta,nzeta,nsmin:nsmax), stat=istat)

       IF (istat .NE. 0) STOP 'ALLOCATION ERROR IN ADDPARDAMPING'
!
!      compute (sin-parity) v * B: note the -signJac factor is in mupar_norm
!
       IF (nsmin .EQ. 1) bsubuijf(:,:,1) = 0
       vdotB(:,:,nsmin:nsmax) =                                         &
           jvsupsijf(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)         &
         + jvsupuijf(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)         &
         + jvsupvijf(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)

       DO js = nsmin,nsmax
          vdotB(:,:,js) = muPar*mupar_norm(js)*vdotB(:,:,js)
       END DO

       tmp = bsubsijf(:,:,nsmin:nsmax)*vdotB                             !COS SYM
       CALL tomnsp_par(tmp, work, fourier_parity_odd, fourier_full_grid, 1, ns_span)
       DO js=nsmin,nsmax
          fsubsmnsf(:,:,js) = fsubsmnsf(:,:,js) + work(:,:,js)
       END DO

       tmp = bsubuijf(:,:,nsmin:nsmax)*vdotB                             !SIN SYM
       CALL tomnsp_par(tmp, work, fourier_parity_even, fourier_full_grid, 1, ns_span)
       DO js=nsmin,nsmax
          fsubumncf(:,:,js) = fsubumncf(:,:,js) + work(:,:,js)
       END DO

       tmp = bsubvijf(:,:,nsmin:nsmax)*vdotB                             !SIN SYM
       CALL tomnsp_par(tmp, work, fourier_parity_even, fourier_full_grid, 1, ns_span)
       DO js=nsmin,nsmax
          fsubvmncf(:,:,js) = fsubvmncf(:,:,js) + work(:,:,js)
       END DO

       DEALLOCATE (vdotB, tmp)

       END SUBROUTINE AddParDamping_par_other
#endif
#endif

!*******************************************************************************
!>  \brief Computed harmonics of the MHD force serially.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = K X B - Grad(p)                                                   (1)
!>
!>  where K is the plasma current density. In curalinear coodinates, the forces
!>  components becomes
!>
!>     F_s = J(K^u*B^v - K^v*B^u) - dpds                                    (2a)
!>
!>     F_u = J(K^v*B^s - K^s*B^v) - dpdu                                    (2b)
!>
!>     F_v = J(K^s*B^u - K^u*B^s) - dpdv                                    (2c)
!>
!>  where J is the coordinate system Jacobian. Note this routine only works on a
!>  single parity.
!>
!>  \param[in]  pmncf       Fourier pressure.
!>  \param[in]  pmnc_ds     Radial derivative of fourier pressure.
!>  \param[in]  jkxbsubsmnc Fourier magnetic pressure and radial pressure
!>                          gradient in the radial direction.
!>  \param[in]  jkxbsubumns Fourier magnetic pressure in the poloidal direction.
!>  \param[in]  jkxbsubvmns Fourier magnetic pressure in the toroidal direction.
!>  \param[out] fsubsmncf   F_s component of the MHD force.
!>  \param[out] fsubumnsf   F_u component of the MHD force.
!>  \param[out] fsubvmnsf   F_v component of the MHD force.
!*******************************************************************************
       SUBROUTINE get_force_harmonics(pmncf, pmncf_ds,                              &
                                      jkxbsubsmnc, jkxbsubumns, jkxbsubvmns,        &
                                      fsubsmncf, fsubumnsf, fsubvmnsf)
       USE quantities, dum1=>fsubsmncf, dum2=>fsubumnsf, dum3=>fsubvmnsf
       USE hessian, ONLY: muPar, l_Compute_Hessian
       USE descriptor_mod, ONLY: iam
       USE evolution, ONLY: l_push_s, l_push_u, l_push_v, l_push_edge
       USE shared_data, ONLY: gc, nprecon, l_PrintOriginForces, col_scale
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(in)  ::     &
           pmncf, pmncf_ds, jkxbsubsmnc, jkxbsubumns, jkxbsubvmns
       REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(out)  ::    &
           fsubsmncf, fsubumnsf, fsubvmnsf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER     :: m0=0, m1=1, m2=2, moff=1
       REAL(rprec), PARAMETER :: p5 = 0.5_dp
       INTEGER                :: m, n, js, noff 
       REAL(rprec)            :: f1, f2, f3, r0

       INTEGER :: i, j, k, nsmin, nsmax, ISTAT 
!-----------------------------------------------
!
!      NOTE: there is an m=1 part of pmncf retained at js=1, which is REALLY at js=hs/2
!      It contributes a spurious piece to fsubvmnsf, which is correctly zeroed
!
       gc = 0                                                            ! zero all forces

       fsubsmncf = -pmncf_ds + jkxbsubsmnc                              ! J(K^u*B^v - K^v*B^u) - dpds
       DO m = 0, mpol
          DO n = -ntor, ntor       
             fsubumnsf(m,n,:) = jkxbsubumns(m,n,:) + m*pmncf(m,n,:)     ! J(K^v*B^s - K^s*B^v) - dpdu  m -> -(d/du)COS(mu+nv)
             fsubvmnsf(m,n,:) = jkxbsubvmns(m,n,:) + n*nfp*pmncf(m,n,:) ! J(K^s*B^u - K^u*B^s) - dpdv  n*nfp -> -(d/dv)COS(mu+n*NFP*v)
           ENDDO
       ENDDO

 !      DO js=1,ns
 !      WRITE (60+iam, *) 'JS: ',js
 !      WRITE (60+iam, '(1p,2e12.4)') fsubumnsf(:,:,js), dum2(:,:,js)
 !      END DO
 !      CALL FLUSH(60+iam)
 !      STOP
!
!      ADD SYMMETRIC DAMPING TERMS -B v dot B TO SUPPRESS NULL SPACE
!      WORKS BEST: comment out IF statement below, allowing -mu v|| term
!      to be added to BOTH linear and nonlinear forces. However, must
!      guarantee that v->0 at end for this to work well

!CALL funct_island WITH vsupXmn...=0 to get un-augmented MHD force
       IF (l_Compute_Hessian .AND. muPar.NE.zero) THEN                  !Always (lin/nonlin) use mu|| = 0 EXCEPT for Hessian calculation
          CALL AddParDamping(fsubsmncf, fsubumnsf, fsubvmnsf, jkxbsubvmns)
       END IF

       fsubsmncf(m0,-ntor:-1,:) = 0                                     !redundant COS modes
       fsubumnsf(m0,-ntor:0,:) = 0                                      !redundant SIN modes
       fsubvmnsf(m0,-ntor:0,:) = 0                                      !"                 "


      IF (iam .eq. 0) THEN
       IF (l_PrintOriginForces) THEN
          fbdy(1) = SQRT(SUM(fsubsmncf(m1,:,1)**2))
          fbdy(2) = SQRT(SUM(fsubsmncf(m1,:,2)**2))
          fbdy(3) = SQRT(SUM(fsubsmncf(m0,:,2)**2) + SUM(fsubsmncf(m2:,:,2)**2))
          WRITE(*,121) fbdy(1),fbdy(2),fbdy(3)
          fbdy(4) = SQRT(SUM(fsubumnsf(m1,:,1)**2))
          fbdy(5) = SQRT(SUM(fsubumnsf(m1,:,2)**2))
          fbdy(6) = SQRT(SUM(fsubumnsf(m0,:,2)**2) + SUM(fsubumnsf(m2:,:,2)**2))
          WRITE(*,122) fbdy(4),fbdy(5),fbdy(6)
          fbdy(7) = SQRT(SUM(fsubvmnsf(m0,:,1)**2))
          fbdy(8) = SQRT(SUM(fsubvmnsf(m0,:,2)**2))
          fbdy(9) = SQRT(SUM(fsubvmnsf(m1:,:,2)**2))
          WRITE(*,123) fbdy(7),fbdy(8),fbdy(9)
          fbdy(10) = SQRT(SUM(fsubumnsf(:,:,ns)**2))
          fbdy(11) = SQRT(SUM(fsubumnsf(:,:,ns-1)**2))
          WRITE(*,124) fbdy(10),fbdy(11)
          fbdy(12) = SQRT(SUM(fsubvmnsf(:,:,ns)**2))
          fbdy(13) = SQRT(SUM(fsubvmnsf(:,:,ns-1)**2))
          WRITE(*,125) fbdy(12),fbdy(13)
         END IF
       END IF
 121   FORMAT(/,' fs(1,m=1): ',1p,e12.3,' fs(2,m=1): ',1pe12.3,' fs(2,m!=1):',1pe12.3)
 122   FORMAT(  ' fu(1,m=1): ',1p,e12.3,' fu(2,m=1): ',1pe12.3,' fu(2,m!=1):',1pe12.3)
 123   FORMAT(  ' fv(1,m=0): ',1p,e12.3,' fv(2,m=0): ',1pe12.3,' fv(2,m>0): ',1pe12.3)
 124   FORMAT(  ' fu(ns)   : ',1p,e12.3,' fu(ns-1):  ',1pe12.3)
 125   FORMAT(  ' fv(ns)   : ',1p,e12.3,' fv(ns-1):  ',1pe12.3,/)

!      ORIGIN BOUNDARY CONDITIONS
       fsubsmncf(m0,:,1) = 0; fsubsmncf(m2:,:,1) = 0
!       r0 = hs_i/2
!       IF (.not.l_natural) fsubumnsf(m1,:,1) = fsubumnsf(m1,:,1)-r0*fsubsmncf(m1,:,1) 
       IF (.not.l_push_s) THEN 
          fsubsmncf(m1,:,1) = 0
       ELSE
          fsubsmncf(m1,:,1) = p5*fsubsmncf(m1,:,1)
       END IF

       fsubumnsf(m0,:,1) = 0; fsubumnsf(m2:,:,1) = 0
       IF (.not.l_push_u) THEN
          fsubumnsf(m1,:,1) = 0
       ELSE
          fsubumnsf(m1,:,1) = p5*fsubumnsf(m1,:,1)
       END IF

       fsubvmnsf(m1:,:,1) = 0

       IF (.not.l_push_v) THEN
          fsubvmnsf(m0,:,1) = 0
       ELSE
          fsubvmnsf(m0,:,1) = p5*fsubvmnsf(m0,:,1)
       END IF
 
!      EDGE BOUNDARY CONDITIONS (factor of p5 due to one-sided average at edge in var principle)
       fsubsmncf(:,:,ns) = 0
       IF (.not.l_push_edge) THEN
          fsubumnsf(:,:,ns) = 0
          fsubvmnsf(:,:,ns) = 0
       ELSE
          fsubumnsf(:,:,ns) = p5*fsubumnsf(:,:,ns)
          fsubvmnsf(:,:,ns) = p5*fsubvmnsf(:,:,ns)
       END IF

       fsubsmncf = signjac*fsubsmncf*col_scale(:,:,1,:)
       fsubumnsf = signjac*fsubumnsf*col_scale(:,:,2,:)
       fsubvmnsf = signjac*fsubvmnsf*col_scale(:,:,3,:)
       
       END SUBROUTINE get_force_harmonics

#if defined(MRC)
!*******************************************************************************
!>  \brief Computed harmonics of the MHD force serially for the other parity.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = K X B - Grad(p)                                                   (1)
!>
!>  where K is the plasma current density. In curalinear coodinates, the forces
!>  components becomes
!>
!>     F_s = J(K^u*B^v - K^v*B^u) - dpds                                    (2a)
!>
!>     F_u = J(K^v*B^s - K^s*B^v) - dpdu                                    (2b)
!>
!>     F_v = J(K^s*B^u - K^u*B^s) - dpdv                                    (2c)
!>
!>  where J is the coordinate system Jacobian. Note this routine only works on a
!>  single parity.
!>
!>  \param[in]  pmnsf       Fourier pressure.
!>  \param[in]  pmns_ds     Radial derivative of fourier pressure.
!>  \param[in]  jkxbsubsmns Fourier magnetic pressure and radial pressure
!>                          gradient in the radial direction.
!>  \param[in]  jkxbsubumnc Fourier magnetic pressure in the poloidal direction.
!>  \param[in]  jkxbsubvmnc Fourier magnetic pressure in the toroidal direction.
!>  \param[out] fsubsmnsf   F_s component of the MHD force.
!>  \param[out] fsubumncf   F_u component of the MHD force.
!>  \param[out] fsubvmncf   F_v component of the MHD force.
!*******************************************************************************
       SUBROUTINE get_force_harmonics_other(pmnsf, pmnsf_ds,                        &
                                      jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc,        &
                                      fsubsmnsf, fsubumncf, fsubvmncf)
       USE quantities, dum1=>fsubsmnsf, dum2=>fsubumncf, dum3=>fsubvmncf
       USE hessian, ONLY: muPar, l_Compute_Hessian
       USE descriptor_mod, ONLY: iam
       USE evolution, ONLY: l_push_s, l_push_u, l_push_v, l_push_edge
       USE shared_data, ONLY: gc, nprecon, l_PrintOriginForces, col_scale
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(in)  ::     &
           pmnsf, pmnsf_ds, jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc
       REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(out)  ::    &
           fsubsmnsf, fsubumncf, fsubvmncf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER     :: m0=0, m1=1, m2=2, moff=1
       REAL(rprec), PARAMETER :: p5 = 0.5_dp
       INTEGER                :: m, n, js, noff
       REAL(rprec)            :: f1, f2, f3, r0

       INTEGER :: i, j, k, nsmin, nsmax, ISTAT
!-----------------------------------------------
!
!      NOTE: there is an m=1 part of pmncf retained at js=1, which is REALLY at js=hs/2
!      It contributes a spurious piece to fsubvmnsf, which is correctly zeroed
!
       gc = 0                                                            ! zero all forces

       fsubsmnsf = -pmnsf_ds + jkxbsubsmns                              ! J(K^u*B^v - K^v*B^u) - dpds
       DO m = 0, mpol
          DO n = -ntor, ntor
             fsubumncf(m,n,:) = jkxbsubumnc(m,n,:) - m*pmnsf(m,n,:)     ! J(K^v*B^s - K^s*B^v) - dpdu  m -> (d/du)SIN(mu+nv)
             fsubvmncf(m,n,:) = jkxbsubvmnc(m,n,:) - n*nfp*pmnsf(m,n,:) ! J(K^s*B^u - K^u*B^s) - dpdv  n*nfp -> (d/dv)SIN(mu+n*NFP*v)
          ENDDO
       ENDDO

!      DO js=1,ns
!      WRITE (60+iam, *) 'JS: ',js
!      WRITE (60+iam, '(1p,2e12.4)') fsubumnsf(:,:,js), dum2(:,:,js)
!      END DO
!      CALL FLUSH(60+iam)
!      STOP
!
!      ADD SYMMETRIC DAMPING TERMS -B v dot B TO SUPPRESS NULL SPACE
!      WORKS BEST: comment out IF statement below, allowing -mu v|| term
!      to be added to BOTH linear and nonlinear forces. However, must
!      guarantee that v->0 at end for this to work well

!CALL funct_island WITH vsupXmn...=0 to get un-augmented MHD force
       IF (l_Compute_Hessian .AND. muPar.NE.zero) THEN                  !Always (lin/nonlin) use mu|| = 0 EXCEPT for Hessian calculation
          CALL AddParDamping_other(fsubsmnsf, fsubumncf, fsubvmncf, jkxbsubvmnc)
       END IF

       fsubsmnsf(m0,-ntor:0,:) = 0                                     !redundant SIN modes
       fsubumncf(m0,-ntor:-1,:) = 0                                    !redundant COS modes
       fsubvmncf(m0,-ntor:-1,:) = 0                                    !"                 "


       IF (iam .eq. 0) THEN
          IF (l_PrintOriginForces) THEN
             fbdy(1) = SQRT(SUM(fsubsmnsf(m1,:,1)**2))
             fbdy(2) = SQRT(SUM(fsubsmnsf(m1,:,2)**2))
             fbdy(3) = SQRT(SUM(fsubsmnsf(m0,:,2)**2) + SUM(fsubsmnsf(m2:,:,2)**2))
             WRITE(*,121) fbdy(1),fbdy(2),fbdy(3)
             fbdy(4) = SQRT(SUM(fsubumncf(m1,:,1)**2))
             fbdy(5) = SQRT(SUM(fsubumncf(m1,:,2)**2))
             fbdy(6) = SQRT(SUM(fsubumncf(m0,:,2)**2) + SUM(fsubumncf(m2:,:,2)**2))
             WRITE(*,122) fbdy(4),fbdy(5),fbdy(6)
             fbdy(7) = SQRT(SUM(fsubvmncf(m0,:,1)**2))
             fbdy(8) = SQRT(SUM(fsubvmncf(m0,:,2)**2))
             fbdy(9) = SQRT(SUM(fsubvmncf(m1:,:,2)**2))
             WRITE(*,123) fbdy(7),fbdy(8),fbdy(9)
             fbdy(10) = SQRT(SUM(fsubumncf(:,:,ns)**2))
             fbdy(11) = SQRT(SUM(fsubumncf(:,:,ns-1)**2))
             WRITE(*,124) fbdy(10),fbdy(11)
             fbdy(12) = SQRT(SUM(fsubvmncf(:,:,ns)**2))
             fbdy(13) = SQRT(SUM(fsubvmncf(:,:,ns-1)**2))
             WRITE(*,125) fbdy(12),fbdy(13)
          END IF
       END IF
 121   FORMAT(/,' fs(1,m=1): ',1p,e12.3,' fs(2,m=1): ',1pe12.3,' fs(2,m!=1):',1pe12.3)
 122   FORMAT(  ' fu(1,m=1): ',1p,e12.3,' fu(2,m=1): ',1pe12.3,' fu(2,m!=1):',1pe12.3)
 123   FORMAT(  ' fv(1,m=0): ',1p,e12.3,' fv(2,m=0): ',1pe12.3,' fv(2,m>0): ',1pe12.3)
 124   FORMAT(  ' fu(ns)   : ',1p,e12.3,' fu(ns-1):  ',1pe12.3)
 125   FORMAT(  ' fv(ns)   : ',1p,e12.3,' fv(ns-1):  ',1pe12.3,/)

!      ORIGIN BOUNDARY CONDITIONS
       fsubsmnsf(m0,:,1) = 0; fsubsmnsf(m2:,:,1) = 0
!       r0 = hs_i/2
!       IF (.not.l_natural) fsubumnsf(m1,:,1) = fsubumnsf(m1,:,1)-r0*fsubsmncf(m1,:,1)
       IF (.not.l_push_s) THEN
          fsubsmnsf(m1,:,1) = 0
       ELSE
          fsubsmnsf(m1,:,1) = p5*fsubsmnsf(m1,:,1)
       END IF

       fsubumncf(m0,:,1) = 0; fsubumncf(m2:,:,1) = 0
       IF (.not.l_push_u) THEN
          fsubumncf(m1,:,1) = 0
       ELSE
          fsubumncf(m1,:,1) = p5*fsubumncf(m1,:,1)
       END IF

       fsubvmncf(m1:,:,1) = 0

       IF (.not.l_push_v) THEN
          fsubvmncf(m0,:,1) = 0
       ELSE
          fsubvmncf(m0,:,1) = p5*fsubvmncf(m0,:,1)
       END IF

!      EDGE BOUNDARY CONDITIONS (factor of p5 due to one-sided average at edge in var principle)
       fsubsmnsf(:,:,ns) = 0
       IF (.not.l_push_edge) THEN
          fsubumncf(:,:,ns) = 0
          fsubvmncf(:,:,ns) = 0
       ELSE
          fsubumncf(:,:,ns) = p5*fsubumncf(:,:,ns)
          fsubvmncf(:,:,ns) = p5*fsubvmncf(:,:,ns)
       END IF

       fsubsmnsf = signjac*fsubsmnsf*col_scale(:,:,1,:)
       fsubumncf = signjac*fsubumncf*col_scale(:,:,2,:)
       fsubvmncf = signjac*fsubvmncf*col_scale(:,:,3,:)

       END SUBROUTINE get_force_harmonics_other
#endif

!*******************************************************************************
!>  \brief Add parallel damping force.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = -B*v.B*mu_||*mu_norm                                              (1)
!>
!>  where mu is the parallel
!>
!>     F_s = B_s*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>     F_u = B_u*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>     F_v = B_v*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>  where J is the coordinate system Jacobian. Note this routine only works on a
!>  single parity.
!>
!>  \param[inout] fsubsmncf F_s component of the MHD force.
!>  \param[inout] fsubumnsf F_u component of the MHD force.
!>  \param[inout] fsubvmnsf F_v component of the MHD force.
!>  \param        work      Temporary work buffer.
!*******************************************************************************
       SUBROUTINE addpardamping (fsubsmncf, fsubumnsf, fsubvmnsf, work)
       USE quantities, dum1=>fsubsmncf, dum2=>fsubumnsf, dum3=>fsubvmnsf
       USE hessian, ONLY: muPar, mupar_norm
       USE fourier, ONLY: tomnsp
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec) :: work(0:mpol,-ntor:ntor,ns)
       REAL(rprec), INTENT(INOUT), DIMENSION(0:mpol,-ntor:ntor,ns) ::   &
           fsubsmncf, fsubumnsf, fsubvmnsf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: js, istat
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: vdotB, tmp
       INTEGER :: i, j, k
!-----------------------------------------------
!
!      Compute parallel (to B) damping force and add it to (linearized) forces
!
       ALLOCATE (vdotB(ntheta,nzeta,ns), tmp(ntheta,nzeta,ns), stat=istat)

       IF (istat .ne. 0) STOP 'ALLOCATION ERROR IN ADDPARDAMPING'

!
!      compute (sin-parity) v * B: note the -signJac factor is in mupar_norm
!
       bsubuijf(:,:,1) = 0
       vdotB = jvsupsijf*bsubsijf + jvsupuijf*bsubuijf + jvsupvijf*bsubvijf

       DO js = 1,ns
          vdotB(:,:,js) = muPar*mupar_norm(js)*vdotB(:,:,js)
       END DO

       tmp = bsubsijf*vdotB                                             !COS SYM
       CALL tomnsp(tmp, work, 0, 0)
       fsubsmncf = fsubsmncf + work

       tmp = bsubuijf*vdotB                                             !SIN SYM
       CALL tomnsp(tmp, work, 1, 0)
       fsubumnsf = fsubumnsf + work

       tmp = bsubvijf*vdotB                                             !SIN SYM
       CALL tomnsp(tmp, work, 1, 0)
       fsubvmnsf = fsubvmnsf + work

       DEALLOCATE (vdotB, tmp)

       END SUBROUTINE addpardamping

#if defined(MRC)
!*******************************************************************************
!>  \brief Add parallel damping force for the other parity.
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = -B*v.B*mu_||*mu_norm                                              (1)
!>
!>  where mu is the parallel
!>
!>     F_s = B_s*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>     F_u = B_u*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>     F_v = B_v*J*(v^s*B_s + v^u*B_u + v^v*B_v)                            (2a)
!>
!>  where J is the coordinate system Jacobian. Note this routine only works on a
!>  single parity.
!>
!>  \param[inout] fsubsmnsf F_s component of the MHD force.
!>  \param[inout] fsubumncf F_u component of the MHD force.
!>  \param[inout] fsubvmncf F_v component of the MHD force.
!>  \param        work      Temporary work buffer.
!*******************************************************************************
       SUBROUTINE addpardamping_other(fsubsmnsf, fsubumncf, fsubvmncf, work)
       USE quantities, dum1=>fsubsmnsf, dum2=>fsubumncf, dum3=>fsubvmncf
       USE hessian, ONLY: muPar, mupar_norm
       USE fourier, ONLY: tomnsp, fourier_parity_odd, fourier_parity_even, &
                          fourier_full_grid
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec) :: work(0:mpol,-ntor:ntor,ns)
       REAL(rprec), INTENT(INOUT), DIMENSION(0:mpol,-ntor:ntor,ns) ::   &
           fsubsmnsf, fsubumncf, fsubvmncf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: js, istat
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: vdotB, tmp
       INTEGER :: i, j, k
!-----------------------------------------------
!
!      Compute parallel (to B) damping force and add it to (linearized) forces
!
       ALLOCATE (vdotB(ntheta,nzeta,ns), tmp(ntheta,nzeta,ns), stat=istat)

       IF (istat .ne. 0) STOP 'ALLOCATION ERROR IN ADDPARDAMPING'

!
!      compute (sin-parity) v * B: note the -signJac factor is in mupar_norm
!
       bsubuijf(:,:,1) = 0
       vdotB = jvsupsijf*bsubsijf + jvsupuijf*bsubuijf + jvsupvijf*bsubvijf

       DO js = 1,ns
          vdotB(:,:,js) = muPar*mupar_norm(js)*vdotB(:,:,js)
       END DO

       tmp = bsubsijf*vdotB                                             !SIN SYM
          CALL tomnsp(tmp, work, fourier_parity_odd, fourier_full_grid)
       fsubsmnsf = fsubsmnsf + work

       tmp = bsubuijf*vdotB                                             !COS SYM
          CALL tomnsp(tmp, work, fourier_parity_even, fourier_full_grid)
       fsubumncf = fsubumncf + work

       tmp = bsubvijf*vdotB                                             !COS SYM
          CALL tomnsp(tmp, work, fourier_parity_even, fourier_full_grid)
       fsubvmncf = fsubvmncf + work

       DEALLOCATE (vdotB, tmp)

       END SUBROUTINE addpardamping_other
#endif
