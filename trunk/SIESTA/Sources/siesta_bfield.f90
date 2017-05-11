!>  \brief Module contained subroutines for updating half-grid magnetic fields as part of the SIESTA project
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 21, 2006
       MODULE siesta_bfield

       USE stel_kinds
       USE stel_constants
       USE fourier 
       USE quantities
       USE island_params, ONLY: nu_i
       USE perturbation, ONLY: lresistive, eta_factor, buv_res
       USE timer_mod
       USE shared_data, ONLY: delta_t, fsq_total, l_natural, niter, fsq_res

       CONTAINS
!#if defined(SKS)      
!*******************************************************************************
!>  \brief Parallel routine for updating B^s, B^u, B^v for ideal (and resistive)
!>  perturbations.
!>
!>  Advances half-mesh magnetic field components in nm space from time t to
!>  t + dt. Used to update values of BSUPMNF and BSUPIJF in \ref quantities
!>  module.
!*******************************************************************************
       SUBROUTINE update_bfield_par (l_add_res)
!     
!     PURPOSE: ADVANCES HALF-MESH MAGNETIC FIELD COMPONENTS IN MN SPACE FROM TIME T TO T+DELTA_T
!              (Used to compute dBSUPXMNH for updating BSUPXMNH) 
!              DISCRETIZED VERSION OF FARADAY'S LAW dB/dt = -curl(E)
!
       USE nscalingtools, ONLY: startglobrow, endglobrow 
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)     :: l_add_res
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m0=0, m1=1, m2=2
       INTEGER :: m, n, l, istat, js, nblock
       LOGICAL :: ladd_pert 
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                        &
         esubsijf, esubuijf, esubvijf,                                   &
         esubsmnsf, esubumncf, esubvmncf,                                &
         esubsmnsh, esubumnch, esubvmnch,                                &
         esubsmncf, esubumnsf, esubvmnsf,                                &
         esubsmnch, esubumnsh, esubvmnsh,                                &
         de_uds, de_vds, resistivity
       REAL(dp) :: ton, toff, rho
       REAL(dp) :: delt_cfl, eta_prof, r0
       INTEGER :: nsmin, nsmax, ns_span
 !-----------------------------------------------
       CALL second0(ton)

       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+1)
       ns_span = nsmax-nsmin+1
       ALLOCATE(esubsijf(ntheta,nzeta,nsmin:nsmax),                      &
                esubuijf(ntheta,nzeta,nsmin:nsmax),                      &
                esubvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)   
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_BFIELD' 
!
!      Calculate ideal covariant components of electric field E-ideal = -v X B on full-space mesh
!  
       esubsijf =-(jvsupuijf(:,:,nsmin:nsmax)*bsupvijf0(:,:,nsmin:nsmax) &
                 - jvsupvijf(:,:,nsmin:nsmax)*bsupuijf0(:,:,nsmin:nsmax))  ! -E_s = V^uB^v - V^vB^u  (SINE): jv -> jac*v
       esubuijf =-(jvsupvijf(:,:,nsmin:nsmax)*bsupsijf0(:,:,nsmin:nsmax) &
                 - jvsupsijf(:,:,nsmin:nsmax)*bsupvijf0(:,:,nsmin:nsmax))  ! -E_u = V^vB^s - V^sB^v  (COSINE)
       esubvijf =-(jvsupsijf(:,:,nsmin:nsmax)*bsupuijf0(:,:,nsmin:nsmax) &
                 - jvsupuijf(:,:,nsmin:nsmax)*bsupsijf0(:,:,nsmin:nsmax))  ! -E_v = V^sB^u - V^uB^s  (COSINE)

!      confirm boundary condition: esubu,v(s=1) = 0 (tangential E vanishes at bdy => dB^s = 0)
       IF (nsmax .EQ. ns) THEN
         IF (ANY(esubuijf(:,:,ns) .NE. zero)) THEN
           STOP 'esubuijf(ns) != 0 in UPDATE_BFIELD_par'
         ELSE IF (ANY(esubvijf(:,:,ns) .NE. zero)) THEN
           STOP 'esubvijf(ns) != 0 in UPDATE_BFIELD_par'
         END IF
       END IF

!
!      Note this will LOWER the energy due to eta*|J|||**2 heating
!      esubX(resistive) = eta(JdotB/B^2)*BsubX
!      so the magnetic energy DECREASES due to this term. Note ksubX=jac*JsubX are the
!      covariant components of the current TIMES the jac factor
!
       ladd_pert = ALLOCATED(buv_res)
       IF (lresistive .AND. (l_add_res .or. ladd_pert)) THEN
          delt_cfl = hs_i**2*ABS(eta_factor)
          IF (fsq_total < fsq_res) delt_cfl = delt_cfl*SQRT(fsq_total/fsq_res)

          ALLOCATE(resistivity(ntheta,nzeta,nsmin:nsmax), stat=istat)
          IF (istat .NE. 0) STOP 'ALLOCATION ERROR IN update_bfield'

          DO js = nsmin, nsmax
             rho = hs_i*(js-1)
             eta_prof = rho*rho*(1-rho)
             resistivity(:,:,js) = delt_cfl*eta_prof
          END DO

!Note: need to divide current by 1/jac
          resistivity = resistivity/jacobf(:,:,nsmin:nsmax)
          IF (ladd_pert) THEN
             DO js = nsmin, nsmax
                resistivity(:,:,js) = resistivity(:,:,js)*buv_res(:,:,js)
             END DO
          END IF

!Isotropic resistivity, E ~ eta*J; for ladd_pert, K = B in init_state so this is A|| B perturbation
          DO js = nsmin, nsmax
             esubsijf(:,:,js) = esubsijf(:,:,js) + resistivity(:,:,js)*ksubsijf(:,:,js)
             esubuijf(:,:,js) = esubuijf(:,:,js) + resistivity(:,:,js)*ksubuijf(:,:,js)
             esubvijf(:,:,js) = esubvijf(:,:,js) + resistivity(:,:,js)*ksubvijf(:,:,js)
          END DO

          DEALLOCATE (resistivity)
       ENDIF 

!
!Calculate harmonics of electric field
! 
       ALLOCATE(esubsmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                esubumncf(0:mpol,-ntor:ntor,nsmin:nsmax),                & 
                esubvmncf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'

       CALL tomnsp_par(esubsijf, esubsmnsf, f_sin, f_full, 1, ns_span)
       CALL tomnsp_par(esubuijf, esubumncf, f_cos, f_full, 1, ns_span)
       CALL tomnsp_par(esubvijf, esubvmncf, f_cos, f_full, 1, ns_span)

       
       ALLOCATE(de_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                de_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'

       ALLOCATE(esubsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                esubumnch(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                esubvmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'  
       
       !!On exit, esubXh valid on nsmin+1:nsmax grid
       nblock = SIZE(esubsmnsh,1)*SIZE(esubsmnsh,2)
       CALL to_half_mesh_par(esubsmnsf, esubsmnsh, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(esubumncf, esubumnch, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(esubvmncf, esubvmnch, nblock, nsmin, nsmax)

!      Gradients of esubu,v on half grid j-half=2,ns: note esubu,v(1) used here!
       !!On exit, work1,2 valid on nsmin+1:nsmax grid
       CALL GradientHalf_par(de_uds, esubumncf, nblock, nsmin, nsmax)
       CALL GradientHalf_par(de_vds, esubvmncf, nblock, nsmin, nsmax)

!      Calculate r.h.s. of the evolution equations of all components of jacobh*BSUP?
!      Compute contravariant components of -curl(E) for right side of Faraday's Equation
!
       r0 = hs_i/2
       nsmin=MAX(1,startglobrow)

!IMPOSE NATURAL BOUNDARY CONDITION AT FIRST HALF-GRID POINT (SPH 04-28-16)
       IF (nsmin.EQ.1 .AND. .NOT.l_natural) THEN
!       IF (.FALSE.) THEN
          djbsupsmnsh(m1,:,1) = (esubsmnsh(m1,:,2)*r0-esubumnch(m1,:,2))/2
          esubsmnsh(m1,:,2) = esubsmnsh(m1,:,2) - djbsupsmnsh(m1,:,1)/r0
          esubumnch(m1,:,2) = esubumnch(m1,:,2) + djbsupsmnsh(m1,:,1)
          de_uds(m1,:,2) = de_uds(m1,:,2) - djbsupsmnsh(m1,:,1)/r0         !Preserves div B = 0
       END IF

       DO m = 0, mpol
         DO n = -ntor, ntor
           djbsupsmnsh(m,n,nsmin:nsmax) =      m*esubvmnch(m,n,nsmin:nsmax) & ! -m and -n*NFP from deriv. COS(m*theta+n*NFP*phi)
                                          -n*nfp*esubumnch(m,n,nsmin:nsmax)
           djbsupumnch(m,n,nsmin:nsmax) = -n*nfp*esubsmnsh(m,n,nsmin:nsmax) + de_vds(m,n,nsmin:nsmax) ! n*NFP from deriv. SIN(m*theta+n*NFP*phi)             
           djbsupvmnch(m,n,nsmin:nsmax) =      m*esubsmnsh(m,n,nsmin:nsmax) - de_uds(m,n,nsmin:nsmax) ! m from deriv. SIN(m*theta+n*NFP*phi)
          ENDDO
       ENDDO

!CONSTRAINT AT ORIGIN [jbsups(m=1) = r0*jbsupu(m=1) <=> esubumnsh = r0*esubsmnch]
!Origin bc for m=1: e_u = r0*e_s (Follows form v X B at origin)
!This suffices to make B^s = r0*B^u
       IF (nsmin .EQ. 1) THEN

          djbsupsmnsh(:,:,1) = 0
          djbsupumnch(:,:,1) = 0
          djbsupvmnch(:,:,1) = 0
      
!      SYMMETRIC FORM (THESE AREN'T USED ANYWHERE)
          djbsupumnch(m1,:,1) = djbsupumnch(m1,:,2)
          IF (l_natural) THEN
             djbsupsmnsh(m1,:,1) = r0*djbsupumnch(m1,:,1)
          ELSE
             djbsupsmnsh(m1,:,1) = djbsupsmnsh(m1,:,2)
          END IF

          djbsupvmnch(m0,:,1) = djbsupvmnch(m0,:,2)
 
       END IF
 
       DEALLOCATE (esubsmnsf, esubumncf, esubvmncf)
       DEALLOCATE (esubsmnsh, esubumnch, esubvmnch)    

!
!      Calculate INCREMENT of magnetic field harmonics: use Euler scheme with DT given by VSUBMN advance equations.
!      Results computed here for HALF-GRID perturbations in MN space are used in the calculation of the force and
!      to advance the B's in UPDATE_STATE (in siesta_state module)
!
!      B^S (SIN) 
       djbsupsmnsh(:,:,nsmin:nsmax) = delta_t*djbsupsmnsh(:,:,nsmin:nsmax)
!
!      B^U (COS) 
       djbsupumnch(:,:,nsmin:nsmax) = delta_t*djbsupumnch(:,:,nsmin:nsmax)
!
!      B^V (COS) 
       djbsupvmnch(:,:,nsmin:nsmax) = delta_t*djbsupvmnch(:,:,nsmin:nsmax)

       IF (.NOT.l_asym) GOTO 1000

       nsmin=MAX(1,startglobrow-1)

       ALLOCATE(esubsmncf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                esubumnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                esubvmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'

       CALL tomnsp_par(esubsijf, esubsmncf, f_cos, f_full, 1, ns_span)
       CALL tomnsp_par(esubuijf, esubumnsf, f_sin, f_full, 1, ns_span)
       CALL tomnsp_par(esubvijf, esubvmnsf, f_sin, f_full, 1, ns_span)

       ALLOCATE(esubsmnch(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                esubumnsh(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                esubvmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'  
       
       !!On exit, esubXh valid on nsmin+1:nsmax grid
       nblock = SIZE(esubsmnch,1)*SIZE(esubsmnch,2)
       CALL to_half_mesh_par(esubsmncf, esubsmnch, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(esubumnsf, esubumnsh, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(esubvmnsf, esubvmnsh, nblock, nsmin, nsmax)

!      Gradients of esubu,v on half grid j-half=2,ns: note esubu,v(1) used here!
       !!On exit, work1,2 valid on nsmin+1:nsmax grid
       CALL GradientHalf_par(de_uds, esubumnsf, nblock, nsmin, nsmax)
       CALL GradientHalf_par(de_vds, esubvmnsf, nblock, nsmin, nsmax)

!      Calculate r.h.s. of the evolution equations of all components of jacobh*BSUP?
!      Compute contravariant components of -curl(E) for right side of Faraday's Equation
!
       nsmin=MAX(1,startglobrow)
       r0 = -r0

!IMPOSE NATURAL BOUNDARY CONDITION AT FIRST HALF-GRID POINT (SPH 04-28-16)
       IF (nsmin.EQ.1 .AND. .NOT.l_natural) THEN
!       IF (.FALSE.) THEN
          djbsupsmnch(m1,:,1) = (esubsmnch(m1,:,2)*r0-esubumnsh(m1,:,2))/2
          esubsmnch(m1,:,2) = esubsmnch(m1,:,2) - djbsupsmnch(m1,:,1)/r0
          esubumnsh(m1,:,2) = esubumnsh(m1,:,2) + djbsupsmnch(m1,:,1)
          de_uds(m1,:,2) = de_uds(m1,:,2) - djbsupsmnch(m1,:,1)/r0         !Preserves div B = 0
       END IF

       DO m = 0, mpol
         DO n = -ntor, ntor
           djbsupsmnch(m,n,nsmin:nsmax) =     -m*esubvmnsh(m,n,nsmin:nsmax) & ! -m and -n*NFP from deriv. SIN(m*theta+n*NFP*phi)
                                          +n*nfp*esubumnsh(m,n,nsmin:nsmax)
           djbsupumnsh(m,n,nsmin:nsmax) =  n*nfp*esubsmnch(m,n,nsmin:nsmax) + de_vds(m,n,nsmin:nsmax) ! n*NFP from deriv. COS(m*theta+n*NFP*phi)             
           djbsupvmnsh(m,n,nsmin:nsmax) =     -m*esubsmnch(m,n,nsmin:nsmax) - de_uds(m,n,nsmin:nsmax) ! m from deriv. COS(m*theta+n*NFP*phi) 
          ENDDO
       ENDDO

!CONSTRAINT AT ORIGIN [jbsups(m=1) = r0*jbsupu(m=1) <=> esubumnch = r0*esubsmnsh]
!Origin bc for m=1: e_u = r0*e_s (Follows form v X B at origin)
!This suffices to make B^s = r0*B^u
       IF (nsmin .EQ. 1) THEN

          djbsupsmnch(:,:,1) = 0
          djbsupumnsh(:,:,1) = 0
          djbsupvmnsh(:,:,1) = 0
      
!      SYMMETRIC FORM (THESE AREN'T USED ANYWHERE)
          djbsupumnsh(m1,:,1) = djbsupumnsh(m1,:,2)
          IF (l_natural) THEN
             djbsupsmnch(m1,:,1) = r0*djbsupumnsh(m1,:,1)
          ELSE
             djbsupsmnch(m1,:,1) = djbsupsmnch(m1,:,2)
          END IF

          djbsupvmnsh(m0,:,1) = djbsupvmnsh(m0,:,2)
 
       END IF
 
       DEALLOCATE (esubsmncf, esubumnsf, esubvmnsf)
       DEALLOCATE (esubsmnch, esubumnsh, esubvmnsh)    
!
!      Calculate INCREMENT of magnetic field harmonics: use Euler scheme with DT given by VSUBMN advance equations.
!      Results computed here for HALF-GRID perturbations in MN space are used in the calculation of the force and
!      to advance the B's in UPDATE_STATE (in siesta_state module)
!
!      B^S (COS)
       djbsupsmnch(:,:,nsmin:nsmax) = delta_t*djbsupsmnch(:,:,nsmin:nsmax)
!
!      B^U (SIN)
       djbsupumnsh(:,:,nsmin:nsmax) = delta_t*djbsupumnsh(:,:,nsmin:nsmax)
!
!      B^V (SIN)
       djbsupvmnsh(:,:,nsmin:nsmax) = delta_t*djbsupvmnsh(:,:,nsmin:nsmax)

 1000  CONTINUE

       DEALLOCATE(esubsijf, esubuijf, esubvijf, stat=istat)
       DEALLOCATE(de_uds, de_vds)    

       CALL second0(toff)
       time_update_bfield = time_update_bfield + (toff-ton)

       END SUBROUTINE update_bfield_par
!#endif

!>  \brief Serial routine for updating B^s, B^u, B^v for ideal (and resistive) perturbations
       SUBROUTINE update_bfield (l_add_res)
!     
!     PURPOSE: ADVANCES HALF-MESH MAGNETIC FIELD COMPONENTS IN MN SPACE FROM TIME T TO T+DELTA_T
!              (Used to updates values of BSUPMNF and BSUPIJF in QUANTITIES MODULE)
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)     :: l_add_res

       CALL update_bfield_par (l_add_res)

       END SUBROUTINE update_bfield

       END MODULE siesta_bfield
