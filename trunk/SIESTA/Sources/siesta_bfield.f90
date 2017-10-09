!>  \brief Module contained subroutines for updating half-grid magnetic fields as part of the SIESTA project
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 21, 2006
       MODULE siesta_bfield

       USE stel_kinds
       USE stel_constants
       USE fourier 
       USE nscalingtools, ONLY: startglobrow, endglobrow
       USE timer_mod
       USE utilities, ONLY: GradientHalf, to_half_mesh
       IMPLICIT NONE
       
       PRIVATE
       INTEGER, PARAMETER :: m0=0, m1=1, m2=2, isym=0, iasym=1
       INTEGER            :: nsmin, nsmax

       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                        &
         esubsijf, esubuijf, esubvijf,                                   &
         esubsmnf, esubumnf, esubvmnf,                                   &
         esubsmnh, esubumnh, esubvmnh, de_uds, de_vds

       PUBLIC :: update_bfield

       CONTAINS
!*******************************************************************************
!>  \brief Routine (both parallel, serial) for updating B^s, B^u, B^v 
!>  for both ideal and resistive perturbations.
!>
!>  Advances half-mesh magnetic field components in nm space from time t to
!>  t + dt. Used to update values of BSUPMNF and BSUPIJF in \ref quantities
!>  module.
!*******************************************************************************
       SUBROUTINE update_bfield (l_add_res)
!     
!     PURPOSE: ADVANCES HALF-MESH MAGNETIC FIELD COMPONENTS IN MN SPACE FROM TIME T TO T+DELTA_T
!              (Used to compute dBSUPXMNH for updating BSUPXMNH) 
!              DISCRETIZED VERSION OF FARADAY'S LAW dB/dt = -curl(E)
!
       USE quantities
       USE descriptor_mod, ONLY: INHESSIAN
       USE perturbation, ONLY: lresistive, eta_factor, buv_res
       USE shared_data, ONLY: lasym, fsq_res, fsq_total
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)     :: l_add_res
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: istat, js
       LOGICAL  :: ladd_pert 
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: resistivity
       REAL(dp) :: ton, toff, rho
       REAL(dp) :: delt_cfl, eta_prof, r0
 !-----------------------------------------------
       CALL second0(ton)

       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+1)

!
!Calculate covariant components of (non-resistive) electric field E = -v X B on full-space mesh
!  
       ALLOCATE(esubsijf(ntheta,nzeta,nsmin:nsmax),                      &
                esubuijf(ntheta,nzeta,nsmin:nsmax),                      &
                esubvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)   
       IF (istat .NE. 0) CALL ErrorAbort('Allocation1 failed in UPDATE_BFIELD')

       ALLOCATE(de_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                de_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) CALL ErrorAbort('Allocation2 failed in UPDATE_BFIELD')

       esubsijf =-(jvsupuijf(:,:,nsmin:nsmax)*bsupvijf0(:,:,nsmin:nsmax) &
                 - jvsupvijf(:,:,nsmin:nsmax)*bsupuijf0(:,:,nsmin:nsmax))  ! -E_s = V^uB^v - V^vB^u  (SINE): jv -> jac*v
       esubuijf =-(jvsupvijf(:,:,nsmin:nsmax)*bsupsijf0(:,:,nsmin:nsmax) &
                 - jvsupsijf(:,:,nsmin:nsmax)*bsupvijf0(:,:,nsmin:nsmax))  ! -E_u = V^vB^s - V^sB^v  (COSINE)
       esubvijf =-(jvsupsijf(:,:,nsmin:nsmax)*bsupuijf0(:,:,nsmin:nsmax) &
                 - jvsupuijf(:,:,nsmin:nsmax)*bsupsijf0(:,:,nsmin:nsmax))  ! -E_v = V^sB^u - V^uB^s  (COSINE)

!confirm boundary condition: esubu,v(s=1) = 0 (tangential E vanishes at bdy => dB^s = 0)
       IF (nsmax .EQ. ns) THEN
         IF (ANY(esubuijf(:,:,ns) .NE. zero)) THEN
           CALL ErrorAbort('esubuijf(ns) != 0 in UPDATE_BFIELD_par')
         ELSE IF (ANY(esubvijf(:,:,ns) .NE. zero)) THEN
           CALL ErrorAbort('esubvijf(ns) != 0 in UPDATE_BFIELD_par')
         END IF
       END IF

!
!Note this will LOWER the energy due to eta*|J|||**2 heating
!      esubX(resistive) = eta(JdotB/B^2)*BsubX
!so the magnetic energy DECREASES due to this term. Note ksubX=JsubX are the
!covariant components of the current
!
       ladd_pert = ALLOCATED(buv_res)
       IF (lresistive .AND. (l_add_res .or. ladd_pert)) THEN
          delt_cfl = hs_i**2*ABS(eta_factor)
          IF (fsq_total < fsq_res) delt_cfl = delt_cfl*SQRT(fsq_total/fsq_res)

          ALLOCATE(resistivity(ntheta,nzeta,nsmin:nsmax), stat=istat)
          IF (istat .NE. 0) CALL ErrorAbort('Allocation5 failed in UPDATE_BFIELD')

          DO js = nsmin, nsmax
             rho = hs_i*(js-1)
             eta_prof = rho*rho*(1-rho)
             resistivity(:,:,js) = delt_cfl*eta_prof
          END DO

          IF (ladd_pert) THEN
             resistivity = resistivity*buv_res(:,:,nsmin:nsmax)
          END IF

!Isotropic resistivity, E ~ eta*J; for ladd_pert=TRUE, K ~ B in init_state
          esubsijf = esubsijf + resistivity*ksubsijf(:,:,nsmin:nsmax)
          esubuijf = esubuijf + resistivity*ksubuijf(:,:,nsmin:nsmax)
          esubvijf = esubvijf + resistivity*ksubvijf(:,:,nsmin:nsmax)

          DEALLOCATE (resistivity)
       ENDIF 

!Allocate working arrays used in Faraday
       ALLOCATE(esubsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
                esubumnf(0:mpol,-ntor:ntor,nsmin:nsmax),                 & 
                esubvmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
                stat=istat)  
       IF (istat .NE. 0) CALL ErrorAbort('Allocation3 failed in UPDATE_BFIELD')

       ALLOCATE(esubsmnh(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
                esubumnh(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
                esubvmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) CALL ErrorAbort('Allocation4 failed in UPDATE_BFIELD')

!Update Bfield using Faraday's Law
       CALL Faraday(djbsupsmnsh, djbsupumnch, djbsupvmnch, isym)
       IF (lasym)                                                        &
       CALL Faraday(djbsupsmnch, djbsupumnsh, djbsupvmnsh, iasym)

       DEALLOCATE (esubsijf, esubuijf, esubvijf, de_uds, de_vds,         &    
                   esubsmnf, esubumnf, esubvmnf,                         & 
                   esubsmnh, esubumnh, esubvmnh, stat=istat)    

       CALL second0(toff)
       time_update_bfield = time_update_bfield + (toff-ton)

       END SUBROUTINE update_bfield

!*******************************************************************************
!>  \brief Use Faraday's law dB = -curl(E) to compute magnitic field perturbation
!*******************************************************************************
       SUBROUTINE Faraday(djbsupsmnh, djbsupumnh, djbsupvmnh, ipar)
       USE shared_data, ONLY: l_natural, delta_t 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: ipar
       REAL(dp),DIMENSION(:,:,:), ALLOCATABLE ::                         &
            djbsupsmnh, djbsupumnh, djbsupvmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER   :: fours, fouru, fourv, nblock, m, n, nmin, nmax
       REAL(dp)  :: sparity, r0, mp, np
!-----------------------------------------------
       nmin=nsmin; nmax=nsmax

       IF (ipar .EQ. isym) THEN                                            !e_s (sin), e_u (cos), e_v (cos)
          fours=f_sin; fouru=f_cos; fourv=f_cos
          r0 = hs_i/2; sparity = 1
       ELSE                                                                !e_s (cos), e_u (sin), e_v (sin)
          fours=f_cos; fouru=f_sin; fourv=f_sin
          r0 = -hs_i/2; sparity = -1
       END IF

!Calculate harmonics of electric field
       CALL tomnsp_par(esubsijf, esubsmnf, fours)
       CALL tomnsp_par(esubuijf, esubumnf, fouru)
       CALL tomnsp_par(esubvijf, esubvmnf, fourv)
       
!On exit, esubXh valid on [nsmin+1:nsmax]
       nblock = SIZE(esubsmnh,1)*SIZE(esubsmnh,2)
       CALL to_half_mesh(esubsmnf, esubsmnh)
       CALL to_half_mesh(esubumnf, esubumnh)
       CALL to_half_mesh(esubvmnf, esubvmnh)

!Gradients of esubu,v on half grid j-half=2,ns: note esubu,v(1) used here!
!On exit, work1,2 valid on [nsmin+1:nsmax]
       CALL GradientHalf(de_uds, esubumnf)
       CALL GradientHalf(de_vds, esubvmnf)


!IMPOSE NATURAL BOUNDARY CONDITION AT FIRST HALF-GRID POINT (SPH 04-28-16)
       IF (nmin.EQ.1 .AND. .NOT.l_natural) THEN
          djbsupsmnh(m1,:,1) = (esubsmnh(m1,:,2)*r0-esubumnh(m1,:,2))/2
          esubsmnh(m1,:,2) = esubsmnh(m1,:,2) - djbsupsmnh(m1,:,1)/r0
          esubumnh(m1,:,2) = esubumnh(m1,:,2) + djbsupsmnh(m1,:,1)
          de_uds(m1,:,2) = de_uds(m1,:,2) - djbsupsmnh(m1,:,1)/r0         !Preserves div B = 0
       END IF

!Calculate -curl(E) = r.h.s. of Faraday's equation for all components of jacobh*BSUPX
!Compute contravariant components of -jacobian*curl(E). Need nmax=endglobrow+1 to
!compute perturbed currents on full mesh in siesta_forces
       nmin=startglobrow
       DO m = 0, mpol
         mp = sparity*m
         DO n = -ntor, ntor
           np = sparity*n*nfp
           djbsupsmnh(m,n,nmin:nmax) = mp*esubvmnh(m,n,nmin:nmax)        & ! -mp  ~ d/du
                                     - np*esubumnh(m,n,nmin:nmax)          ! -np  ~ d/dv
           djbsupumnh(m,n,nmin:nmax) =-np*esubsmnh(m,n,nmin:nmax)        &
		                             + de_vds(m,n,nmin:nmax)
           djbsupvmnh(m,n,nmin:nmax) = mp*esubsmnh(m,n,nmin:nmax)        &
		                             - de_uds(m,n,nmin:nmax)
          ENDDO
       ENDDO

!CONSTRAINT AT ORIGIN [jbsups(m=1) = r0*jbsupu(m=1) <=> esubumnch = r0*esubsmnsh]
!Origin bc for m=1: e_u = r0*e_s (Follows from v X B at origin)
!This suffices to make B^s = r0*B^u
       IF (nmin .EQ. 1) THEN
          djbsupsmnh(:,:,1) = 0
          djbsupumnh(:,:,1) = 0
          djbsupvmnh(:,:,1) = 0
      
!      SYMMETRIC FORM (THESE AREN'T USED ANYWHERE)
          djbsupumnh(m1,:,1) = djbsupumnh(m1,:,2)
          IF (l_natural) THEN
             djbsupsmnh(m1,:,1) = r0*djbsupumnh(m1,:,1)
          ELSE
             djbsupsmnh(m1,:,1) = djbsupsmnh(m1,:,2)
          END IF

          djbsupvmnh(m0,:,1) = djbsupvmnh(m0,:,2)
       END IF

!
!      Calculate INCREMENT of magnetic field harmonics: use Euler scheme with DT given by VSUBMN advance equations.
!      Results computed here for HALF-GRID perturbations in MN space are used in the calculation of the force and
!      to advance the B's in UPDATE_STATE (in siesta_state module)
!
       djbsupsmnh = delta_t*djbsupsmnh
       djbsupumnh = delta_t*djbsupumnh
       djbsupvmnh = delta_t*djbsupvmnh

       END SUBROUTINE Faraday

       END MODULE siesta_bfield

