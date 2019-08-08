!>  \brief Module contained subroutines for updating half-grid magnetic fields as part of the SIESTA project
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 21, 2006
       MODULE siesta_bfield
       USE v3_utilities, ONLY: assert
       USE stel_kinds
       USE stel_constants
       USE fourier 
       USE nscalingtools, ONLY: startglobrow, endglobrow
       USE timer_mod
       IMPLICIT NONE
       
       PRIVATE
       INTEGER            :: nsmin, nsmax

       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                       &
         esubsijf, esubuijf, esubvijf
       REAL(dp), DIMENSION(:,:,:), POINTER     ::                       &
         esubsmnf, esubumnf, esubvmnf

       PUBLIC :: update_bfield

       CONTAINS
!*******************************************************************************
!>  \brief Routine (both parallel, serial) for updating B^s, B^u, B^v 
!>  for both ideal and resistive perturbations.
!>
!>  Advances half-mesh magnetic field components in nm space from time t to
!>  t + dt. Used to update values of BSUPMNH and BSUPIJH stored in \ref quantities
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
       USE siesta_namelist, ONLY: lresistive, eta_factor
       USE shared_data, ONLY: lasym, fsq_res, fsq_total, buv_res
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(IN)     :: l_add_res
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: istat, js
       LOGICAL  :: ladd_pert0, ltest=.TRUE.
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: resistivity
       REAL(dp) :: ton, toff, rho
       REAL(dp) :: delt_cfl, eta_prof
 !-----------------------------------------------
       CALL second0(ton)

       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+1)

!
!Calculate covariant components of (non-resistive) electric field E = -v X B on full-space mesh
!  
       ALLOCATE(esubsijf(ntheta,nzeta,nsmin:nsmax),                     &
                esubuijf(ntheta,nzeta,nsmin:nsmax),                     &
                esubvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)   
       CALL ASSERT(istat.EQ.0,'Allocation1 failed in UPDATE_BFIELD')


       esubsijf =-(jvsupuijf(:,:,nsmin:nsmax)*bsupvijf0(:,:,nsmin:nsmax) &
                 - jvsupvijf(:,:,nsmin:nsmax)*bsupuijf0(:,:,nsmin:nsmax))  ! -E_s = V^uB^v - V^vB^u  (SINE): jv -> jac*v
       esubuijf =-(jvsupvijf(:,:,nsmin:nsmax)*bsupsijf0(:,:,nsmin:nsmax) &
                 - jvsupsijf(:,:,nsmin:nsmax)*bsupvijf0(:,:,nsmin:nsmax))  ! -E_u = V^vB^s - V^sB^v  (COSINE)
       esubvijf =-(jvsupsijf(:,:,nsmin:nsmax)*bsupuijf0(:,:,nsmin:nsmax) &
                 - jvsupuijf(:,:,nsmin:nsmax)*bsupsijf0(:,:,nsmin:nsmax))  ! -E_v = V^sB^u - V^uB^s  (COSINE)

!verify boundary condition: esubu,v(s=1) = 0 (tangential E vanishes at bdy => dB^s = 0)
       IF (nsmax .EQ. ns) ltest = ALL(esubuijf(:,:,ns).EQ.zero)
       CALL ASSERT(ltest,'esubuijf(ns) != 0 in UPDATE_BFIELD')
       IF (nsmax .EQ. ns) ltest = ALL(esubvijf(:,:,ns).EQ.zero)
       CALL ASSERT(ltest,'esubvijf(ns) != 0 in UPDATE_BFIELD')

!
!Note this will LOWER the energy due to eta*|J|||**2 heating
!      esubX(resistive) = eta(JdotB/B^2)*BsubX
!so the magnetic energy DECREASES due to this term. Note ksubX=JsubX are the
!covariant components of the current
!
       ladd_pert0 = ALLOCATED(buv_res)
       IF (lresistive .AND. (l_add_res .or. ladd_pert0)) THEN
          delt_cfl = hs_i**2*ABS(eta_factor)
          IF (fsq_total < fsq_res) delt_cfl = delt_cfl*SQRT(fsq_total/fsq_res)

          ALLOCATE(resistivity(ntheta,nzeta,nsmin:nsmax), stat=istat)
          CALL ASSERT(istat.EQ.0,'Allocation2 failed in UPDATE_BFIELD')

          DO js = nsmin, nsmax
             rho = hs_i*(js-1)
             eta_prof = rho*rho*(1-rho)
             resistivity(:,:,js) = delt_cfl*eta_prof
          END DO

          IF (ladd_pert0) THEN
             resistivity = resistivity*buv_res(:,:,nsmin:nsmax)
          ELSE
             resistivity = resistivity/jacobf(:,:,nsmin:nsmax)            !divide out jacobf factor from cv_currents
          ENDIF

!Isotropic resistivity, E ~ eta*J; for ladd_pert0=TRUE, K ~ B in init_state
          esubsijf = esubsijf + resistivity*ksubsijf(:,:,nsmin:nsmax)
          esubuijf = esubuijf + resistivity*ksubuijf(:,:,nsmin:nsmax)
          esubvijf = esubvijf + resistivity*ksubvijf(:,:,nsmin:nsmax)

          DEALLOCATE (resistivity)
       ENDIF 
        

!Allocate working arrays used in Faraday
       ALLOCATE(esubsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                esubumnf(0:mpol,-ntor:ntor,nsmin:nsmax),                & 
                esubvmnf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       CALL ASSERT(istat.EQ.0,'Allocation3 failed in UPDATE_BFIELD')

!Update Bfield using Faraday's Law
       CALL Faraday(djbsupsmnsh, djbsupumnch, djbsupvmnch, isym)
       IF (lasym)                                                       &
       CALL Faraday(djbsupsmnch, djbsupumnsh, djbsupvmnsh, iasym)

       DEALLOCATE (esubsijf, esubuijf, esubvijf,                        &
                   esubsmnf, esubumnf, esubvmnf, stat=istat)    
       CALL ASSERT(istat.EQ.0,'Deallocation failed in UPDATE_BFIELD')

       CALL second0(toff)
       time_update_bfield = time_update_bfield + (toff-ton)

       END SUBROUTINE update_bfield
       

!*******************************************************************************
!>  \brief Use Faraday's law dB = -curl(E) to compute magnitic field perturbation
!*******************************************************************************
       SUBROUTINE Faraday(djbsupsmnh, djbsupumnh, djbsupvmnh, ipar)
       USE shared_data, ONLY: l_natural, delta_t 
       USE hessian, ONLY: l_compute_hessian
       USE utilities, ONLY: Curl_FtoH
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: ipar
       REAL(dp),DIMENSION(:,:,:), ALLOCATABLE ::                        &
            djbsupsmnh, djbsupumnh, djbsupvmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER   :: fours, fouru, fourv, m, n, nmin, nmax
       REAL(dp)  :: sparity, r12, mp, np
!-----------------------------------------------
       nmin=nsmin; nmax=nsmax

       IF (ipar .EQ. isym) THEN                                           !e_s (sin), e_u (cos), e_v (cos)
          fours=f_sin; fouru=f_cos; fourv=f_cos; sparity = 1
       ELSE                                                               !e_s (cos), e_u (sin), e_v (sin)
          fours=f_cos; fouru=f_sin; fourv=f_sin; sparity = -1
       END IF
       r12 = sparity*hs_i/2

!Calculate harmonics of electric field
       CALL tomnsp_par(esubsijf, esubsmnf, fours)
       CALL tomnsp_par(esubuijf, esubumnf, fouru)
       CALL tomnsp_par(esubvijf, esubvmnf, fourv)
       
!IMPOSE BOUNDARY CONDITION AT FIRST HALF-GRID POINT (SPH 04-28-16)
       IF (startglobrow .EQ. 1) THEN
!These conditions are needed so that delta-W jogs (CheckForces) agree with forces at origin
          esubsmnf(m2:,:,1) = 0                                          !f_v(1) wrong w/0 this
          esubumnf(m0,:,1) = 0; esubumnf(m2:,:,1) = 0                    !f_v(1) wrong w/0 this
          esubvmnf(m1:,:,1) = 0
          
          IF (.NOT.l_natural) THEN
             djbsupsmnh(m1,:,2) = (esubsmnf(m1,:,1) + esubsmnf(m1,:,2))/2
             djbsupumnh(m1,:,2) = (esubumnf(m1,:,1) + esubumnf(m1,:,2))/2
             djbsupsmnh(m1,:,1) = r12*djbsupsmnh(m1,:,2) - djbsupumnh(m1,:,2)          !-> 0
             !This constrains [esubs(2)*r12 - esubu(2)] = 0 at s=r12 (first half grid pt)
             !needed to make djb^s ~ r12*djb^u there
             esubsmnf(m1,:,1) = esubsmnf(m1,:,1) - djbsupsmnh(m1,:,1)/r12
             esubumnf(m1,:,1) = esubumnf(m1,:,1) + djbsupsmnh(m1,:,1)
          END IF
       END IF
       
       CALL CURL_FtoH( esubsmnf,   esubumnf,   esubvmnf,                &
                     djbsupsmnh, djbsupumnh, djbsupvmnh, isym, nmin, nmax)

!BOUNDARY CONDITIONS AT ORIGIN (used to compute bfields at origin in siesta_init, siesta_forces)
       IF (nmin .EQ. 1) THEN
          djbsupsmnh(:,:,1) = 0; djbsupsmnh(m1,:,1) = djbsupsmnh(m1,:,2)
          djbsupumnh(:,:,1) = 0; djbsupumnh(m0:m1,:,1) = djbsupumnh(m0:m1,:,2)
          djbsupvmnh(:,:,1) = 0; djbsupvmnh(m0,:,1) = djbsupvmnh(m0,:,2);
       END IF
!
!      Calculate INCREMENT of magnetic field harmonics: use Euler scheme with DT given by VSUBMN advance equations.
!      Results computed here for HALF-GRID perturbations in MN space are used in the calculation of the force and
!      to advance the B's in UPDATE_STATE (in siesta_state module)
!      Minus sign since Faraday is dB = -CURL(E)
!

       djbsupsmnh = -delta_t*djbsupsmnh
       djbsupumnh = -delta_t*djbsupumnh
       djbsupvmnh = -delta_t*djbsupvmnh
       
       END SUBROUTINE Faraday

       END MODULE siesta_bfield


