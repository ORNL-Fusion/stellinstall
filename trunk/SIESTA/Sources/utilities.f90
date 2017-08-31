!#if defined(SKS)
!---------------------------------------------------------------------------------------------
      SUBROUTINE GradientHalf_par(gradienth, vectorf, n2, nsmin, nsmax)
      USE stel_kinds
      USE quantities, ONLY: ns, ohs
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)      :: n2, nsmin, nsmax
      REAL(dp), INTENT(IN)  :: vectorf(1+n2*(nsmin-1):n2*nsmax)   
      REAL(dp), INTENT(OUT) :: gradienth(1+n2*(nsmin-1):n2*nsmax)      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                  :: lstart, lend
!-----------------------------------------------
!     Contraction map
!     INPUT (vectorf)   : [nsmin:nsmax]
!     OUTPUT(gradienth) : [nsmin+1:nsmax]

      lstart=1+n2*nsmin; lend=n2*nsmax
      gradienth(lstart:lend) = ohs*(vectorf(lstart:lend) - vectorf(lstart-n2:lend-n2))

      END SUBROUTINE GradientHalf_par
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE GradientFull_par(gradientf, vectorh, Mblk, nsmin, nsmax)
      USE stel_kinds
      USE quantities, ONLY: ns, ohs
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)   :: Mblk, nsmin, nsmax
      REAL(dp), INTENT(in)  :: vectorh(1+Mblk*(nsmin-1):Mblk*nsmax)
      REAL(dp), INTENT(out) :: gradientf(1+Mblk*(nsmin-1):Mblk*nsmax)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: lstart, lend
!-----------------------------------------------
!     Contraction map
!     INPUT (vectorh)   : [nsmin:nsmax]
!     OUTPUT(gradientf) : [nsmin:nsmax-1]

      lstart=1+Mblk*(nsmin-1); lend=Mblk*(nsmax-1)
      gradientf(lstart:lend) = ohs*(vectorh(lstart+Mblk:lend+Mblk) - vectorh(lstart:lend))

      IF (nsmin .EQ. 1)  gradientf(1:Mblk) = 0

      IF (nsmax.EQ.ns .AND. nsmax.NE.nsmin) THEN
         lstart=1+lend;  lend = lstart+Mblk-1
         gradientf(lstart:lend) = gradientf(lstart-Mblk:lend-Mblk)
      END IF

      END SUBROUTINE GradientFull_par
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE to_half_mesh_par(vecf, vech, n2, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ns=>ns_i
      IMPLICIT NONE        
!  
!       This subroutine migrates a full mesh quantity to the half mesh
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)   :: n2, nsmin, nsmax
      REAL(dp), INTENT(IN)  :: vecf(1+n2*(nsmin-1):n2*nsmax)   
      REAL(dp), INTENT(OUT) :: vech(1+n2*(nsmin-1):n2*nsmax)      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: lstart, lend
      REAL(dp), PARAMETER   :: p5 = 0.5_dp
!-----------------------------------------------
!     Contraction map
!     INPUT (vecf)   : [nsmin:nsmax]
!     OUTPUT(vech)   : [nsmin+1:nsmax]

      lstart=1+n2*nsmin; lend=n2*nsmax

      vech(lstart:lend) = p5*(vecf(lstart:lend) + vecf(lstart-n2:lend-n2))
      IF (nsmin .EQ. 1) vech(1:n2) = 0

      END SUBROUTINE to_half_mesh_par 
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh_par(vech, vecf, n2, nsmin, nsmax)
      USE stel_kinds
      USE stel_constants, ONLY: zero
      USE island_params, ns=>ns_i
      USE nscalingtools, ONLY: ErrorAbort
      IMPLICIT NONE        
!  
!       This subroutine migrates a half mesh quantity to the full mesh.
!       First and last point are extrapolated linearly.
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: n2, nsmin, nsmax
      REAL(dp), INTENT(IN)  :: vech(1+n2*(nsmin-1):n2*nsmax)      
      REAL(dp), INTENT(OUT) :: vecf(1+n2*(nsmin-1):n2*nsmax)   
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
      INTEGER     :: lstart, lend
!-----------------------------------------------
!     Contraction map
!     INPUT (vech)   : [nsmin:nsmax]
!     OUTPUT(vecf)   : [nsmin:nsmax-1]

      lstart=1+n2*(nsmin-1); lend=n2*(nsmax-1)
      vecf(lstart:lend) = p5*(vech(lstart:lend) + vech(lstart+n2:lend+n2))

      IF (nsmin .EQ. 1) THEN
         vecf(1:n2)  = vech(1+n2:2*n2)
         IF (nsmax .EQ. nsmin) CALL ErrorAbort('TO_FULL_MESH_PAR needs to broadcast vech(2)')
      END IF
      IF (nsmax .EQ. ns) THEN
         lstart=lend+1;  lend = lstart+n2-1
         vecf(lstart:lend) = vech(lstart:lend)
      END IF

      END SUBROUTINE to_full_mesh_par
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE toupper_forces_par
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,            &
                         mpol=>mpol_i, ntor=>ntor_i
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf,             &
                            fsubsmnsf, fsubumncf, fsubvmncf,             &
                            fsupsmnsf, fsupumncf, fsupvmncf,             &
                            fsupsmncf, fsupumnsf, fsupvmnsf, wint 
      USE fourier, ONLY: toijsp_par, tomnsp_par, f_none, f_cos, f_sin,   &
                         f_sum
      USE nscalingtools, ONLY: startglobrow, endglobrow, ErrorAbort
      USE shared_data, ONLY: lasym
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                         &
                                    fsubsijf, fsubuijf, fsubvijf,        &
                                    fsupsijf, fsupuijf, fsupvijf
!-----------------------------------------------
      INTEGER :: nsmin, nsmax, ns_span, istat
!-----------------------------------------------
      nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow)
      ns_span=nsmax-nsmin+1

      ALLOCATE(fsubsijf(ntheta,nzeta,nsmin:nsmax),                       &
               fsubuijf(ntheta,nzeta,nsmin:nsmax),                       &
               fsubvijf(ntheta,nzeta,nsmin:nsmax),                       &
               fsupsijf(ntheta,nzeta,nsmin:nsmax),                       &
               fsupuijf(ntheta,nzeta,nsmin:nsmax),                       &
               fsupvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort(' Allocation error in toupper_forces_par')

      CALL toijsp_par(fsubsmncf(:,:,nsmin:nsmax), fsubsijf, f_none, f_cos, 1, ns_span)
      CALL toijsp_par(fsubumnsf(:,:,nsmin:nsmax), fsubuijf, f_none, f_sin, 1, ns_span)
      CALL toijsp_par(fsubvmnsf(:,:,nsmin:nsmax), fsubvijf, f_none, f_sin, 1, ns_span)
 
      IF (lasym) THEN

      CALL toijsp_par(fsubsmnsf(:,:,nsmin:nsmax), fsubsijf, f_sum, f_sin, 1, ns_span)
      CALL toijsp_par(fsubumncf(:,:,nsmin:nsmax), fsubuijf, f_sum, f_cos, 1, ns_span)
      CALL toijsp_par(fsubvmncf(:,:,nsmin:nsmax), fsubvijf, f_sum, f_cos, 1, ns_span)
 
      END IF

!      
!     COMPUTE CONTRAVARIANT COMPONENTS FOR CHECKING OUTPUT
!
      CALL toupper_par(fsubsijf, fsubuijf, fsubvijf, fsupsijf, fsupuijf, fsupvijf) 
      DEALLOCATE(fsubsijf, fsubuijf, fsubvijf)

      CALL tomnsp_par(fsupsijf, fsupsmncf(:,:,nsmin:nsmax), f_cos, 1, ns_span)
      CALL tomnsp_par(fsupuijf, fsupumnsf(:,:,nsmin:nsmax), f_sin, 1, ns_span)
      CALL tomnsp_par(fsupvijf, fsupvmnsf(:,:,nsmin:nsmax), f_sin, 1, ns_span)

      IF (lasym) THEN

      CALL tomnsp_par(fsupsijf, fsupsmnsf(:,:,nsmin:nsmax), f_sin, 1, ns_span)
      CALL tomnsp_par(fsupuijf, fsupumncf(:,:,nsmin:nsmax), f_cos, 1, ns_span)
      CALL tomnsp_par(fsupvijf, fsupvmncf(:,:,nsmin:nsmax), f_cos, 1, ns_span)

      END IF

      DEALLOCATE(fsupsijf, fsupuijf, fsupvijf, stat=istat)

      END SUBROUTINE toupper_forces_par
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE tolowerh_par(xsupsij, xsupuij, xsupvij,                 &
                              xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: nuv=>nuv_i
      USE metrics, ONLY: gss, gsu, gsv, guu, guv, gvv                     ! LOWER METRIC ELEMENTS (ON half mesh)
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

      END SUBROUTINE tolowerh_par
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE tolowerf_par(xsupsij, xsupuij, xsupvij,                 &
                              xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: ns=>ns_i, nuv=>nuv_i
      USE metrics, ONLY: gssf, guuf, gvvf, gsuf, gsvf, guvf
      IMPLICIT NONE        
!
!     This subroutine moves from lower to upper magnetic field components.
!     on the full mesh. To avoid interpolating, inverts directly the upper metric elements
!
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

      END SUBROUTINE tolowerf_par
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
        SUBROUTINE toupper_par(xsubsij, xsubuij, xsubvij,                &
                               xsupsij, xsupuij, xsupvij)
        USE stel_kinds
        USE island_params, ns=>ns_i, nuv=>nuv_i
        USE metrics, ONLY: hss, hsu, hsv, huu, huv, hvv                   ! UPPER METRIC ELEMENTS (ON full mesh)
        USE nscalingtools, ONLY: startglobrow, endglobrow
        IMPLICIT NONE        
!
!     This subroutine converts from lower to upper components on the FULL mesh
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER :: istat
        REAL(dp), DIMENSION(nuv,ns), INTENT(IN)  ::                      &
        xsubsij, xsubuij, xsubvij 
        REAL(dp), DIMENSION(nuv,ns), INTENT(OUT) ::                      &
        xsupsij, xsupuij, xsupvij        
!-----------------------------------------------
        INTEGER :: nsmin, nsmax, js, js1
!-----------------------------------------------

        nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow)
        DO js=nsmin, nsmax
           js1 = js+1-nsmin
        xsupsij(:,js1) = hss(:,js)*xsubsij(:,js1)                        &
                       + hsu(:,js)*xsubuij(:,js1)                        &
                       + hsv(:,js)*xsubvij(:,js1)
        xsupuij(:,js1) = hsu(:,js)*xsubsij(:,js1)                        &
                       + huu(:,js)*xsubuij(:,js1)                        &
                       + huv(:,js)*xsubvij(:,js1)
        xsupvij(:,js1) = hsv(:,js)*xsubsij(:,js1)                        &
                       + huv(:,js)*xsubuij(:,js1)                        &
                       + hvv(:,js)*xsubvij(:,js1)       
        END DO

        END SUBROUTINE toupper_par
!---------------------------------------------------------------------------------------------
!#endif
        SUBROUTINE to_half_mesh(vecf, vech, n2)
        USE stel_kinds
        USE island_params, ns=>ns_i
        IMPLICIT NONE        
!  
!       This subroutine migrates a full mesh quantity to the half mesh.
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(IN)                        :: n2
        REAL(dp), DIMENSION(n2,ns), INTENT(IN)  :: vecf    
        REAL(dp), DIMENSION(n2,ns), INTENT(OUT) :: vech  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        REAL(dp)  :: p5 = 0.5_dp
!-----------------------------------------------
     
        vech(:,1) = 0
        vech(:,2:ns) = p5*(vecf(:,2:ns) + vecf(:,1:nsh))
        
        END SUBROUTINE to_half_mesh        


        SUBROUTINE to_full_mesh(vech, vecf, n2)
#if defined(SKS)
        USE timer_mod, ONLY: to_full_mesh_time
        USE descriptor_mod, ONLY: DIAGONALDONE, INHESSIAN
#endif        
        USE stel_kinds
        USE stel_constants, ONLY: zero
        USE island_params, ns=>ns_i, mpol=>mpol_i, ntor=>ntor_i
        IMPLICIT NONE        
!  
!       This subroutine migrates a half mesh quantity to the full mesh.
!       First and last point are extrapolated linearly.
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(IN)                     :: n2
        REAL(dp), DIMENSION(n2,ns), INTENT(IN)  :: vech      
        REAL(dp), DIMENSION(n2,ns), INTENT(OUT) :: vecf   
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER     :: i1=2, i2=i1+1
        REAL(dp), PARAMETER :: p5 = 0.5_dp
        REAL(dp) :: skston, skstoff
!-----------------------------------------------
        CALL second0(skston)

        vecf(:,i1:nsh) = p5*(vech(:,i1:nsh) + vech(:,i2:ns))
!        vecf(ns,:) = 2*vech(ns,:) - vecf(nsh,:)
!        vecf(1,:)  = 2*vech(2,:) - vecf(2,:)
        vecf(:,ns) = vech(:,ns)
        vecf(:,1) = vech(:,2)
                           
        CALL second0(skstoff)
#if defined(SKS)
        IF(DIAGONALDONE.AND.INHESSIAN) to_full_mesh_time=to_full_mesh_time+(skstoff-skston)
#endif
        END SUBROUTINE to_full_mesh

        
        SUBROUTINE toupper(xsubsij, xsubuij, xsubvij, xsupsij,           &
                           xsupuij, xsupvij)
        USE stel_kinds
        USE island_params, ns=>ns_i, nuv=>nuv_i
        USE metrics, ONLY: hss, hsu, hsv, huu, huv, hvv                   ! UPPER METRIC ELEMENTS (ON full mesh)
        IMPLICIT NONE        
!
!     This subroutine moves from lower to upper components on the FULL mesh
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER :: istat
        REAL(dp), DIMENSION(nuv,ns), INTENT(IN)  ::                      &
          xsubsij, xsubuij, xsubvij
        REAL(dp), DIMENSION(nuv,ns), INTENT(OUT) ::                      &
          xsupsij, xsupuij, xsupvij        
!-----------------------------------------------

        xsupsij = hss*xsubsij + hsu*xsubuij + hsv*xsubvij
        xsupuij = hsu*xsubsij + huu*xsubuij + huv*xsubvij
        xsupvij = hsv*xsubsij + huv*xsubuij + hvv*xsubvij       
        
        END SUBROUTINE toupper

      SUBROUTINE tolowerh(xsupsij, xsupuij, xsupvij, xsubsij,            &
                          xsubuij, xsubvij)
      USE stel_kinds
      USE island_params, ns=>ns_i, nuv=>nuv_i
      USE metrics, ONLY: gss, gsu, gsv, guu, guv, gvv                     ! LOWER METRIC ELEMENTS (ON half mesh)
      IMPLICIT NONE        
!
!     Moves contravariant (upper) to covariant (lower) magnetic field
!     components on the half radial mesh.
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!----------------------------------------------
      INTEGER :: istat
      REAL(dp), DIMENSION(nuv,ns), INTENT(IN)  ::                        &
          xsupsij, xsupuij, xsupvij 
      REAL(dp), DIMENSION(nuv,ns), INTENT(OUT) ::                        &
          xsubsij, xsubuij, xsubvij
!-----------------------------------------------

      xsubsij = gss*xsupsij + gsu*xsupuij + gsv*xsupvij
      xsubuij = gsu*xsupsij + guu*xsupuij + guv*xsupvij
      xsubvij = gsv*xsupsij + guv*xsupuij + gvv*xsupvij       
        
      END SUBROUTINE tolowerh


      SUBROUTINE tolowerf(xsupsij, xsupuij, xsupvij, xsubsij,            &
                          xsubuij, xsubvij)
      USE stel_kinds
      USE island_params, ns=>ns_i, nuv=>nuv_i
      USE metrics, ONLY: gssf, guuf, gvvf, gsuf, gsvf, guvf
      IMPLICIT NONE        
!
!     This subroutine moves from lower to upper magnetic field components.
!     on the full mesh. To avoid interpolating, inverts directly the upper metric elements
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: istat
      REAL(dp), DIMENSION(nuv,ns), INTENT(IN)  ::                        &
          xsupsij, xsupuij, xsupvij 
      REAL(dp), DIMENSION(nuv,ns), INTENT(OUT) ::                        &
          xsubsij, xsubuij, xsubvij
!-----------------------------------------------

      xsubsij = gssf*xsupsij + gsuf*xsupuij + gsvf*xsupvij
      xsubuij = gsuf*xsupsij + guuf*xsupuij + guvf*xsupvij
      xsubvij = gsvf*xsupsij + guvf*xsupuij + gvvf*xsupvij       
        
      END SUBROUTINE tolowerf


      SUBROUTINE toupper_forces

      CALL toupper_forces_par

      END SUBROUTINE toupper_forces

      SUBROUTINE GradientHalf(gradienth, vectorf, n2)
      USE stel_kinds
      USE quantities, ONLY: ns, nsh, ohs
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: n2
      REAL(dp), INTENT(IN)  :: vectorf(n2,ns)
      REAL(dp), INTENT(OUT) :: gradienth(n2,ns)
!-----------------------------------------------

      gradienth(:,2:ns) = ohs*(vectorf(:,2:ns) - vectorf(:,1:nsh))  

      END SUBROUTINE GradientHalf

      SUBROUTINE GradientFull(gradientf, vectorh, n2)
      USE stel_kinds
      USE quantities, ONLY: ns, nsh, ohs
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)      :: n2
      REAL(dp), INTENT(in)  :: vectorh(n2,ns)
      REAL(dp), INTENT(out) :: gradientf(n2,ns)
!-----------------------------------------------

      gradientf(:,2:nsh) = ohs*(vectorh(:,3:ns) - vectorh(:,2:nsh))  
      gradientf(:,ns) = gradientf(:,nsh)
!!    SPECIAL m values will be assigned OUTSIDE this call
      gradientf(:,1)  = 0                

      END SUBROUTINE GradientFull


      SUBROUTINE ZeroHalfPoint(array, work, jspoint, mcut, iparity)
      USE stel_kinds
      USE quantities, ONLY: ns, ntheta, nzeta, mpol, ntor, jacobh
      USE fourier, ONLY: toijsp2, tomnsp2
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: jspoint, iparity, mcut
      REAL(dp), INTENT(INOUT) :: array(0:mpol,-ntor:ntor,ns)
      REAL(dp)                :: work(ntheta,nzeta,ns)
      REAL(dp), ALLOCATABLE   :: temp(:,:)
!-----------------------------------------------
      ALLOCATE (temp(0:mpol,-ntor:ntor))

      temp = array(:,:,jspoint)       
      CALL toijsp2(temp, work, iparity, jspoint)
      work(jspoint,:,:) = work(:,:,jspoint)/jacobh(:,:,jspoint)
      CALL tomnsp2(work, temp, iparity, jspoint)
      temp(mcut:,:) = 0
      CALL toijsp2(temp, work, iparity, jspoint)
      work(jspoint,:,:) = work(:,:,jspoint)*jacobh(:,:,jspoint)
      CALL tomnsp2(work, temp, iparity, jspoint)
      array(:,:,jspoint) = temp

      DEALLOCATE (temp)

      END SUBROUTINE ZeroHalfPoint


      SUBROUTINE SurfAverage(average, q3d, nsmin, nsmax)
      USE stel_kinds, ONLY: dp
	  USE quantities, ONLY: wint
      USE island_params, ONLY: ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: nsmin, nsmax
      REAL(dp), INTENT(OUT) :: average(nsmin:nsmax)
      REAL(dp), INTENT(IN)  :: q3d(ntheta,nzeta,nsmin:nsmax)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER    :: m0=0
      INTEGER               :: js
!-----------------------------------------------

      DO js = nsmin, nsmax
		 average(js) = SUM(q3d(:,:,js)*wint(:,:,js))
      END DO

      END SUBROUTINE SurfAverage

#if defined(MRC)
#else
      SUBROUTINE Truncate(num, iprec0)
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)     :: iprec0
      REAL(dp), INTENT(INOUT) :: num
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*24 :: chnum, tchnum
!-----------------------------------------------
!
!     TRUNCATES double-precision to precision iprec0 digits, keeping exponent range of double
!     WRITE TO INTERNAL FILE TO DO TRUNCATION
!
!     RETURN

      WRITE (chnum, '(a,i2,a,i2,a)') '(1p,e',iprec0+7,'.',iprec0,')'
      WRITE (tchnum, chnum) num

      READ (tchnum, chnum) num

      END SUBROUTINE Truncate
#endif

