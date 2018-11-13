      MODULE utilities
      USE stel_kinds
      USE island_params, ONLY: ohs=>ohs_i, ns=>ns_i, mpol=>mpol_i,      &
                               ntor=>ntor_i, nfp=>nfp_i
      USE v3_utilities, ONLY: assert
      IMPLICIT NONE

      INTEGER, PARAMETER :: m0=0, m1=1, m2=2, isym=0, iasym=1
     
      INTERFACE to_full_mesh
         MODULE PROCEDURE to_full_mesh2, to_full_mesh3
      END INTERFACE
      
      INTERFACE CURL_HtoF
         MODULE PROCEDURE CURL_HTOF, CURL_HtoF_TARG
      END INTERFACE

#ifndef WIN32
      INTERFACE to_half_mesh
         MODULE PROCEDURE to_half_mesh, to_half_mesh_3P
      END INTERFACE

      INTERFACE GradientHalf
         MODULE PROCEDURE GradientHalf, GradientHalf3P
      END INTERFACE
#endif

      REAL(dp), PARAMETER, PRIVATE :: p5 = 0.5_dp
      INTEGER, PRIVATE             :: nsmin, nsmax

      CONTAINS
      
      SUBROUTINE GradientHalf(gradienth, vecf)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), ALLOCATABLE :: vecf(:,:,:)
      REAL(dp), ALLOCATABLE :: gradienth(:,:,:)
!-----------------------------------------------
!     Contraction map
!     INPUT (vecf)   : [nsmin:nsmax]
!     OUTPUT(gradienth) : [nsmin+1:nsmax]
      nsmin = LBOUND(gradienth,3); nsmax = UBOUND(gradienth,3)

      gradienth(:,:,nsmin+1:nsmax) = ohs*(vecf(:,:,nsmin+1:nsmax)    & 
                                   -      vecf(:,:,nsmin:nsmax-1))
      gradienth(:,:,nsmin) = 0

      END SUBROUTINE GradientHalf
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
#ifndef WIN32      
	  SUBROUTINE GradientHalf3P(gradienth, vecf)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), POINTER     :: vecf(:,:,:)
      REAL(dp), ALLOCATABLE :: gradienth(:,:,:)
!-----------------------------------------------
!     Contraction map
!     INPUT (vecf)      : [nsmin:nsmax]
!     OUTPUT(gradienth) : [nsmin+1:nsmax]

      nsmin = LBOUND(gradienth,3); nsmax = UBOUND(gradienth,3)

      CALL ASSERT(nsmin+1.GE.LBOUND(vecf,3),'LBOUND WRONG IN GradientHalf3P')
      CALL ASSERT(nsmax.LE.UBOUND(vecf,3),'UBOUND WRONG IN GradientHalf3P')

      gradienth(:,:,nsmin+1:nsmax) = ohs*(vecf(:,:,nsmin+1:nsmax)       & 
                                   -      vecf(:,:,nsmin:nsmax-1))
      gradienth(:,:,nsmin) = 0

      END SUBROUTINE GradientHalf3P
#endif
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE GradientFull(gradientf, vech)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), ALLOCATABLE  :: vech(:,:,:)
      REAL(dp), ALLOCATABLE :: gradientf(:,:,:)
!-----------------------------------------------
!     Contraction map
!     INPUT (vech)   : [nsmin:nsmax]
!     OUTPUT(gradientf) : [nsmin:nsmax-1]
      nsmin = LBOUND(vech,3); nsmax = UBOUND(vech,3)

      gradientf(:,:,nsmin:nsmax-1) = ohs*(vech(:,:,nsmin+1:nsmax)    &
                                   -      vech(:,:,nsmin:nsmax-1))

      IF (nsmin .EQ. 1) gradientf(:,:,nsmin) = 0
!      gradientf(:,:,nsmax) = 0                              

      IF (nsmax.EQ.ns .AND. nsmax.NE.nsmin) THEN
         gradientf(:,:,nsmax) = gradientf(:,:,nsmax-1)
      END IF

      END SUBROUTINE GradientFull
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE to_half_mesh(vecf, vech)
!  
!     This subroutine migrates a full mesh quantity to the half mesh
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), ALLOCATABLE  :: vecf(:,:,:)   
      REAL(dp), ALLOCATABLE :: vech(:,:,:)
!-----------------------------------------------
!     Contraction map
!     INPUT (vecf)   : [nsmin:nsmax]
!     OUTPUT(vech)   : [nsmin+1:nsmax]

      nsmin = LBOUND(vech,3); nsmax = UBOUND(vech,3)

      vech(:,:,nsmin+1:nsmax) = p5*(vecf(:,:,nsmin:nsmax-1)             &
                              +     vecf(:,:,nsmin+1:nsmax))
      vech(:,:,nsmin) = 0

      END SUBROUTINE to_half_mesh 
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
#ifndef WIN32
      SUBROUTINE to_half_mesh_3P(vecf, vech)
!  
!     This subroutine migrates a full mesh quantity to the half mesh
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), POINTER  :: vecf(:,:,:)   
      REAL(dp), ALLOCATABLE :: vech(:,:,:)
!-----------------------------------------------
!     Contraction map
!     INPUT (vecf)   : [nsmin:nsmax]
!     OUTPUT(vech)   : [nsmin+1:nsmax]

      nsmin = LBOUND(vech,3); nsmax = UBOUND(vech,3)
      CALL ASSERT(nsmin+1.GE.LBOUND(vecf,3),'LBOUND WRONG IN to_half_mesh_3P')
      CALL ASSERT(nsmax.LE.UBOUND(vecf,3),'UBOUND WRONG IN to_half_mesh_3P')

      vech(:,:,nsmin+1:nsmax) = p5*(vecf(:,:,nsmin+1:nsmax)             &
                              +     vecf(:,:,nsmin:nsmax-1))
!      vech(:,:,nsmin) = 0

      END SUBROUTINE to_half_mesh_3P
#endif
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh3(vech, vecf)
!  
!     This subroutine migrates a half mesh quantity to the full mesh.
!     First and last point are extrapolated linearly.
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), ALLOCATABLE :: vech(:,:,:)
      REAL(dp), ALLOCATABLE :: vecf(:,:,:)
!-----------------------------------------------
!     Contraction map
!     INPUT (vech)   : [nsmin:nsmax]
!     OUTPUT(vecf)   : [nsmin:nsmax-1]
 
      nsmin = LBOUND(vecf,3); nsmax = UBOUND(vecf,3)
      CALL ASSERT(nsmin.GE.LBOUND(vech,3),'LBOUND WRONG IN to_full_mesh3')
      CALL ASSERT(nsmax.LE.UBOUND(vech,3),'UBOUND WRONG IN to_full_mesh3')

      vecf(:,:,nsmin:nsmax-1) = p5*(vech(:,:,nsmin+1:nsmax)             &
                              +     vech(:,:,nsmin:nsmax-1))
      vecf(:,:,nsmax) = 0                              
      IF (nsmin .EQ. 1)  vecf(:,:,1) = vech(:,:,2)
      IF (nsmax .EQ. ns) vecf(:,:,ns)= vech(:,:,ns)
      
      END SUBROUTINE to_full_mesh3
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh2(vech, vecf)
!  
!     This subroutine migrates a half mesh quantity to the full mesh.
!     First and last point are extrapolated linearly.
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), ALLOCATABLE :: vech(:,:)
      REAL(dp), ALLOCATABLE :: vecf(:,:)
!-----------------------------------------------
!     Contraction map
!     INPUT (vech)   : [nsmin:nsmax]
!     OUTPUT(vecf)   : [nsmin:nsmax-1]

      nsmin = LBOUND(vecf,2); nsmax = UBOUND(vecf,2)
      CALL ASSERT(nsmin.GE.LBOUND(vech,2),'LBOUND WRONG IN to_full_mesh2')
      CALL ASSERT(nsmax.LE.UBOUND(vech,2),'UBOUND WRONG IN to_full_mesh2')

      vecf(:,nsmin:nsmax-1) = p5*(vech(:,nsmin:nsmax-1)                 &
                            +     vech(:,nsmin+1:nsmax))

      vecf(:,nsmax) = 0                              
      IF (nsmin .EQ. 1)  vecf(:,1) = vech(:,2)
      IF (nsmax .EQ. ns) vecf(:,ns)= vech(:,ns)

      END SUBROUTINE to_full_mesh2

!>    \brief Compute real-space components of contravariant jac*B^i 
!>     on half radial grid from curl(A) (drop "j" on B's for notational
!>     simplicity)
      SUBROUTINE CURL_FtoH(Asubsmnf, Asubumnf, Asubvmnf,                &
                           Bsupsmnh, Bsupumnh, Bsupvmnh, ipar, nsmin0, nsmax0)
      USE nscalingtools, ONLY: startglobrow
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ipar, nsmin0, nsmax0
      REAL(dp), POINTER, DIMENSION(:,:,:)     ::                        &
                Asubsmnf, Asubumnf, Asubvmnf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
                bsupsmnh, bsupumnh, bsupvmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: m, n, mp, np, sparity, istat, moff, nmin, nmax
!----------------------------------------------- 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
         asubsmnh, asubumnh, asubvmnh, da_uds, da_vds
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         sparity=1
      ELSE
         sparity=-1
      END IF

      nmin = nsmin0; nmax = nsmax0

      ALLOCATE(asubsmnh(0:mpol,-ntor:ntor,nmin:nmax),                   &
               asubumnh(0:mpol,-ntor:ntor,nmin:nmax),                   &
               asubvmnh(0:mpol,-ntor:ntor,nmin:nmax), stat=istat)   
      CALL ASSERT(istat.EQ.0,'Allocation1 failed in CURL_FtoH')

      CALL to_half_mesh(asubsmnf, asubsmnh)                       
      CALL to_half_mesh(asubumnf, asubumnh)
      CALL to_half_mesh(asubvmnf, asubvmnh)

!
!     Compute sqrt(g)B^X contravariant B-field components in Fourier space
!                   
      ALLOCATE(da_uds(0:mpol,-ntor:ntor,nmin:nmax),                     &
               da_vds(0:mpol,-ntor:ntor,nmin:nmax), stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation2 failed in CURL_FtoH')

!     Radial gradients of asubmnu,v on half grid (excluding origin and edge: treat separately)
      CALL GradientHalf(da_uds, asubumnf)
      CALL GradientHalf(da_vds, asubvmnf)

!
!     (sqrt(g)*B^X)mn
!
      nmin = MAX(2, startglobrow)
      CALL ASSERT(nmin.GE.LBOUND(bsupsmnh,3),'LBOUND WRONG IN CURL_FTOH')
      CALL ASSERT(nmax.LE.UBOUND(bsupsmnh,3),'UBOUND WRONG IN CURL_FTOH')
      
      DO m = 0, mpol
         mp = sparity*m
         DO n = -ntor, ntor
           np = sparity*n*nfp
           bsupsmnh(m,n,nmin:nmax) = -mp*asubvmnh(m,n,nmin:nmax)        &
                                   +  np*asubumnh(m,n,nmin:nmax)
           bsupumnh(m,n,nmin:nmax) =  np*asubsmnh(m,n,nmin:nmax)        &
                                   -  da_vds(m,n,nmin:nmax)                     
           bsupvmnh(m,n,nmin:nmax) = -mp*asubsmnh(m,n,nmin:nmax)        &
                                   +  da_uds(m,n,nmin:nmax)
         ENDDO
      ENDDO

      DEALLOCATE (asubsmnh, asubumnh, asubvmnh, da_vds, da_uds, stat=istat)   
      CALL ASSERT(istat.EQ.0,'Deallocation failed in CURL_FtoH')

      IF (startglobrow .EQ. 1) THEN
         bsupsmnh(:,:,1) = 0
         bsupumnh(:,:,1) = 0
         bsupvmnh(:,:,1) = 0
      END IF

      bsupsmnh(m0,-ntor:-1,nmin:nmax) = 0
      bsupumnh(m0,-ntor:-1,nmin:nmax) = 0
      bsupvmnh(m0,-ntor:-1,nmin:nmax) = 0

      END SUBROUTINE CURL_FtoH


!>    \brief Compute real-space components of contravariant jac*J^i == K^i 
!>     on half radial grid from curl(B)
!>     simplicity)
      SUBROUTINE CURL_HtoF(Bsubsmnh, Bsubumnh, Bsubvmnh,                &
                           Ksupsmnf, Ksupumnf, Ksupvmnf, ipar, nsmin0, nsmax0)
      USE nscalingtools, ONLY: endglobrow
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ipar, nsmin0, nsmax0
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
                Bsubsmnh, Bsubumnh, Bsubvmnh
      REAL(dp), DIMENSION(:,:,:), POINTER ::                            &
                Ksupsmnf, Ksupumnf, Ksupvmnf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: m, n, mp, np, mo, no, sparity, istat, nmin, nmax, moff, noff
!----------------------------------------------- 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
         bsubsmnf, bsubumnf, bsubvmnf, db_uds, db_vds
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         sparity=1
      ELSE
         sparity=-1
      END IF
      
      nmin = nsmin0; nmax = nsmax0
      moff = LBOUND(Ksupsmnf,1); noff = LBOUND(Ksupsmnf,2)+ntor

      ALLOCATE(bsubsmnf(0:mpol,-ntor:ntor,nmin:nmax),                   & 
               bsubumnf(0:mpol,-ntor:ntor,nmin:nmax),                   &
               bsubvmnf(0:mpol,-ntor:ntor,nmin:nmax), stat=istat)   
      CALL ASSERT(istat.EQ.0,'Allocation1 failed in CURL_HtoF')
 
      CALL to_full_mesh(bsubsmnh, bsubsmnf)                             ! Compute BSUBXMNF
      CALL to_full_mesh(bsubumnh, bsubumnf)
      CALL to_full_mesh(bsubvmnh, bsubvmnf)

!     Compute sqrt(g)J^X == K^X contravariant current components in Fourier space
!                   
      ALLOCATE(db_uds(0:mpol,-ntor:ntor,nmin:nmax),                     &
               db_vds(0:mpol,-ntor:ntor,nmin:nmax), stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation2 failed in CURL_HtoF')

!     Radial gradients of bsubmnu,v on full grid (excluding origin and edge: treat separately)
      CALL GradientFull(db_uds, bsubumnh)
      CALL GradientFull(db_vds, bsubvmnh)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
      IF (nmin .EQ. 1) THEN
         bsubsmnf(m2:,:,1) = 0
         bsubumnf(:,:,1) = 0;  bsubumnf(m1,:,1) = bsubumnh(m1,:,2)
         db_uds(:,:,1) = 0;  db_uds(m0:m1,:,1) = 2*ohs*bsubumnh(m0:m1,:,2)
         db_vds(:,:,1) = 0;  db_vds(m1,:,1) = 2*ohs*bsubvmnh(m1,:,2)
         bsubvmnf(m2:,:,1) = 0  
      END IF

      nmax = endglobrow
      CALL ASSERT(nmin.GE.LBOUND(ksupsmnf,3),'LBOUND WRONG IN CURL_HTOF')
      CALL ASSERT(nmax.LE.UBOUND(ksupsmnf,3),'UBOUND WRONG IN CURL_HTOF')
!
!     (sqrt(g)*J^X)mn
!
      DO m = 0, mpol
         mp = sparity*m
         mo = m+moff
         DO n = -ntor, ntor
           np = sparity*n*nfp
           no = n+noff
           ksupsmnf(mo,no,nmin:nmax) = -mp*bsubvmnf(m,n,nmin:nmax)      &
                                     +  np*bsubumnf(m,n,nmin:nmax)
           ksupumnf(mo,no,nmin:nmax) =  np*bsubsmnf(m,n,nmin:nmax)      &
                                     -  db_vds(m,n,nmin:nmax)                     
           ksupvmnf(mo,no,nmin:nmax) = -mp*bsubsmnf(m,n,nmin:nmax)      &
                                     +  db_uds(m,n,nmin:nmax)
         ENDDO
      ENDDO
      
!Avoid fp failure in debug mode      
      ksupsmnf(:,:,nmax+1:UBOUND(ksupsmnf,3)) = 0
      ksupumnf(:,:,nmax+1:UBOUND(ksupsmnf,3)) = 0
      ksupvmnf(:,:,nmax+1:UBOUND(ksupsmnf,3)) = 0

      ksupsmnf(:,:,LBOUND(ksupsmnf,3):nmin-1) = 0
      ksupumnf(:,:,LBOUND(ksupsmnf,3):nmin-1) = 0
      ksupvmnf(:,:,LBOUND(ksupsmnf,3):nmin-1) = 0

      ksupsmnf(m0+moff,-ntor+noff:-1+noff,nmin:nmax) = 0
      ksupumnf(m0+moff,-ntor+noff:-1+noff,nmin:nmax) = 0
      ksupvmnf(m0+moff,-ntor+noff:-1+noff,nmin:nmax) = 0

      DEALLOCATE(bsubsmnf, bsubumnf,  bsubvmnf, db_uds, db_vds)   


      END SUBROUTINE CURL_HtoF


      SUBROUTINE CURL_HtoF_TARG(Bsubsmnh, Bsubumnh, Bsubvmnh,           &
                 Ksupsmnf, Ksupumnf, Ksupvmnf, ipar, nsmin0, nsmax0)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ipar, nsmin0, nsmax0
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
                Bsubsmnh, Bsubumnh, Bsubvmnh
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::                &
                Ksupsmnf, Ksupumnf, Ksupvmnf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), POINTER, DIMENSION(:,:,:) :: ksp, kup, kvp
!----------------------------------------------- 
      ksp=>ksupsmnf; kup=>ksupumnf; kvp=>ksupvmnf
      CALL CURL_HtoF(Bsubsmnh, Bsubumnh, Bsubvmnh,                      &
                     ksp, kup, kvp, ipar, nsmin0, nsmax0)

      END SUBROUTINE CURL_HtoF_TARG

          
      END MODULE utilities
          
