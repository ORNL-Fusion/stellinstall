      MODULE utilities

      INTERFACE to_full_mesh
         MODULE PROCEDURE to_full_mesh2, to_full_mesh3
      END INTERFACE

      CONTAINS
      
      SUBROUTINE GradientHalf(gradienth, vectorf)
      USE stel_kinds
      USE island_params, ONLY: ohs=>ohs_i
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), ALLOCATABLE :: vectorf(:,:,:)
      REAL(dp), ALLOCATABLE :: gradienth(:,:,:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: nsmin, nsmax
!-----------------------------------------------
!     Contraction map
!     INPUT (vectorf)   : [nsmin:nsmax]
!     OUTPUT(gradienth) : [nsmin+1:nsmax]
      nsmin = LBOUND(vectorf,3); nsmax = UBOUND(vectorf,3)

      gradienth(:,:,nsmin+1:nsmax) = ohs*(vectorf(:,:,nsmin+1:nsmax)    & 
                                   -      vectorf(:,:,nsmin:nsmax-1))

      END SUBROUTINE GradientHalf
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE GradientFull(gradientf, vectorh)
      USE stel_kinds
      USE island_params, ONLY: ohs=>ohs_i, ns=>ns_i
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), ALLOCATABLE  :: vectorh(:,:,:)
      REAL(dp), ALLOCATABLE :: gradientf(:,:,:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: nsmin, nsmax
!-----------------------------------------------
!     Contraction map
!     INPUT (vectorh)   : [nsmin:nsmax]
!     OUTPUT(gradientf) : [nsmin:nsmax-1]
      nsmin = LBOUND(vectorh,3); nsmax = UBOUND(vectorh,3)

      gradientf(:,:,nsmin:nsmax-1) = ohs*(vectorh(:,:,nsmin+1:nsmax)    &
                                   -      vectorh(:,:,nsmin:nsmax-1))

      IF (nsmin .EQ. 1) gradientf(:,:,nsmin) = 0

      IF (nsmax.EQ.ns .AND. nsmax.NE.nsmin) THEN
         gradientf(:,:,nsmax) = gradientf(:,:,nsmax-1)
      END IF

      END SUBROUTINE GradientFull
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE to_half_mesh(vecf, vech)
      USE stel_kinds
      USE island_params, ns=>ns_i
      IMPLICIT NONE        
!  
!     This subroutine migrates a full mesh quantity to the half mesh
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), ALLOCATABLE  :: vecf(:,:,:)   
      REAL(dp), ALLOCATABLE :: vech(:,:,:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: nsmin, nsmax
      REAL(dp), PARAMETER   :: p5 = 0.5_dp
!-----------------------------------------------
!     Contraction map
!     INPUT (vecf)   : [nsmin:nsmax]
!     OUTPUT(vech)   : [nsmin+1:nsmax]

      nsmin = LBOUND(vech,3); nsmax = UBOUND(vech,3)

      vech(:,:,nsmin+1:nsmax) = p5*(vecf(:,:,nsmin+1:nsmax)             &
                              +     vecf(:,:,nsmin:nsmax-1))
      IF (nsmin .EQ. 1) vech(:,:,nsmin) = 0

      END SUBROUTINE to_half_mesh 
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh3(vech, vecf)
      USE stel_kinds
      USE island_params, ns=>ns_i
      IMPLICIT NONE        
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
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
      INTEGER     :: nsmin, nsmax
!-----------------------------------------------
!     Contraction map
!     INPUT (vech)   : [nsmin:nsmax]
!     OUTPUT(vecf)   : [nsmin:nsmax-1]

      nsmin = LBOUND(vech,3)
      nsmax = UBOUND(vech,3)

      vecf(:,:,nsmin:nsmax-1) = p5*(vech(:,:,nsmin:nsmax-1)             &
                              +     vech(:,:,nsmin+1:nsmax))
      IF (nsmin .EQ. 1)  vecf(:,:,1) = vech(:,:,2)
      IF (nsmax .EQ. ns) vecf(:,:,ns)= vech(:,:,ns)

      END SUBROUTINE to_full_mesh3
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh2(vech, vecf)
      USE stel_kinds
      USE island_params, ns=>ns_i
      IMPLICIT NONE        
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
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
      INTEGER     :: nsmin, nsmax
!-----------------------------------------------
!     Contraction map
!     INPUT (vech)   : [nsmin:nsmax]
!     OUTPUT(vecf)   : [nsmin:nsmax-1]

      nsmin = LBOUND(vech,2)
      nsmax = UBOUND(vech,2)

      vecf(:,nsmin:nsmax-1) = p5*(vech(:,nsmin:nsmax-1)                 &
                            +     vech(:,nsmin+1:nsmax))
      IF (nsmin .EQ. 1)  vecf(:,1) = vech(:,2)
      IF (nsmax .EQ. ns) vecf(:,ns)= vech(:,ns)

      END SUBROUTINE to_full_mesh2

          
      END MODULE utilities
