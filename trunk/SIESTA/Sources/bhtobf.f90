      MODULE bhtobf
      USE v3_utilities, ONLY: assert      
      CONTAINS
      
      SUBROUTINE bhalftobfull (bsupsijh, bsupuijh, bsupvijh,            &
                               bsupsijf, bsupuijf, bsupvijf,            &
                               pijh, pijf)
!     
!     WRITTEN 01-04-10 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Store BsupX (and p if persent) on full mesh
!     NOTE:    On Entry, bsups,u_h at js=1 contain the "evolved" values of 
!              bsups,u_f at s=0 (js=1), as computed in UPDATE_BFIELD
!
       USE stel_kinds
       USE quantities, ONLY: ns, jacobf, jacobh
       USE timer_mod
       USE utilities, ONLY: to_full_mesh
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT)    ::     &
           bsupsijh, bsupvijh
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) ::        &
           bsupuijh, bsupsijf, bsupuijf, bsupvijf
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), OPTIONAL ::             &
           pijh, pijf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(dp) :: skston, skstoff
       INTEGER  :: nsmin, nsmax
!-----------------------------------------------
#undef _OLD_AV
!#define _OLD_AV
       CALL second0(skston)
       
       nsmin = LBOUND(bsupsijh,3)
       nsmax = UBOUND(bsupsijh,3)
       
!Interpolate b^X finite at the origin
!On exit, bsupXF = [nsmin, nsmax-1]
!Radial (s)
#if defined(_OLD_AV)
       bsupsijh = bsupsijh/jacobh(:,:,nsmin:nsmax)
#endif       
       CALL to_full_mesh(bsupsijh, bsupsijf)
       IF (nsmin .EQ. 1) bsupsijf(:,:,1) = bsupsijh(:,:,1)
#if !defined(_OLD_AV)       
       bsupsijf = bsupsijf/jacobf(:,:,nsmin:nsmax)
       bsupsijh = bsupsijh/jacobh(:,:,nsmin:nsmax)
#endif       
!Poloidal(u): on input, bsupuijh = bsupuijh*jacobh
       CALL to_full_mesh(bsupuijh, bsupuijf)
       IF (nsmin .EQ. 1) bsupuijf(:,:,1) = bsupuijh(:,:,1)
       bsupuijf = bsupuijf/jacobf(:,:,nsmin:nsmax)
       bsupuijh = bsupuijh/jacobh(:,:,nsmin:nsmax)

!Toroidal (v)
#if defined(_OLD_AV)
       bsupvijh = bsupvijh/jacobh(:,:,nsmin:nsmax)
#endif
       CALL to_full_mesh(bsupvijh, bsupvijf)
       IF (nsmin .EQ. 1) bsupvijf(:,:,1) = bsupvijh(:,:,1)

#if !defined(_OLD_AV)       
       bsupvijf = bsupvijf/jacobf(:,:,nsmin:nsmax)
       bsupvijh = bsupvijh/jacobh(:,:,nsmin:nsmax)
#endif       
!pressure
       IF (PRESENT(pijh)) THEN
          CALL ASSERT(SIZE(bsupsijh,3).EQ.SIZE(pijh,3),                 &
                      'bhtobf pijh SIZE WRONG')
#if defined(_OLD_AV)
          pijh = pijh/jacobh(:,:,nsmin:nsmax)
#endif          
          CALL to_full_mesh(pijh, pijf)
          IF (nsmin .EQ. 1) pijf(:,:,1) = pijh(:,:,1)
#if !defined(_OLD_AV)       
          pijf = pijf/jacobf(:,:,nsmin:nsmax)
          pijh = pijh/jacobh(:,:,nsmin:nsmax)
#endif          
       END IF

       IF (nsmax .EQ. ns) bsupsijf(:,:,ns) = 0
       
       CALL second0(skstoff)
#if defined(SKS)
       bhtobf_time=bhtobf_time+(skstoff-skston)
#endif
       END SUBROUTINE bhalftobfull

       END MODULE bhtobf
