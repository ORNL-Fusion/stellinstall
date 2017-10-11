      MODULE bhtobf
      
      CONTAINS
      
      SUBROUTINE bhalftobfull (bsupsijh, bsupuijh, bsupvijh,            &
                               bsupsijf, bsupuijf, bsupvijf,            &
                               pijh, pijf)
!     
!     WRITTEN 01-04-10 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Store BsupX on full mesh
!     NOTE:    On Entry, bsups,u_h at js=1 contain the "evolved" values of 
!              bsups,u_f at s=0 (js=1), as computed in UPDATE_BFIELD
!
       USE stel_kinds
       USE quantities, ONLY: ns, ntheta, nzeta, jacobf, jacobh

       USE descriptor_mod,ONLY: iam, nprocs
       USE nscalingtools, ONLY: startglobrow, endglobrow
       USE timer_mod, ONLY: bhtobf_time
       USE utilities, ONLY: to_full_mesh
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!       INTEGER, INTENT(in) :: nsmin, nsmax
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(IN)    ::        &
           bsupsijh, bsupvijh
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) ::        &
           bsupuijh, bsupsijf, bsupuijf, bsupvijf
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), OPTIONAL ::             &
           pijh, pijf
!-----------------------------------------------
       REAL(dp) :: skston, skstoff
       INTEGER  :: ns_span, nsmin, nsmax
!-----------------------------------------------
       CALL second0(skston)
       
       nsmin = LBOUND(bsupsijh,3)
       nsmax = UBOUND(bsupsijh,3)
       
       ns_span=nsmax-nsmin+1

!Interpolate b^X finite at the origin
!On exit, bsupXF = [nsmin, nsmax-1]
!Radial (s)
       CALL to_full_mesh(bsupsijh, bsupsijf)
       IF (nsmax .EQ. ns) bsupsijf(:,:,ns) = 0

!Poloidal(u): on input, bsupuijh = bsupuijh*jacobh
       CALL to_full_mesh(bsupuijh, bsupuijf)
       bsupuijf = bsupuijf/jacobf(:,:,nsmin:nsmax)
       bsupuijh = bsupuijh/jacobh(:,:,nsmin:nsmax)

!Toroidal (v)
       CALL to_full_mesh(bsupvijh, bsupvijf)
       
!pressure
       IF (PRESENT(pijh)) THEN
          IF (nsmin .NE. LBOUND(pijh,3)) CALL ErrorAbort('LBOUND WRONG!')
          IF (nsmax .NE. UBOUND(pijh,3)) CALL ErrorAbort('UBOUND WRONG!')
          CALL to_full_mesh(pijh, pijf)
       END IF

       CALL second0(skstoff)
       bhtobf_time=bhtobf_time+(skstoff-skston)

       END SUBROUTINE bhalftobfull

       END MODULE bhtobf
