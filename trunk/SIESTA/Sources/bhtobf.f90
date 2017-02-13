!#if defined(SKS)
      SUBROUTINE bhalftobfull_par (bsupsijh, bsupuijh, bsupvijh,     &
                                   bsupsijf, bsupuijf, bsupvijf,     &
                                   nsmin, nsmax)
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
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(in) :: nsmin, nsmax
       REAL(dp), DIMENSION(ntheta,nzeta,nsmin:nsmax), INTENT(IN)    ::  &
           bsupsijh, bsupvijh
       REAL(dp), DIMENSION(ntheta,nzeta,nsmin:nsmax), INTENT(INOUT) ::  &
           bsupuijh 
       REAL(dp), DIMENSION(ntheta,nzeta,nsmin:nsmax), INTENT(OUT)   ::  &
           bsupsijf, bsupuijf, bsupvijf
!-----------------------------------------------
       REAL(dp) :: skston, skstoff
       INTEGER  :: ns_span, nblock
!-----------------------------------------------
       CALL second0(skston)
       ns_span=nsmax-nsmin+1
       nblock = SIZE(bsupsijh,1)*SIZE(bsupsijh,2)

!Interpolate b^X finite at the origin
!On exit, bsupXF = [nsmin, nsmax-1]
!Radial (s)
       CALL to_full_mesh_par(bsupsijh, bsupsijf, nblock, nsmin, nsmax)
       IF (nsmax .EQ. ns) bsupsijf(:,:,ns) = 0

!Poloidal(u): on input, bsupuijh = bsupuijh*jacobh
       CALL to_full_mesh_par(bsupuijh, bsupuijf, nblock, nsmin, nsmax)

!Toroidal (v)
       CALL to_full_mesh_par(bsupvijh, bsupvijf, nblock, nsmin, nsmax)

       bsupuijf = bsupuijf/jacobf(:,:,nsmin:nsmax)
       bsupuijh = bsupuijh/jacobh(:,:,nsmin:nsmax)

       CALL second0(skstoff)
       bhtobf_time=bhtobf_time+(skstoff-skston)

      END SUBROUTINE bhalftobfull_par
!#endif

      SUBROUTINE bhalftobfull (bsupsijh, bsupuijh, bsupvijh,         &
                               bsupsijf, bsupuijf, bsupvijf)
!     
!     WRITTEN 01-04-10 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Store BsupX on full mesh
!     NOTE:    On Entry, bsups,u_h at js=1 contain the "evolved" values of 
!              bsups,u_f at s=0 (js=1), as computed in UPDATE_BFIELD
!
       USE stel_kinds
       USE quantities, ONLY: ns, ntheta, nzeta, jacobf, jacobh
       USE descriptor_mod,ONLY: DIAGONALDONE, INHESSIAN 
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(in)    ::           &
           bsupsijh, bsupvijh
       REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(inout) ::           &
           bsupuijh
       REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(out)   ::           &
           bsupsijf, bsupuijf, bsupvijf
!-----------------------------------------------
       CALL bhalftobfull_par (bsupsijh, bsupuijh, bsupvijh,     &
                                   bsupsijf, bsupuijf, bsupvijf,     &
                                   1, ns)

       END SUBROUTINE bhalftobfull
