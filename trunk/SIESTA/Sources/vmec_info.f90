!>
!!  \brief Module for reading in VMEC input and computing metric elements
!!   on the SIESTA (sqrt-flux) radial grid.
!<
      MODULE vmec_info
      USE v3_utilities, ONLY: assert
      USE stel_kinds

      IMPLICIT NONE
!     
!     WRITTEN 06-27-06 BY S. P. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: COMPUTES AND STORES THE REAL-SPACE METRIC ELEMENTS, JACOBIAN BASED 
!     ON A SQRT(FLUX) SPLINED VMEC - COORDINATE SYSTEM
!

!     VARIABLE DECLARATIONS (3D for now; will convert to 2D or 1D as needed)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), DIMENSION(:,:), ALLOCATABLE ::                          &
                   rmnc_spline, zmns_spline, lmns_spline,               &
                   rmns_spline, zmnc_spline, lmnc_spline

      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                        &
                   rmnc_i, zmns_i, lmns_i,                              &
                   rmns_i, zmnc_i, lmnc_i

      END MODULE vmec_info
