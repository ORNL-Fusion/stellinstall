!  Small equilibria resolutions can be slower when run in parallel. Allow V3FIT
!  to be compiled with OpenMP support while disabling VMEC OpenMP support.
#ifdef NO_VMEC_OMP
#define OMP NO_OMP
#endif

      SUBROUTINE belicu(bx, by, bz, cos1, sin1, rp, zp)
      USE vacmod
      USE biotsavart
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nuv2), INTENT(in)  :: cos1, sin1, rp, zp
      REAL(rprec), DIMENSION(nuv2), INTENT(out) :: bx, by, bz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j
      REAL(rprec), DIMENSION(3) :: xpt, bvec
!-----------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(j,xpt,bvec) 
      DO j = 1,nuv2
         xpt(1) = rp(j) * cos1(j)
         xpt(2) = rp(j) * sin1(j)
         xpt(3) = zp(j) 

         CALL bsc_b (single_coil, xpt, bvec)

         bx(j) = bvec(1);  by(j) = bvec(2);  bz(j) = bvec(3)
      END DO
!$OMP END PARALLEL DO
      CALL cleanup_biotsavart

      END SUBROUTINE belicu
