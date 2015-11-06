      SUBROUTINE allocate_nunv
      USE vmec_main
      USE vmec_params, ONLY: ntmax

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1
C-----------------------------------------------

      CALL free_mem_nunv

      ALLOCATE (bsubu0(nznt), rbsq(nznt), dbsq(nznt), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #1 in allocate_nunv'

#ifdef _ANIMEC
      ALLOCATE (pperp_ns(nznt), stat=istat1)
#endif
      ALLOCATE (rmn_bdy(0:ntor,0:mpol1,ntmax),
     1          zmn_bdy(0:ntor,0:mpol1,ntmax), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #2 in allocate_nunv'

!     PERSISTENT ARRAYS (DURATION OF PROGRAM)
      IF (lfreeb) CALL allocate_freeb_nunv

      END SUBROUTINE allocate_nunv

!  Separate out the freeb allocations so these can be allocated with out
!  needing to delete the other arrays.
      SUBROUTINE allocate_freeb_nunv
      USE vmec_main
      USE vmec_params, ONLY: ntmax
#ifdef _VACUUM2
      USE vac2_vacmod
#else
      USE vacmod
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1
C-----------------------------------------------

#ifdef _VACUUM2
      ALLOCATE (bsqsav(nznt,3), potvac(2*mnpd), stat=istat1)
#else
      ALLOCATE (amatsav(mnpd2*mnpd2),bvecsav(mnpd2),
     2          bsqsav(nznt,3), potvac(2*mnpd), raxis_nestor(nv),
     3          zaxis_nestor(nv), stat=istat1)

      ALLOCATE (bsubu_sur(nuv2), bsubv_sur(nuv2))
      ALLOCATE (bsupu_sur(nuv2), bsupv_sur(nuv2))
#endif
      IF (istat1.ne.0) STOP 'allocation error #1 in allocate_freeb_nunv'

      END SUBROUTINE allocate_freeb_nunv
