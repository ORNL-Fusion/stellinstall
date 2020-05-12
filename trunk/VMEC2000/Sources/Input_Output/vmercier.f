      MODULE vmercier
      USE vparams, ONLY: nsd, rprec, dp
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nsd) ::
     1   Dshear, Dwell, Dcurr, Dmerc, Dgeod,
     2   E1_GL, E2pH2_GL, F1_GL, F2_GL, F3_GL, H1_GL, Dmerc_GL
      END MODULE vmercier
