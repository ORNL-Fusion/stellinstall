      SUBROUTINE deallocate_all
      USE parambs
      USE vmec0
      IMPLICIT NONE

      IF (ALLOCATED(amnfit)) DEALLOCATE (amnfit)
      IF (ALLOCATED(dibs)) DEALLOCATE (dibs, aibs, dibst, aibst,
     1    bsdense, bsdensi, bstempe, bstempi, bsdenste, bsdensti,
     2    bstempte, bstempti, qsafety, capr, caps, h2,
     3    ftrapped, fpassing, epsttok, fttok, gbsnorm, aiterm1,
     4    other1, rhoar, bsnorm, fptok, amain,
     5    bmax1, thetamax, zetahmax, d_rho, b2avg,
     6    ajBbs, phip ,theta, zetah)
      IF (ALLOCATED(flux)) DEALLOCATE (flux, aiogar, aipsi,
     1    gpsi, pres1, betar, dense, densi, tempe1, tempi1, lsurf,
     2    jlist, jlist_idx)

      DEALLOCATE(idx)
      DEALLOCATE(bfield, b2obm, gsqrt_b, sinmi, sinnj, cosmi, cosnj,
     1  bfieldm, sinmim, sinnjm, cosmim, cosnjm)

      END SUBROUTINE deallocate_all
