      SUBROUTINE allocate_radial
      USE parambs
      USE vmec0
      USE read_boozer_mod

      IMPLICIT NONE
      INTEGER :: istat

      ALLOCATE (flux(irdim), qsafety(irdim), aiogar(irdim),
     1    idx(irup), aipsi(irdim), gpsi(irdim), pres1(irdim),
     2    betar(irdim), dense(irdim), densi(irdim),
     3    tempe1(irdim), tempi1(irdim), lsurf(irdim),
     4    jlist(irdim), jlist_idx(irdim), stat=istat)
      IF (istat .ne. 0) RETURN

      ALLOCATE (amnfit(irup,-mboz_b:mboz_b,0:nboz_b), stat=istat)
      IF (istat .ne. 0) RETURN

      ALLOCATE (dibs(irup), aibs(irup), dibst(irup), aibst(irup),
     1    bsdense(irup), bsdensi(irup), bstempe(irup), bstempi(irup),
     2    bsdenste(irup), bsdensti(irup), bstempte(irup),
     3    bstempti(irup), capr(irup),
     4    caps(irup), h2(irup), ftrapped(irup), fpassing(irup),
     5    epsttok(irup), fttok(irup), gbsnorm(irup), aiterm1(irup),
     6    other1(irup),
     7    rhoar(irup), bsnorm(irup), fptok(irup), amain(irup),
     8    bmax1(irup), thetamax(irup), zetahmax(irup),
     9    ajBbs(irup), phip(irup), d_rho(irdim),
     A    b2avg(irup), stat = istat)

      IF (istat .ne. 0) STOP 'allocation error in allocate_radial'

      END SUBROUTINE allocate_radial
