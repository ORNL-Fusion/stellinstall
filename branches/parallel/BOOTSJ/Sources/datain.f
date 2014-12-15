      SUBROUTINE datain(extension, iunit_in, iunit, ians)
c--
c  READ and initialize data for BOOTSJ
c--
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE parambs
      USE vmec0
      USE read_boozer_mod
      USE bootsj_input
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iunit, ians, iunit_in
      CHARACTER*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, n, m, i1, ntheta_min, nzeta_min,
     1   mn, ir, idum
      REAL(rprec) :: temp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: work
      REAL(rprec) :: status, tempe0, tempi0, pres10, pres0, a, a1
C-----------------------------------------------
!
!     Read in boozer |B|, iota, phip, pres AND ALLOCATE needed arrays
!
      CALL read_boozer(extension)

      jlist = 0
      READ (iunit_in, *,iostat=idum) jlist
      IF (idum .gt. 0) STOP ' Error reading bootsj input file in DATAIN'

      CLOSE (iunit_in)


!
!
!- give default values of control data otherwise read from the file 'bootin'
!
!  nrho is an obsolete variable that is no longer used.  it is set to ns-1
!  where ns is the number of whole mesh points used in VMEC
!
!  nrho is an obsolete variable that is no longer used.  it is set to ns-1
!  where ns is the number of whole mesh points used in VMEC
!
      nrho = 30                !number of rho values to use
      mbuse = 6                !number of m (poloidal) terms in B field.
      nbuse = 4                !number of nzetah (toroidal) terms in B field.
      zeff1 = 1.0_dp           !effective ion charge
      dens0 = 0.3_dp           !central electron density in 10**20 m-3
      teti = 2.0_dp            !ratio of Te/Ti for electron to ion
                               !temperature profiles
      tempres = -one           !tempe1(i)=pres(ir)**tempres
                               !IF(tempres.eq.-1.0_dp) THEN
                               !tempe1(i)=SQRT(pres(i))
      damp = -0.01_dp          !superceded by damp_bs
      damp_bs = -0.01_dp       !damping factor for resonances
      isymm0 = 0               !IF ne.0 THEN force a symmetric-device calculation,
      tempe1 = one             !1 keV DEFAULT temperature
      tempi1 = zero
      ate    = zero
      ati    = zero
      ate(0) = -one
      ati(0) = -one
!
!- Read control data
!
      CALL read_namelist (iunit, i, 'bootin')
      IF (i .ne. 0 ) STOP 'Error reading bootin NAMELIST'

!  If the new damping variable is not READ in, set to to the value this gives
!  an unchanged damping factor, namely damp_bs**2 = damp.  If, in addition,
!  damp is not READ in, THEN set damp_bs = 0.001

      IF(damp_bs .lt. zero) THEN !in this CASE no damp_bs was READ in
        IF(damp .gt. zero) THEN
           damp_bs = SQRT(damp)
        ELSE
           damp_bs = 0.001_dp      !in this CASE no damp was READ in
        ENDIF
      ENDIF
      teti = ABS(teti)

!
!     CHECK DIMENSION SIZE CONSISTENCY
!
      IF(nboz_b .lt. nbuse) THEN
         IF (lscreen) WRITE(*,*) 'nbuse > nbos_b, nbuse = nboz_b'
         nbuse = nboz_b
      ENDIF
      IF(mboz_b .lt. mbuse) THEN
         IF (lscreen) WRITE(*,*) 'mbuse > mbos_b, mbuse = mboz_b'
         mbuse = mboz_b
      ENDIF

      nzeta_min = 2*nbuse + 1
      ntheta_min = 2*mbuse + 1

      DO i = 0, 6
         nzetah = 4*2**i
         IF(nzetah .gt. nzeta_min) EXIT
         nzetah = 2*2**i * 3
         IF(nzetah .gt. nzeta_min) EXIT
      ENDDO

      DO i = 0, 6
         nthetah = 4*2**i
         IF(nthetah .gt. ntheta_min) EXIT
         nthetah = 2*2**i * 3
         IF(nthetah .gt. ntheta_min) EXIT
      ENDDO

      IF(lscreen) PRINT *, 'mbuse = ',mbuse,'nbuse = ',
     1    nbuse,'nthetah = ',nthetah, ' nzetah = ',nzetah

   90 FORMAT(5e16.8)

!     convert bmn's to amnfit's (positive m's ONLY to positive n's only)

      amnfit = zero

      lsurf = .false.
      status = TINY(a1)
      DO ir = 1, irup
         DO mn = 1,mnboz_b
            m = ixm_b(mn)
            n = ixn_b(mn)/nfp_b
            IF (m .gt. mboz_b) STOP 'boozmn indexing conflict, m'
            IF (ABS(n) .gt. nboz_b) STOP 'boozmn indexing conflict, n'
            IF (n.lt.0) THEN
               m = -m
               n = -n
            END IF
            IF (m.eq.0 .and. n.eq.0 .and. bmnc_b(mn,ir).gt.status)
     1         lsurf(ir) = .true.
            amnfit(ir,m,n) = bmnc_b(mn,ir+1)     !!2nd half grid == 1st grid pt. here
         END DO
      END DO

      CALL read_boozer_deallocate

      zeff1 = MAX(one,zeff1)
c  setup s grid.  A nonuniform mesh is allowed.

      psimax = MAXVAL(ABS(flux))

c  we need to keep the sign of psimax
      IF(flux(irup) .lt. zero) psimax = -psimax

c  first normalize the mesh to 1
c  and THEN shift from the full to half mesh with the first point on 1, not 2
c  also need to calulate the deltas for differencing on the vmec grid
c  special values at the ENDs will be calculated in bongrid

      DO ir = 1, irup
        rhoar(ir) = 0.5_dp*(flux(ir) + flux(ir+1))/psimax
        d_rho(ir) = (flux(ir+1) - flux(ir))/psimax
      END DO

c      WRITE (ijbs, 130) irup, psimax
c  130 FORMAT(' Last flux surface used is number ',i2,' with PSIMAX =',1p
c     1   e11.4,' WB')
c
c- Switch sign of q, IF desired.
c
      IF (iotasign .lt. 0) THEN
         qsafety(:irup) = iotasign*qsafety(:irup)
      ENDIF

c
c                                                !Boozer I/g values
      aiogar(:irup) = aipsi(:irup)/(gpsi(:irup)+1.0e-36_dp)
c
      CALL positiv (pres1, irup, 2) !to be sure that arrays are positive
      CALL positiv (betar, irup, 2)
c  evaluate electron and ion temperature profiles on vmec mesh
      IF (ANY(ate .ne. zero)) THEN
         DO ir = 1, irup
            tempe1(ir) = temp(rhoar(ir), ate)
         END DO
      END IF
      IF (ANY(ati .ne. zero)) THEN
         DO ir = 1, irup
            tempi1(ir) = temp(rhoar(ir), ati)
         END DO
      END IF

      tempe0 = tempe1(1)            !central electron temperature in keV
      tempi0 = tempi1(1)                 !central ion temperature in keV
      pres10 = pres1(1)                    !central total pressure in Pa
      pres0 = 1.6022E4_DP
c  If the leading coefficient on either te or ti is negative, the intent is
c  to assume a maximum density, and then calcuate the profiles using a
c  power law given by ABS(tempres).  this coefficient must originally be
c  be negative and ca not be greater than one.  That is to say, the profiles
c  are determined by ne(0) and tempres.  The ratio of Te to Ti is given by
c  teti.
      IF (tempe0.le.zero .or. tempi0.le.zero) tempres = ABS(tempres)
      tempres = MIN(one,tempres)
c
      IF (tempres .ge. zero) THEN   !in that CASE, calculate the temperature
         teti = teti + 1.E-36_dp
         a = one + one/(zeff1*teti)
c                                   !central electron temperature in keV
         tempe0 = pres10/(a*pres0*dens0)
         tempi0 = tempe0/teti            !central ion temperature in keV
         tempe1(:irup) = pres1(:irup)**tempres
c                             !suggested Te and Ti profiles are the same
         tempi1(:irup) = tempe1(:irup)
         a = tempe0/tempe1(1)
         a1 = tempi0/tempi1(1)
         tempe1(:irup) = tempe1(:irup)*a
         tempi1(:irup) = tempi1(:irup)*a1
      ENDIF
c
      CALL positiv (tempe1, irup, 2)
      CALL positiv (tempi1, irup, 2)

      dense(:irup) = pres1(:irup)/(pres0*(tempe1(:irup)+tempi1(:irup)/
     1   zeff1)+1.E-36_dp)
      ALLOCATE(work(irdim))
      CALL smooth1 (dense, 1, irup, work, zero)
      CALL positiv (dense, irup, 2)
      i1 = irup - 1
      a = tempe1(irup) + tempi1(irup)/zeff1
      a1 = tempe1(i1) + tempi1(i1)/zeff1
      dense(irup) = dense(i1)*a1*betar(irup)/(a*betar(i1)+1.E-36_dp)
c  the above game is apparently to control the spline derivatives at the
c  edge.
      densi(:irup) = dense(:irup)/zeff1
      dens0 = dense(1)                         !central electron density


c
c- Echo input data to output file.
c
      WRITE (ians, bootin)

      DEALLOCATE (work)

      END SUBROUTINE datain
