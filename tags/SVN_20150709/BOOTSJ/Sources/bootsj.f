      SUBROUTINE bootsj(aibstot, extension, iunit_in)
c
!  The code BOOTSJ calculates the bootstrap current for 3D configurations.
!  The main part of the code, calculation of the geometrical factor GBSNORM,
!  was written by Johnny S. Tolliver of ORNL on the basis of
!  Ker-Chung Shaing's formulas. Other parts of the code, as well as modifications
!  to the Tolliver's part, were written by Paul Moroz of UW-Madison.
c
!  References:
c
!  1. K.C. Shaing, B.A. Carreras, N. Dominguez, V.E. Lynch, J.S. Tolliver
!    "Bootstrap current control in stellarators", Phys. Fluids B1, 1663 (1989).
!  2. K.C. Shaing, E.C. Crume, Jr., J.S. Tolliver, S.P. Hirshman, W.I. van Rij
!     "Bootstrap current and parallel viscosity in the low collisionality
!     regime in toroidal plasmas", Phys. Fluids B1, 148 (1989).
!  3. K.C. Shaing, S.P. Hirshman, J.S. Tolliver "Parallel viscosity-driven
!     neoclassical fluxes in the banana regime in nonsymmetri! toroidal
!     plasmas", Phys. Fluids 29, 2548 (1986).
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE parambs
      USE safe_open_mod
      USE trig
      USE vmec0
      USE fftpack, ONLY : fftfax, cftfax
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER     :: iunit_in
      REAL(rprec) :: aibstot
      CHARACTER*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nfax = 13
      INTEGER, PARAMETER :: indata0 = 7
      INTEGER, PARAMETER :: jbs_file=59, ans_file=18, ans_dat_file=19
      REAL(rprec) :: one = 1, p5 = .5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(nfax) :: ifaxu, ifaxv
      INTEGER :: ntrigu, ntrigv
      INTEGER :: irho, irho1, ierr, iunit, ijbs, ians, ians_plot
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: cputimes
      REAL(rprec) :: time1, timecpu, time2, r, x, al31t, gradbs1, 
     1               gradbs2, gradbs3, gradbs4,  al31s
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: al31
C-----------------------------------------------
c
c  define the starting CPU time
c
      CALL second0 (time1)
      timecpu = time1
      IF (lscreen) WRITE (*, 4) version_
    4 FORMAT(/,' Start BOOTSJ: Version ', a)

c  OPEN files

      iunit = indata0
      CALL safe_open(iunit, ierr, 'input.' // TRIM(extension), 'old',
     1     'formatted')
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening input file: input.', TRIM(extension)
         RETURN
      END IF

      ijbs = jbs_file
      CALL safe_open(ijbs, ierr, 'jBbs.'//TRIM(extension), 'replace',
     1     'formatted')

      ians = ans_file
      CALL safe_open(ians, ierr, 'answers.'//TRIM(extension), 'replace',
     1     'formatted')

      ians_plot = ans_dat_file
      CALL safe_open(ians_plot, ierr, 'answers_plot.'//TRIM(extension),
     1     'replace', 'formatted')

c  READ and initialize data

      CALL datain(TRIM(extension), iunit_in, iunit, ians)
      CLOSE (iunit)

      ntrigu = 3*nthetah/2 + 1
      ntrigv = 2*nzetah
      ALLOCATE (cputimes(irup))
      ALLOCATE (
     1  dmn(-mbuse:mbuse,0:nbuse), fmn(-mbuse:mbuse,0:nbuse),
     2  rfmn(-mbuse:mbuse,0:nbuse),alpha1mn(-mbuse:mbuse,0:nbuse),
     3  trigsv(ntrigv),trigsu(ntrigu), stat=irho)

      IF (irho .ne. 0) STOP 'allocation error in bootsj main'
c     convert jlist to idx form, THEN "and" the two lists
c     SPH: Corrected error here, writing off END of jlist error when jlist <= 1
c

      jlist_idx = 0
      DO irho = 1, irup
        IF (jlist(irho) .gt. 1) jlist_idx(jlist(irho)-1) = 1
        idx(irho) = idx(irho)*jlist_idx(irho)
      ENDDO


      l_boot_all = .true.

!  IF ANY of the boozer evaluation surfaces are missing, or not requested
!  THEN set l_boot_all=false, total current can not be calculated

      DO irho=1, irup
         IF(idx(irho) .eq. 0) l_boot_all = .false.
      ENDDO
      IF(.not.l_boot_all) THEN
         IF (lscreen) WRITE (*,*) 'partial surface evaluation'
         WRITE (ians,*) 'partial surface evaluation'
      ENDIF


      CALL fftfax (nthetah, ifaxu, trigsu)
      CALL cftfax (nzetah, ifaxv, trigsv)


c start main radial loop

      DO irho = 1, irup

c  initialize timing for radius step

         CALL second0 (time2)
         timecpu = time2 - time1
         cputimes(irho) = timecpu

c  IF there is no boozer information available, skip radial point

         IF(idx(irho) .eq. 0)  CYCLE
         irho1 = irho - 1
         r = SQRT(rhoar(irho) + 1.E-36_dp)

c  initialize  angle grids grid for first radial evaluation.  For this
c  and ALL subsequent radial evaluate B and related and quantites as well
c  as plasma derivatives and the tokamak trapped fraction.

         CALL bongrid(irho, ians)

c  calculate bootstrap current for equivalent tokamak

         x = fttok(irho)/(fptok(irho)+1.E-36_dp)

         al31t = al31(x,zeff1,alphae,alphai)

c  calculate gradient factors overALL normalization inlcuding q and
c  boozer g.

         CALL grad (gradbs1, gradbs2, gradbs3, gradbs4, irho)

         bsdenste(irho) = gradbs1*al31t          !due to dens gradient
         bsdensti(irho) = gradbs2*al31t          !due to dens gradient
         bstempte(irho) = gradbs3*al31t          !due to temp gradient
         bstempti(irho) = gradbs4*al31t          !due to temp gradient

         dibst(irho) = bsdenste(irho) + bsdensti(irho) + bstempte(irho)
     1       + bstempti(irho)                    !total Jbst

         IF (l_boot_all) THEN
            IF (irho .eq. 1) THEN
               aibst(1) = dibst(1)*d_rho(1)
            ELSE
               aibst(irho) = aibst(irho1)+dibst(irho)*d_rho(irho)
            ENDIF
         END IF

c  Now start the general evaluation.
c  Find coefficients d(n,m) and evaluate the fraction trapped.

         CALL denmf (trigsu, trigsv, ifaxu, ifaxv, irho)

c  Evaluate R, S, and H2.

         CALL caprsh2(irho)

c  evaluate the integral term.

         CALL woflam (trigsu, trigsv, ifaxu, ifaxv, irho)

c  evaluate the summation term in w(lambda) that does not depend on lambda.
c  note that while the paper shows an itegral, it is the fpassing integral
c  that cancels the fpassing in the over all multiplier

         CALL othersums(irho)

c  Calculate the final answer

         amain(irho) = p5*(one - aiogar(irho)/qsafety(irho)) + p5*(one +
     1      aiogar(irho)/qsafety(irho))*h2(irho)

         gbsnorm(irho) = amain(irho) + other1(irho) +  aiterm1(irho)

c- derivative of the enclosed bootstrap current over the normalized
c  toroidal flux, dIbs/ds, and the enclosed Ibs (in MA)

         x = ftrapped(irho)/(fpassing(irho)+1.E-36_dp)
         al31s = al31(x,zeff1,alphae,alphai)
         CALL grad (gradbs1, gradbs2, gradbs3, gradbs4, irho)
c
c                                                !due to dens gradient
         bsdense(irho) = gbsnorm(irho)*gradbs1*al31s
c                                                !due to dens gradient
         bsdensi(irho) = gbsnorm(irho)*gradbs2*al31s
c                                                !due to temp gradient
         bstempe(irho) = gbsnorm(irho)*gradbs3*al31s
c                                                !due to temp gradient
         bstempi(irho) = gbsnorm(irho)*gradbs4*al31s

         dibs(irho) = bsdense(irho) + bsdensi(irho) + bstempe(irho) +
     1      bstempi(irho)                        !total Jbst (dI/ds)

c   convert to j dot b.  2*mu0 from beta, 10**6 from ma, psimax is 1/dpsi/ds
c   and sign_jacobian takes out the sign previously used for dpsi to da
ccccc  flux was changed to real flux in boot vmec so that the psimax now
ccccc  needs an additional sign_jacobian so that the two will cancel
             ajBbs(irho) = (2.0e6_dp)*mu0*dibs(irho)*
     1      (pres1(irho)/betar(irho))/psimax

         IF (l_boot_all) THEN
            IF (irho .eq. 1) THEN
               aibs(1) = dibs(1)*d_rho(1)
            ELSE
               aibs(irho) = aibs(irho1) + dibs(irho)*d_rho(irho)
            ENDIF
         END IF

c- the ratio of bootstrap current to that in ESD (equivalent symmetric device):

         bsnorm(irho) = dibs(irho)/(dibst(irho)+1.E-36_dp)

c  get time at END of loop IF completed

         CALL second0 (time2)
         timecpu = time2 - time1
         cputimes(irho) = timecpu
c

      END DO
c
c- Output answers for BOOTSJ
c
      CALL output (cputimes, aibstot, ijbs, ians, ians_plot)
      CLOSE(ians)
      CLOSE(ians_plot)
      CLOSE(ijbs)

      CALL deallocate_all
      IF (lscreen) WRITE (*,400) (cputimes(irup)-cputimes(1))
  400 FORMAT(1x,'Finished BOOTSJ, time =  ', f8.3, '  sec')

      DEALLOCATE (cputimes, trigsu, trigsv)
      DEALLOCATE (dmn, fmn, rfmn, alpha1mn)

      END SUBROUTINE bootsj
