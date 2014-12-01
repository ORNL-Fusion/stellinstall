      SUBROUTINE output(cputimes, aibstot, ijbs, ians, ians_plot)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE parambs
      USE vmec0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ijbs, ians, ians_plot
      REAL(rprec) :: aibstot
      REAL(rprec), DIMENSION(*) :: cputimes
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1,
     1  D18 = 1.E-18_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: aibsttot
C-----------------------------------------------
      WRITE (ijbs, 87)
      DO i = 1, irup
         IF(idx(i) .eq. 1) THEN
            WRITE (ijbs, *) (i+1), ajBbs(i), bsnorm(i)
         ENDIF
      ENDDO
   87 FORMAT(1x,'surface #  <Jbs_dot_B>    <Jbs>/<Jbs-tok>')

      WRITE (ians, 90)
      DO i=1,irup
         IF(idx(i) .eq. 1) THEN
            WRITE (ians, 100) rhoar(i),capr(i),caps(i),ftrapped(i),
     1         h2(i), amain(i),aiterm1(i),
     2         other1(i), gbsnorm(i)
         ENDIF
      ENDDO
   90 FORMAT(/,1x,
     1   '  s        R           S        ftrapped       H2     ',
     2   '    amain      lam int     other      gnorm')
  100 FORMAT(1x,0p,f5.3,1p,8e12.4)


      WRITE (ians, 101)
      DO i=1,irup
         IF(idx(i) .eq. 1) THEN
            WRITE (ians, 102) rhoar(i),qsafety(i),thetamax(i),
     1          zetahmax(i),bmax1(i),ftrapped(i)/fttok(i),
     2          fptok(i)/fpassing(i),cputimes(i)
         ENDIF
      ENDDO
  101 FORMAT(/,1x,'  s        Q        thetamax     zetamax   ',
     1   '    BMAX    ft/fttok    fptok/fp    cpu secs')
  102 FORMAT(1x,0p,f5.3,1p,8e12.4)

      WRITE(ians,*) '   s    gnorm       jbsnorm  ',
     1   '   dI/ds      I(s)(hm)  '
     2    ,'  j_grad_ne   j_grad_ni   j_grad_Te   j_grad_Ti '
      DO i=1,irup
         IF(idx(i) .eq. 1) THEN
         WRITE(ians,32)rhoar(i), gbsnorm(i),bsnorm(i), dibs(i),
     1   aibs(i),  bsdense(i), bsdensi(i), bstempe(i),
     1   bstempi(i)
         ENDIF
      ENDDO
   32 FORMAT(1x,0p,f5.3,1p,9e12.4)


      WRITE (ians, 104)
      DO i=1,irup
         IF(idx(i) .eq. 1) THEN
            WRITE (ians, 105)i,tempe1(i),tempi1(i),dense(i),
     1      dense(i)/zeff1,   betar(i),ajBbs(i)
         END IF
      ENDDO
  104 FORMAT(/,1x,'     Te          Ti          Ne         Ni
     1Beta          jB')
  105 FORMAT(1x,i3,1p,6e12.4)


      IF(l_boot_all .and. lscreen) THEN
         aibstot = aibs(irup)
         WRITE (*, '(a,f12.7,a,i3,a,i3)') ' Total bootstrap current =',
     1      aibstot , ' MA'
         aibsttot = aibst(irup)
         aibstot = aibstot*1.e6_dp                     !in Amperes
      ENDIF

  107 FORMAT('Ibs = ',f7.4,' MA')
  175 FORMAT(2x,1p,5e14.6)

      DO i = 1, irup
         IF(idx(i) .eq. 1) THEN
           WRITE(ians_plot, *) rhoar(i),  gbsnorm(i), amain(i),
     1      aiterm1(i), other1(i),
     2      dibs(i), bsdense(i), bsdensi(i), bstempe(i), bstempi(i),
     3      qsafety(i), ftrapped(i), bsnorm(i),
     4      tempe1(i),tempi1(i),dense(i), dense(i)/zeff1, betar(i),
     5      ajBbs(i)
         ENDIF
      ENDDO

      END SUBROUTINE output
