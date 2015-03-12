!  Small equilibria resolutions can be slower when run in parallel. Allow V3FIT
!  to be compiled with OpenMP support while disabling VMEC OpenMP support.
#ifdef NO_VMEC_OMP
#define OMP NO_OMP
#endif

      SUBROUTINE alias(gcons, ztemp, gcs, gsc, gcc, gss)
      USE vmec_main
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(out) ::
     1   gcons
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(in)  :: ztemp
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1) :: gcs, gsc, gcc, gss
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, i, ir, jk, jka, n, k, js, l, nsz
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work, gcona
C-----------------------------------------------
      nsz=ns*nzeta
      ALLOCATE (work(nsz,4), gcona(nsz,ntheta3))

!      gcons = 0
!      gcona = 0

!      gcs = 0;  gsc = 0
!      gcc = 0;  gss = 0

      DO m = 1, mpol1-1
!         work = 0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(i,ir,jk,jka)
         DO jk = 1,nsz
            work(jk,:) = 0
            DO i = 1, ntheta2
               work(jk,1) = work(jk,1) + ztemp(jk,i)*cosmui(i,m)
               work(jk,2) = work(jk,2) + ztemp(jk,i)*sinmui(i,m)
               IF (.NOT.lasym) CYCLE
               ir = ntheta1 + 2 - i
               IF (i .EQ. 1) ir = 1
               jka = ireflect(jk)
               work(jk,3) = work(jk,3) + ztemp(jka,ir)*cosmui(i,m)
               work(jk,4) = work(jk,4) + ztemp(jka,ir)*sinmui(i,m)
            END DO
         END DO
!$OMP END DO

!$OMP DO PRIVATE(n)
         DO n=0,ntor
            gcs(:,n,m) = 0; gsc(:,n,m) = 0
            IF (.NOT.lasym) CYCLE
            gss(:,n,m) = 0; gcc(:,n,m) = 0
         END DO
!$OMP END DO

!$OMP DO PRIVATE(k,l,n,js)
         DO n = 0, ntor
            DO k = 1, nzeta
               l = ns*(k-1)
               IF (.NOT.lasym) THEN
               DO js = 2,ns
                  gcs(js,n,m) = gcs(js,n,m) + tcon(js)*work(js+l,1)*
     1               sinnv(k,n)
                  gsc(js,n,m) = gsc(js,n,m) + tcon(js)*work(js+l,2)*
     1               cosnv(k,n)
               END DO
               ELSE
               DO js = 2,ns
                  gcs(js,n,m) = gcs(js,n,m) + p5*tcon(js)*sinnv(k,n)*
     1               (work(js+l,1)-work(js+l,3))
                  gsc(js,n,m) = gsc(js,n,m) + p5*tcon(js)*cosnv(k,n)*
     1               (work(js+l,2)-work(js+l,4))
                  gss(js,n,m) = gss(js,n,m) + p5*tcon(js)*sinnv(k,n)*
     1               (work(js+l,2)+work(js+l,4))
                  gcc(js,n,m) = gcc(js,n,m) + p5*tcon(js)*cosnv(k,n)*
     1               (work(js+l,1)+work(js+l,3))
               END DO
               END IF
            END DO
         END DO
!$OMP END DO
!
!        INVERSE FOURIER TRANSFORM DE-ALIASED GCON
!
!         work = 0
!$OMP DO PRIVATE(jk)
         DO jk = 1,nsz
            work(jk,:) = 0
         END DO
!$OMP END DO

!$OMP DO PRIVATE(n,k,l,js)
         DO k = 1, nzeta
            l = ns*(k-1)
            DO n = 0, ntor
               DO js = 2, ns
                  work(js+l,3) = work(js+l,3) + gcs(js,n,m)*sinnv(k,n)
                  work(js+l,4) = work(js+l,4) + gsc(js,n,m)*cosnv(k,n)
               END DO
               IF (.NOT.lasym) CYCLE
               DO js = 2, ns
                  work(js+l,1) = work(js+l,1) + gcc(js,n,m)*cosnv(k,n)
                  work(js+l,2) = work(js+l,2) + gss(js,n,m)*sinnv(k,n)
               END DO
            END DO
         END DO
!$OMP END DO

!$OMP DO PRIVATE(i,jk)
         DO i = 1, ntheta2
            IF (m .EQ. 1) gcons(:,i) = 0
            DO jk = 1, nsz
               gcons(jk,i) = gcons(jk,i) + (work(jk,3)*cosmu(i,m)
     1                     + work(jk,4)*sinmu(i,m))*faccon(m)
            END DO
            IF (.NOT.lasym) CYCLE
            IF (m .EQ. 1) gcona(:,i) = 0
            DO jk = 1, nsz
               gcona(jk,i) = gcona(jk,i) + (work(jk,1)*cosmu(i,m)
     1                     + work(jk,2)*sinmu(i,m))*faccon(m)
            END DO
         END DO
!$OMP END DO
!$OMP END PARALLEL
      END DO

      IF (lasym) THEN
!
!     EXTEND GCON INTO THETA = PI,2*PI DOMAIN
!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(i, ir, jk, jka)
         DO i = 1 + ntheta2, ntheta1
            ir = ntheta1 + 2 - i
            DO jk = 1, nsz
               jka = ireflect(jk)
               gcons(jk,i) = -gcons(jka,ir) + gcona(jka,ir)
            END DO
         END DO
!$OMP END DO
!
!     ADD SYMMETRIC, ANTI-SYMMETRIC PIECES IN THETA = 0,PI DOMAIN
!
!$OMP WORKSHARE
         gcons(:,:ntheta2) = gcons(:,:ntheta2) + gcona(:,:ntheta2)
!$OMP END WORKSHARE
!$OMP END PARALLEL

      END IF

      DEALLOCATE (work, gcona)

      END SUBROUTINE alias

      SUBROUTINE alias_new(gcons, ztemp, gcs, gsc, gcc, gss)
      USE vmec_main
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(out) :: gcons
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(in)  :: ztemp
      REAL(rprec), DIMENSION(ns,0:ntor) :: gcs, gsc, gcc, gss
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, i, jk, n, k, js
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work, gcona
!-----------------------------------------------
      ALLOCATE (work(ns*nzeta,4), gcona(ns*nzeta,ntheta3), stat=i)
      IF (i .NE. 0) STOP 'Allocation error in ALIAS'

      gcons = 0

      MLOOP: DO m = 1, mpol1-1
         work = 0
         DO i = 1, ntheta3
            DO jk = 1, ns*nzeta
               work(jk,1) = work(jk,1) + ztemp(jk,i)*cosmui3(i,m)
               work(jk,2) = work(jk,2) + ztemp(jk,i)*sinmui(i,m)
            END DO
         END DO

         gcs = 0;  gsc = 0
         gcc = 0;  gss = 0
         DO n = 0, ntor
            jk=0
            DO k = 1, nzeta
               jk=jk+1        !js starts at 2
               DO js = 2,ns
                  jk = jk+1              
                  gsc(js,n) = gsc(js,n) + tcon(js)*work(jk,2)*
     1               cosnv(k,n)
                  IF (lthreed)
     1               gcs(js,n) = gcs(js,n) + tcon(js)*work(jk,1)*
     1               sinnv(k,n)
                  IF (.NOT.lasym) CYCLE
                  gcc(js,n) = gcc(js,n) + tcon(js)*work(jk,1)*
     1               cosnv(k,n)
                  IF (.NOT.lthreed) CYCLE
                  gss(js,n) = gss(js,n) + tcon(js)*work(jk,2)*
     1               sinnv(k,n)
               END DO
            END DO
         END DO

!
!        INVERSE FOURIER TRANSFORM DE-ALIASED GCON
!
         work = 0

         DO n = 0, ntor
            jk=0
            DO k = 1, nzeta
               jk=jk+1
               DO js = 2, ns
                  jk=jk+1
                  work(jk,4) = work(jk,4) + gsc(js,n)*cosnv(k,n)
                  IF (lthreed)
     1               work(jk,3) = work(jk,3) + gcs(js,n)*sinnv(k,n)
                  IF (.NOT.lasym) CYCLE
                  work(jk,1) = work(jk,1) + gcc(js,n)*cosnv(k,n)
                  IF (.NOT.lthreed) CYCLE
                  work(jk,2) = work(jk,2) + gss(js,n)*sinnv(k,n)
               END DO
            END DO
         END DO

         DO i = 1, ntheta3
            DO jk = 1, ns*nzeta
               gcons(jk,i) = gcons(jk,i) + (work(jk,3)*cosmu(i,m)
     1                     + work(jk,4)*sinmu(i,m))*faccon(m)
               IF (.NOT.lasym) CYCLE
               gcons(jk,i) = gcons(jk,i) + (work(jk,1)*cosmu(i,m)
     1                     + work(jk,2)*sinmu(i,m))*faccon(m)
            END DO
         END DO
      END DO MLOOP

      DEALLOCATE (work, gcona)

      END SUBROUTINE alias_new
