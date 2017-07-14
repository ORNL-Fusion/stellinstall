!*******************************************************************************
!>  @file primed_grid.f
!>  @brief Contains module @ref primed_grid.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref primed_grid_class. This contains the
!>  state variables to define the currents and positions of the volumn integral.
!*******************************************************************************
      MODULE primed_grid
      USE stel_kinds, ONLY: rprec
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  primed grid module parameters
!*******************************************************************************
!>  Bit position for the use curl ju response flag.
      INTEGER, PARAMETER :: primed_grid_ju_flag = 1
!>  Bit position for the use curl jv response flag.
      INTEGER, PARAMETER :: primed_grid_jv_flag = 2

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) primed grid base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a primed grid. This is grid the volume integral will
!>  be summed over.
!-------------------------------------------------------------------------------
      TYPE primed_grid_class
!>  Volume integration element.
         REAL (rprec)                            :: dvol
!>  Toroidal grid spacing.
         REAL (rprec)                            :: dv

!>  X position.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: x => null()
!>  Y position.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: y => null()
!>  Z position.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: z => null()

!>  Current density in the X direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: j_x => null()
!>  Current density in the Y direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: j_y => null()
!>  Current density in the Z direction.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: j_z => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object depending
!>  on the option flags.
!>
!>  @param[in] num_v Number of toroidal grid points.
!>  @param[in] flags Number of toroidal grid points.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct(num_v, flags)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), POINTER :: primed_grid_construct
      INTEGER, INTENT(in)               :: num_v
      INTEGER, INTENT(in)               :: flags

!  local variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (BTEST(flags, primed_grid_ju_flag)) THEN
         primed_grid_construct => primed_grid_construct_ju(num_v)
      ELSE IF (BTEST(flags, primed_grid_jv_flag)) THEN
         primed_grid_construct => primed_grid_construct_jv(num_v)
      ELSE
         primed_grid_construct => primed_grid_construct_both(num_v)
      END IF

      WRITE(*,1000)

      CALL profiler_set_stop_time('primed_grid_construct', start_time)

1000  FORMAT('Prime Grid Ready')

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object with an
!>  This computes the currents and positions on the primed grid. Plasma currents
!>  are obtained from Curl(B).
!>
!>  @param[in] num_v Number of toroidal grid points.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct_both(num_v)
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: xm, xn, mnmax, mnmax_nyq, xm_nyq, xn_nyq,       &
     &                         lasym, xm, xn, rmnc, rmns, zmnc, zmns,          &
     &                         currumnc, currvmnc, currumns, currvmns,         &
     &                         isigng, ns

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), POINTER :: primed_grid_construct_both
      INTEGER, INTENT(in)                     :: num_v

!  local variables
      REAL (rprec)                            :: start_time
      INTEGER                                 :: si
      INTEGER                                 :: ui
      INTEGER                                 :: vi
      REAL (rprec)                            :: u
      REAL (rprec)                            :: v
      REAL (rprec)                            :: r
      REAL (rprec)                            :: z
      REAL (rprec)                            :: ru
      REAL (rprec)                            :: zu
      REAL (rprec)                            :: rv
      REAL (rprec)                            :: zv
      REAL (rprec)                            :: ju
      REAL (rprec)                            :: jv
      REAL (rprec)                            :: jr
      REAL (rprec)                            :: jp
      REAL (rprec)                            :: ds
      REAL (rprec)                            :: cosv
      REAL (rprec)                            :: sinv
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: cosmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: sinmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: cosmn_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: sinmn_nyq

!  local parameters
      INTEGER, PARAMETER                      :: num_u = 101
      REAL (rprec)                            :: du = twopi/num_u

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(primed_grid_construct_both)

      ALLOCATE(primed_grid_construct_both%x(ns,num_u,num_v))
      ALLOCATE(primed_grid_construct_both%y(ns,num_u,num_v))
      ALLOCATE(primed_grid_construct_both%z(ns,num_u,num_v))

      ALLOCATE(primed_grid_construct_both%j_x(ns,num_u,num_v))
      ALLOCATE(primed_grid_construct_both%j_y(ns,num_u,num_v))
      ALLOCATE(primed_grid_construct_both%j_z(ns,num_u,num_v))

      ds = 1.0/(ns - 1.0)
      primed_grid_construct_both%dv = twopi/num_v

      primed_grid_construct_both%dvol = isigng*ds*du                           &
     &   * primed_grid_construct_both%dv/(2.0*twopi)

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(si, ui, vi, r, z, ru, rv, zu, zv, ju, jv, jr, jp,               &
!$OMP&         u, v, cosmn, sinmn, cosmn_nyq, sinmn_nyq, cosv, sinv)

      ALLOCATE(cosmn(mnmax))
      ALLOCATE(sinmn(mnmax))
      ALLOCATE(cosmn_nyq(mnmax_nyq))
      IF (lasym) THEN
         ALLOCATE(sinmn_nyq(mnmax_nyq))
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO vi = 1, num_v
         v = (vi - 0.5)*primed_grid_construct_both%dv

         cosv = COS(v)
         sinv = SIN(v)

         DO ui = 1, num_u
            u = (ui - 0.5)*du

            cosmn = COS(xm*u - xn*v)
            sinmn = SIN(xm*u - xn*v)
            cosmn_nyq = COS(xm_nyq*u - xn_nyq*v)
            IF (lasym) THEN
               sinmn_nyq = SIN(xm_nyq*u - xn_nyq*v)
            END IF

            DO si = 1, ns
               r = SUM(rmnc(:,si)*cosmn(:))
               z = SUM(zmns(:,si)*sinmn(:))

               ru = -SUM(xm*rmnc(:,si)*sinmn(:))
               rv =  SUM(xn*rmnc(:,si)*sinmn(:))
               zu =  SUM(xm*zmns(:,si)*cosmn(:))
               zv = -SUM(xn*zmns(:,si)*cosmn(:))

               ju = SUM(currumnc(:,si)*cosmn_nyq(:))*mu0
               jv = SUM(currvmnc(:,si)*cosmn_nyq(:))*mu0

               IF (lasym) THEN
                  r = r + SUM(rmns(:,si)*sinmn(:))
                  z = z + SUM(zmnc(:,si)*cosmn(:))

                  ru = ru + SUM(xm*rmns(:,si)*cosmn(:))
                  rv = rv - SUM(xn*rmns(:,si)*cosmn(:))
                  zu = zu - SUM(xm*zmnc(:,si)*sinmn(:))
                  zv = zv + SUM(xn*zmnc(:,si)*sinmn(:))

                  ju = ju + SUM(currumns(:,si)*sinmn_nyq(:))*mu0
                  jv = jv + SUM(currvmns(:,si)*sinmn_nyq(:))*mu0
               END IF

               jr = ju*ru + jv*rv
               jp = jv*r
               primed_grid_construct_both%j_z(si,ui,vi) = ju*zu + jv*zv

               IF (si .eq. 1 .or. si .eq. ns) THEN
                  jr = jr/2.0
                  jp = jp/2.0
                  primed_grid_construct_both%j_z(si,ui,vi) =                   &
     &               primed_grid_construct_both%j_z(si,ui,vi)/2.0
               END IF

               primed_grid_construct_both%j_x(si,ui,vi) = jr*cosv              &
     &                                                  - jp*sinv
               primed_grid_construct_both%j_y(si,ui,vi) = jr*sinv              &
     &                                                  + jp*cosv

               primed_grid_construct_both%x(si,ui,vi) = r*cosv
               primed_grid_construct_both%y(si,ui,vi) = r*sinv
               primed_grid_construct_both%z(si,ui,vi) = z
            END DO
         END DO
      END DO
!$OMP END DO

      DEALLOCATE(cosmn)
      DEALLOCATE(sinmn)
      DEALLOCATE(cosmn_nyq)
      IF (lasym) THEN
         DEALLOCATE(sinmn_nyq)
      END IF
!$OMP END PARALLEL

      CALL profiler_set_stop_time('primed_grid_construct_both',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object with an
!>  This computes the currents and positions on the primed grid. J^u plasma
!>  current is obtained from force balance.
!>
!>  J^u = (p' + J^v*B^u)/B^v
!>
!>  @param[in] num_v Number of toroidal grid points.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct_ju(num_v)
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: xm, xn, mnmax, mnmax_nyq, xm_nyq, xn_nyq,       &
     &                         lasym, xm, xn, rmnc, rmns, zmnc, zmns,          &
     &                         currvmnc, currvmns,                             &
     &                         bsupumnc, bsupvmnc, bsupumns, bsupvmns,         &
     &                         isigng, ns, presf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), POINTER :: primed_grid_construct_ju
      INTEGER, INTENT(in)                       :: num_v

!  local variables
      REAL (rprec)                              :: start_time
      INTEGER                                   :: si
      INTEGER                                   :: ui
      INTEGER                                   :: vi
      REAL (rprec)                              :: u
      REAL (rprec)                              :: v
      REAL (rprec)                              :: r
      REAL (rprec)                              :: z
      REAL (rprec)                              :: ru
      REAL (rprec)                              :: zu
      REAL (rprec)                              :: rv
      REAL (rprec)                              :: zv
      REAL (rprec)                              :: bu
      REAL (rprec)                              :: bv
      REAL (rprec)                              :: ju
      REAL (rprec)                              :: jv
      REAL (rprec)                              :: jr
      REAL (rprec)                              :: jp
      REAL (rprec)                              :: ds
      REAL (rprec)                              :: cosv
      REAL (rprec)                              :: sinv
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: sinmn
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: rmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: rmnsh
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: zmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: zmnsh
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: currvmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: currvmnsh
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: p_prime
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosmn_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: sinmn_nyq

!  local parameters
      INTEGER, PARAMETER                        :: num_u = 101
      REAL (rprec)                              :: du = twopi/num_u

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(primed_grid_construct_ju)

      ALLOCATE(primed_grid_construct_ju%x(ns - 1,num_u,num_v))
      ALLOCATE(primed_grid_construct_ju%y(ns - 1,num_u,num_v))
      ALLOCATE(primed_grid_construct_ju%z(ns - 1,num_u,num_v))

      ALLOCATE(primed_grid_construct_ju%j_x(ns - 1,num_u,num_v))
      ALLOCATE(primed_grid_construct_ju%j_y(ns - 1,num_u,num_v))
      ALLOCATE(primed_grid_construct_ju%j_z(ns - 1,num_u,num_v))

      ALLOCATE(rmnch(mnmax,ns - 1))
      ALLOCATE(zmnsh(mnmax,ns - 1))
      ALLOCATE(currvmnch(mnmax_nyq,ns - 1))

      IF (lasym) THEN
         ALLOCATE(rmnsh(mnmax,ns - 1))
         ALLOCATE(zmnch(mnmax,ns - 1))
         ALLOCATE(currvmnsh(mnmax_nyq,ns - 1))
      END IF

      ALLOCATE(p_prime(ns - 1))

      ds = 1.0/(ns - 1.0)
      primed_grid_construct_ju%dv = twopi/num_v

      primed_grid_construct_ju%dvol = isigng*ds*du                             &
     &   * primed_grid_construct_ju%dv/(2.0*twopi)

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(si, ui, vi, r, z, ru, rv, zu, zv, ju, jv, jr, jp,               &
!$OMP&         u, v, cosmn, sinmn, cosmn_nyq, sinmn_nyq, cosv, sinv,           &
!$OMP&         bu, bv)

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO si = 1, ns - 1
         rmnch(:,si) = (rmnc(:,si + 1) + rmnc(:,si))/2.0
         zmnsh(:,si) = (zmns(:,si + 1) + zmns(:,si))/2.0

         currvmnch(:,si) = (currvmnc(:,si + 1) + currvmnc(:,si))/2.0

         p_prime(si) = (presf(si + 1) - presf(si))/ds

         IF (lasym) THEN
            rmnsh(:,si) = (rmns(:,si + 1) + rmns(:,si))/2.0
            zmnch(:,si) = (zmnc(:,si + 1) + zmnc(:,si))/2.0

            currvmnsh(:,si) = (currvmns(:,si + 1) + currvmns(:,si))/2.0
         END IF
      END DO
!$OMP END DO

      ALLOCATE(cosmn(mnmax))
      ALLOCATE(sinmn(mnmax))
      ALLOCATE(cosmn_nyq(mnmax_nyq))
      IF (lasym) THEN
         ALLOCATE(sinmn_nyq(mnmax_nyq))
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO vi = 1, num_v
         v = (vi - 0.5)*primed_grid_construct_ju%dv

         cosv = COS(v)
         sinv = SIN(v)

         DO ui = 1, num_u
            u = (ui - 0.5)*du

            cosmn = COS(xm*u - xn*v)
            sinmn = SIN(xm*u - xn*v)
            cosmn_nyq = COS(xm_nyq*u - xn_nyq*v)
            IF (lasym) THEN
               sinmn_nyq = SIN(xm_nyq*u - xn_nyq*v)
            END IF

            DO si = 1, ns - 1
               r = SUM(rmnch(:,si)*cosmn(:))
               z = SUM(zmnsh(:,si)*sinmn(:))

               ru = -SUM(xm*rmnch(:,si)*sinmn(:))
               rv =  SUM(xn*rmnch(:,si)*sinmn(:))
               zu =  SUM(xm*zmnsh(:,si)*cosmn(:))
               zv = -SUM(xn*zmnsh(:,si)*cosmn(:))

               jv = SUM(currvmnch(:,si)*cosmn_nyq(:))

               bu = SUM(bsupumnc(:,si + 1)*cosmn_nyq(:))
               bv = SUM(bsupvmnc(:,si + 1)*cosmn_nyq(:))

               IF (lasym) THEN
                  r = r + SUM(rmnsh(:,si)*sinmn(:))
                  z = z + SUM(zmnch(:,si)*cosmn(:))

                  ru = ru + SUM(xm*rmnsh(:,si)*cosmn(:))
                  rv = rv - SUM(xn*rmnsh(:,si)*cosmn(:))
                  zu = zu - SUM(xm*zmnch(:,si)*sinmn(:))
                  zv = zv + SUM(xn*zmnch(:,si)*sinmn(:))

                  jv = jv + SUM(currvmns(:,si)*sinmn_nyq(:))

                  bu = bu + SUM(bsupumns(:,si + 1)*sinmn_nyq(:))
                  bv = bv + SUM(bsupvmns(:,si + 1)*sinmn_nyq(:))
               END IF

               ju = (p_prime(si) + jv*bu)/bv*mu0
               jv = jv*mu0

               jr = ju*ru + jv*rv
               jp = jv*r
               primed_grid_construct_ju%j_z(si,ui,vi) = ju*zu + jv*zv

               primed_grid_construct_ju%j_x(si,ui,vi) = jr*cosv                &
     &                                                - jp*sinv
               primed_grid_construct_ju%j_y(si,ui,vi) = jr*sinv                &
     &                                                + jp*cosv

               primed_grid_construct_ju%x(si,ui,vi) = r*cosv
               primed_grid_construct_ju%y(si,ui,vi) = r*sinv
               primed_grid_construct_ju%z(si,ui,vi) = z
            END DO
         END DO
      END DO
!$OMP END DO

      DEALLOCATE(cosmn)
      DEALLOCATE(sinmn)
      DEALLOCATE(cosmn_nyq)
      IF (lasym) THEN
         DEALLOCATE(sinmn_nyq)
      END IF
!$OMP END PARALLEL

      DEALLOCATE(rmnch)
      DEALLOCATE(zmnsh)
      DEALLOCATE(currvmnch)
      IF (lasym) THEN
         DEALLOCATE(rmnsh)
         DEALLOCATE(zmnch)
         DEALLOCATE(currvmnsh)
      END IF

      CALL profiler_set_stop_time('primed_grid_construct_ju',                  &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object with an
!>  This computes the currents and positions on the primed grid. J^v plasma
!>  current is obtained from force balance.
!>
!>  J^v = (J^u*B^v - p')/B^u
!>
!>  @param[in] num_v Number of toroidal grid points.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct_jv(num_v)
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: xm, xn, mnmax, mnmax_nyq, xm_nyq, xn_nyq,       &
     &                         lasym, xm, xn, rmnc, rmns, zmnc, zmns,          &
     &                         currumnc, currumns,                             &
     &                         bsupumnc, bsupvmnc, bsupumns, bsupvmns,         &
     &                         isigng, ns, presf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), POINTER :: primed_grid_construct_jv
      INTEGER, INTENT(in)                       :: num_v

!  local variables
      REAL (rprec)                              :: start_time
      INTEGER                                   :: si
      INTEGER                                   :: ui
      INTEGER                                   :: vi
      REAL (rprec)                              :: u
      REAL (rprec)                              :: v
      REAL (rprec)                              :: r
      REAL (rprec)                              :: z
      REAL (rprec)                              :: ru
      REAL (rprec)                              :: zu
      REAL (rprec)                              :: rv
      REAL (rprec)                              :: zv
      REAL (rprec)                              :: bu
      REAL (rprec)                              :: bv
      REAL (rprec)                              :: ju
      REAL (rprec)                              :: jv
      REAL (rprec)                              :: jr
      REAL (rprec)                              :: jp
      REAL (rprec)                              :: ds
      REAL (rprec)                              :: cosv
      REAL (rprec)                              :: sinv
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: sinmn
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: rmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: rmnsh
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: zmnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: zmnsh
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: currumnch
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: currumnsh
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: p_prime
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosmn_nyq
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: sinmn_nyq

!  local parameters
      INTEGER, PARAMETER                        :: num_u = 101
      REAL (rprec)                              :: du = twopi/num_u

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(primed_grid_construct_jv)

      ALLOCATE(primed_grid_construct_jv%x(ns - 1,num_u,num_v))
      ALLOCATE(primed_grid_construct_jv%y(ns - 1,num_u,num_v))
      ALLOCATE(primed_grid_construct_jv%z(ns - 1,num_u,num_v))

      ALLOCATE(primed_grid_construct_jv%j_x(ns - 1,num_u,num_v))
      ALLOCATE(primed_grid_construct_jv%j_y(ns - 1,num_u,num_v))
      ALLOCATE(primed_grid_construct_jv%j_z(ns - 1,num_u,num_v))

      ALLOCATE(rmnch(mnmax,ns - 1))
      ALLOCATE(zmnsh(mnmax,ns - 1))
      ALLOCATE(currumnch(mnmax_nyq,ns - 1))

      IF (lasym) THEN
         ALLOCATE(rmnsh(mnmax,ns - 1))
         ALLOCATE(zmnch(mnmax,ns - 1))
         ALLOCATE(currumnsh(mnmax_nyq,ns - 1))
      END IF

      ALLOCATE(p_prime(ns - 1))

      ds = 1.0/(ns - 1.0)
      primed_grid_construct_jv%dv = twopi/num_v

      primed_grid_construct_jv%dvol = isigng*ds*du                             &
     &   * primed_grid_construct_jv%dv/(2.0*twopi)

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(si, ui, vi, r, z, ru, rv, zu, zv, ju, jv, jr, jp,               &
!$OMP&         u, v, cosmn, sinmn, cosmn_nyq, sinmn_nyq, cosv, sinv,           &
!$OMP&         bu, bv)

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO si = 1, ns - 1
         rmnch(:,si) = (rmnc(:,si + 1) + rmnc(:,si))/2.0
         zmnsh(:,si) = (zmns(:,si + 1) + zmns(:,si))/2.0

         currumnch(:,si) = (currumnc(:,si + 1) + currumnc(:,si))/2.0

         p_prime(si) = (presf(si + 1) - presf(si))/ds

         IF (lasym) THEN
            rmnsh(:,si) = (rmns(:,si + 1) + rmns(:,si))/2.0
            zmnch(:,si) = (zmnc(:,si + 1) + zmnc(:,si))/2.0

            currumnsh(:,si) = (currumns(:,si + 1) + currumns(:,si))/2.0
         END IF
      END DO
!$OMP END DO

      ALLOCATE(cosmn(mnmax))
      ALLOCATE(sinmn(mnmax))
      ALLOCATE(cosmn_nyq(mnmax_nyq))
      IF (lasym) THEN
         ALLOCATE(sinmn_nyq(mnmax_nyq))
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO vi = 1, num_v
         v = (vi - 0.5)*primed_grid_construct_jv%dv

         cosv = COS(v)
         sinv = SIN(v)

         DO ui = 1, num_u
            u = (ui - 0.5)*du

            cosmn = COS(xm*u - xn*v)
            sinmn = SIN(xm*u - xn*v)
            cosmn_nyq = COS(xm_nyq*u - xn_nyq*v)
            IF (lasym) THEN
               sinmn_nyq = SIN(xm_nyq*u - xn_nyq*v)
            END IF

            DO si = 1, ns - 1
               r = SUM(rmnch(:,si)*cosmn(:))
               z = SUM(zmnsh(:,si)*sinmn(:))

               ru = -SUM(xm*rmnch(:,si)*sinmn(:))
               rv =  SUM(xn*rmnch(:,si)*sinmn(:))
               zu =  SUM(xm*zmnsh(:,si)*cosmn(:))
               zv = -SUM(xn*zmnsh(:,si)*cosmn(:))

               jv = SUM(currumnch(:,si)*cosmn_nyq(:))

               bu = SUM(bsupumnc(:,si + 1)*cosmn_nyq(:))
               bv = SUM(bsupvmnc(:,si + 1)*cosmn_nyq(:))

               IF (lasym) THEN
                  r = r + SUM(rmnsh(:,si)*sinmn(:))
                  z = z + SUM(zmnch(:,si)*cosmn(:))

                  ru = ru + SUM(xm*rmnsh(:,si)*cosmn(:))
                  rv = rv - SUM(xn*rmnsh(:,si)*cosmn(:))
                  zu = zu - SUM(xm*zmnch(:,si)*sinmn(:))
                  zv = zv + SUM(xn*zmnch(:,si)*sinmn(:))

                  ju = ju + SUM(currumns(:,si)*sinmn_nyq(:))

                  bu = bu + SUM(bsupumns(:,si + 1)*sinmn_nyq(:))
                  bv = bv + SUM(bsupvmns(:,si + 1)*sinmn_nyq(:))
               END IF

               jv = (ju*bv - p_prime(si))/bu*mu0
               ju = ju*mu0

               jr = ju*ru + jv*rv
               jp = jv*r
               primed_grid_construct_jv%j_z(si,ui,vi) = ju*zu + jv*zv

               primed_grid_construct_jv%j_x(si,ui,vi) = jr*cosv                &
     &                                                - jp*sinv
               primed_grid_construct_jv%j_y(si,ui,vi) = jr*sinv                &
     &                                                + jp*cosv

               primed_grid_construct_jv%x(si,ui,vi) = r*cosv
               primed_grid_construct_jv%y(si,ui,vi) = r*sinv
               primed_grid_construct_jv%z(si,ui,vi) = z
            END DO
         END DO
      END DO
!$OMP END DO

      DEALLOCATE(cosmn)
      DEALLOCATE(sinmn)
      DEALLOCATE(cosmn_nyq)
      IF (lasym) THEN
         DEALLOCATE(sinmn_nyq)
      END IF
!$OMP END PARALLEL

      DEALLOCATE(rmnch)
      DEALLOCATE(zmnsh)
      DEALLOCATE(currumnch)
      IF (lasym) THEN
         DEALLOCATE(rmnsh)
         DEALLOCATE(zmnch)
         DEALLOCATE(currumnsh)
      END IF

      CALL profiler_set_stop_time('primed_grid_construct_jv',                  &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref primed_grid_class object.
!>
!>  Deallocates memory and uninitializes a @ref primed_grid_class object.
!>
!>  @param[inout] this A @ref primed_grid_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE primed_grid_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%x)) THEN
         DEALLOCATE(this%x)
         this%x => null()
      END IF

      IF (ASSOCIATED(this%y)) THEN
         DEALLOCATE(this%y)
         this%y => null()
      END IF

      IF (ASSOCIATED(this%z)) THEN
         DEALLOCATE(this%z)
         this%z => null()
      END IF

      IF (ASSOCIATED(this%j_x)) THEN
         DEALLOCATE(this%j_x)
         this%j_x => null()
      END IF

      IF (ASSOCIATED(this%j_y)) THEN
         DEALLOCATE(this%j_y)
         this%j_y => null()
      END IF

      IF (ASSOCIATED(this%j_z)) THEN
         DEALLOCATE(this%j_z)
         this%j_z => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

      END MODULE
