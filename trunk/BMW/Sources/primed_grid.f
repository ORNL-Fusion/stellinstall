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
!>  Bit position for the use siesta instead of VMEC.
      INTEGER, PARAMETER :: primed_grid_siesta_flag = 3

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
      FUNCTION primed_grid_construct(num_v, flags, siesta_file)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), POINTER :: primed_grid_construct
      INTEGER, INTENT(in)               :: num_v
      INTEGER, INTENT(in)               :: flags
      CHARACTER (len=*), INTENT(in)     :: siesta_file

!  local variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (BTEST(flags, primed_grid_ju_flag)) THEN
         primed_grid_construct => primed_grid_construct_ju(num_v)
      ELSE IF (BTEST(flags, primed_grid_jv_flag)) THEN
         primed_grid_construct => primed_grid_construct_jv(num_v)
      ELSE IF (BTEST(flags, primed_grid_siesta_flag)) THEN
         primed_grid_construct =>                                              &
     &      primed_grid_construct_siesta(num_v, siesta_file)
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
      USE read_wout_mod, ONLY: mnmax, mnmax_nyq, xm, xn, xm_nyq, xn_nyq,       &
     &                         lasym, ns, isigng, rmnc, rmns,                  &
     &                         zmnc, zmns, currumnc, currvmnc,                 &
     &                         currumns, currvmns

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

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref primed_grid_class object.
!>
!>  Allocates memory and initializes a @ref primed_grid_class object with an
!>  This computes the currents and positions on the primed grid. Plasma currents
!>  are obtained from Curl(B) of the siesta solution.
!>
!>  @param[in] num_v Number of toroidal grid points.
!>  @returns A pointer to a constructed @ref primed_grid_class object.
!-------------------------------------------------------------------------------
      FUNCTION primed_grid_construct_siesta(num_v, siesta_file)
      USE stel_constants, ONLY: twopi, mu0
      USE read_wout_mod, ONLY: nfp, ns, rmnc, rmns, zmnc, zmns, isigng,        &
     &                         mnmax, xm, xn, lasym
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (primed_grid_class), POINTER :: primed_grid_construct_siesta
      INTEGER, INTENT(in)                         :: num_v
      CHARACTER (len=*)                           :: siesta_file

!  local variables
      REAL (rprec)                                :: start_time
      INTEGER                                     :: si
      INTEGER                                     :: ui
      INTEGER                                     :: vi
      REAL (rprec)                                :: s_v
      REAL (rprec)                                :: s
      REAL (rprec)                                :: u
      REAL (rprec)                                :: v
      REAL (rprec)                                :: r
      REAL (rprec)                                :: z
      REAL (rprec)                                :: rs
      REAL (rprec)                                :: zs
      REAL (rprec)                                :: ru
      REAL (rprec)                                :: zu
      REAL (rprec)                                :: rv
      REAL (rprec)                                :: zv
      REAL (rprec)                                :: js
      REAL (rprec)                                :: ju
      REAL (rprec)                                :: jv
      REAL (rprec)                                :: jr
      REAL (rprec)                                :: jp
      REAL (rprec)                                :: jz
      REAL (rprec)                                :: ds
      REAL (rprec)                                :: cosv
      REAL (rprec)                                :: sinv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE   :: cosmn
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE   :: sinmn
      REAL (rprec), DIMENSION(:), ALLOCATABLE     :: cosmn_vmec
      REAL (rprec), DIMENSION(:), ALLOCATABLE     :: sinmn_vmec
      INTEGER                                     :: ncid
      INTEGER                                     :: status
      INTEGER                                     :: nrad
      INTEGER                                     :: mpol
      INTEGER                                     :: ntor
      INTEGER                                     :: m
      INTEGER                                     :: n
      INTEGER                                     :: flags
      REAL (rprec)                                :: b_factor
      REAL (rprec)                                :: wlow_r
      INTEGER                                     :: ilow
      INTEGER                                     :: ihigh
      REAL (rprec)                                :: wlow
      REAL (rprec)                                :: whigh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupsmncf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupsmnsf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupumncf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupumnsf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupvmncf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupvmnsf
      REAL (rprec), DIMENSION(:), ALLOCATABLE     :: rmncint
      REAL (rprec), DIMENSION(:), ALLOCATABLE     :: rmnsint
      REAL (rprec), DIMENSION(:), ALLOCATABLE     :: zmncint
      REAL (rprec), DIMENSION(:), ALLOCATABLE     :: zmnsint

!  local parameters
      INTEGER, PARAMETER                          :: num_u = 101
      REAL (rprec)                                :: du = twopi/num_u
      INTEGER, PARAMETER                          :: lasym_f = 31

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL cdf_open(ncid, TRIM(siesta_file), 'r', status)

      CALL cdf_read(ncid, 'nrad', nrad)
      CALL cdf_read(ncid, 'ntor', ntor)
      CALL cdf_read(ncid, 'mpol', mpol)

      CALL cdf_read(ncid, 'state_flags', flags)

      CALL cdf_read(ncid, 'b_factor', b_factor)

      ALLOCATE(jksupsmnsf(0:mpol,-ntor:ntor,nrad))
      ALLOCATE(jksupumncf(0:mpol,-ntor:ntor,nrad))
      ALLOCATE(jksupvmncf(0:mpol,-ntor:ntor,nrad))

      CALL cdf_read(ncid, 'jksupsmnsf(m,n,r)', jksupsmnsf)
      CALL cdf_read(ncid, 'jksupumncf(m,n,r)', jksupumncf)
      CALL cdf_read(ncid, 'jksupvmncf(m,n,r)', jksupvmncf)

      IF (BTEST(flags, lasym_f)) THEN
         ALLOCATE(jksupsmncf(0:mpol,-ntor:ntor,nrad))
         ALLOCATE(jksupumnsf(0:mpol,-ntor:ntor,nrad))
         ALLOCATE(jksupvmnsf(0:mpol,-ntor:ntor,nrad))

         CALL cdf_read(ncid, 'jksupsmncf(m,n,r)', jksupsmncf)
         CALL cdf_read(ncid, 'jksupumnsf(m,n,r)', jksupumnsf)
         CALL cdf_read(ncid, 'jksupvmnsf(m,n,r)', jksupvmnsf)
      END IF

      CALL cdf_close(ncid)

      ALLOCATE(primed_grid_construct_siesta)

      ALLOCATE(primed_grid_construct_siesta%x(nrad,num_u,num_v))
      ALLOCATE(primed_grid_construct_siesta%y(nrad,num_u,num_v))
      ALLOCATE(primed_grid_construct_siesta%z(nrad,num_u,num_v))

      ALLOCATE(primed_grid_construct_siesta%j_x(nrad,num_u,num_v))
      ALLOCATE(primed_grid_construct_siesta%j_y(nrad,num_u,num_v))
      ALLOCATE(primed_grid_construct_siesta%j_z(nrad,num_u,num_v))

      ds = 1.0/(ns - 1.0)
      primed_grid_construct_siesta%dv = twopi/num_v

      primed_grid_construct_siesta%dvol = isigng*ds*du                         &
     &   * primed_grid_construct_siesta%dv/(2.0*twopi)

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(si, ui, vi, r, z, rs, ru, rv, zs, zu, zv, js, ju, jv,           &
!$OMP&         jr, jp, s, u, v, cosmn, sinmn, cosmn_vmec, sinmn_vmec,          &
!$OMP&         cosv, sinv, ilow, ihigh, wlow, whigh, rmncint, zmnsint,         &
!$OMP&         rmnsint, zmncint)

      ALLOCATE(cosmn(0:mpol,-ntor:ntor))
      ALLOCATE(sinmn(0:mpol,-ntor:ntor))

      ALLOCATE(cosmn_vmec(mnmax))
      ALLOCATE(sinmn_vmec(mnmax))

      ALLOCATE(rmncint(mnmax))
      ALLOCATE(zmnsint(mnmax))
      IF (BTEST(flags, lasym_f)) THEN
         ALLOCATE(rmnsint(mnmax))
         ALLOCATE(zmncint(mnmax))
      END IF

!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO vi = 1, num_v
         v = (vi - 0.5)*primed_grid_construct_siesta%dv

         cosv = COS(v)
         sinv = SIN(v)

         DO ui = 1, num_u
            u = (ui - 0.5)*du

            cosmn_vmec = COS(xm*u - xn*v)
            sinmn_vmec = SIN(xm*u - xn*v)

            DO n = -ntor, ntor
               DO m = 0, mpol
                  cosmn(m,n) = COS(m*u + nfp*n*v)
                  sinmn(m,n) = SIN(m*u + nfp*n*v)
               END DO
            END DO

            DO si = 1, nrad
               s = (si - 1.0)/(nrad - 1.0)
               s_v = s*s

               wlow_r = s_v*(ns - 1.0) + 1.0
               ilow = FLOOR(wlow_r)
               IF (ilow .eq. ns) THEN
                  ilow = ns - 1
               END IF
               ihigh = ilow + 1
               wlow = -wlow_r + 1.0 + ilow
               whigh = 1.0 - wlow

               rmncint = wlow*rmnc(:,ilow) + whigh*rmnc(:,ihigh)
               zmnsint = wlow*zmns(:,ilow) + whigh*zmns(:,ihigh)

               r = SUM(rmncint*cosmn_vmec)
               z = SUM(zmnsint*sinmn_vmec)

               ru = -SUM(xm*rmncint*sinmn_vmec)
               rv =  SUM(xn*rmncint*sinmn_vmec)
               zu =  SUM(xm*zmnsint*cosmn_vmec)
               zv = -SUM(xn*zmnsint*cosmn_vmec)

               rmncint = 2.0*s*(rmnc(:,ihigh) - rmnc(:,ilow))*(ns - 1.0)
               zmnsint = 2.0*s*(zmns(:,ihigh) - zmns(:,ilow))*(ns - 1.0)

               rs = SUM(rmncint*cosmn_vmec)
               zs = SUM(zmnsint*sinmn_vmec)

               js = SUM(jksupsmnsf(:,:,si)*sinmn)
               ju = SUM(jksupumncf(:,:,si)*cosmn)
               jv = SUM(jksupvmncf(:,:,si)*cosmn)

               IF (BTEST(flags, lasym_f)) THEN
                  IF (lasym) THEN
                     rmnsint = wlow*rmns(:,ilow) + whigh*rmns(:,ihigh)
                     zmncint = wlow*zmnc(:,ilow) + whigh*zmnc(:,ihigh)

                     r = r + SUM(rmnsint*sinmn_vmec)
                     z = z + SUM(zmncint*cosmn_vmec)

                     ru = ru + SUM(xm*rmnsint*cosmn_vmec)
                     rv = rv - SUM(xn*rmnsint*cosmn_vmec)
                     zu = zu - SUM(xm*zmncint*sinmn_vmec)
                     zv = zv + SUM(xn*zmncint*sinmn_vmec)

                     rmnsint = 2.0*s*(rmns(:,ihigh) - rmns(:,ilow))            &
     &                       * (ns - 1.0)
                     zmncint = 2.0*s*(zmnc(:,ihigh) - zmnc(:,ilow))            &
     &                       * (ns - 1.0)

                     rs = rs + SUM(rmnsint*sinmn_vmec)
                     zs = zs + SUM(zmncint*cosmn_vmec)
                  END IF

                  js = js + SUM(jksupsmncf(:,:,si)*cosmn)
                  ju = ju + SUM(jksupumnsf(:,:,si)*sinmn)
                  jv = jv + SUM(jksupvmnsf(:,:,si)*sinmn)

               END IF

               js = js/(b_factor*mu0)
               ju = ju/(b_factor*mu0)
               jv = jv/(b_factor*mu0)

               jr = js*rs + ju*ru + jv*rv
               jp = jv*r
               jz = ju*zs + ju*zu + jv*zv

               IF (si .eq. 1 .or. si .eq. ns) THEN
                  jr = jr/2.0
                  jp = jp/2.0
                  jz = jz/2.0
               END IF

               primed_grid_construct_siesta%j_x(si,ui,vi) = jr*cosv            &
     &                                                    - jp*sinv
               primed_grid_construct_siesta%j_y(si,ui,vi) = jr*sinv            &
     &                                                    + jp*cosv
               primed_grid_construct_siesta%j_z(si,ui,vi) = jz

               primed_grid_construct_siesta%x(si,ui,vi) = r*cosv
               primed_grid_construct_siesta%y(si,ui,vi) = r*sinv
               primed_grid_construct_siesta%z(si,ui,vi) = z
            END DO
         END DO
       END DO
!$OMP END DO

       DEALLOCATE(cosmn)
       DEALLOCATE(sinmn)

       DEALLOCATE(cosmn_vmec)
       DEALLOCATE(sinmn_vmec)

       DEALLOCATE(rmncint)
       DEALLOCATE(zmnsint)
       IF (BTEST(flags, lasym_f)) THEN
          DEALLOCATE(rmnsint)
          DEALLOCATE(zmncint)
       END IF

!$OMP END PARALLEL

      CALL profiler_set_stop_time('primed_grid_construct_siesta',              &
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
