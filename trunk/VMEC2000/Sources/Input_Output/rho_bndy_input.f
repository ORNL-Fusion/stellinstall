!*******************************************************************************
!>  @file rho_bndy_input.f
!>  @brief Convert between rho boundary and r,z boundary representation.
!
!  Note separating the Doxygen comment block here so detailed description is
!  found in the Module not the file.
!
!>  Defines routines to change between rho and r,z representation.
!*******************************************************************************
      MODULE rho_bndy
      USE stel_kinds, ONLY: rprec

      IMPLICIT NONE

!*******************************************************************************
!  rho_bndy module parameters
!*******************************************************************************
!>  Stellarator Symmetric boundary flag.
      LOGICAL, PARAMETER :: rho_bndy_sym = .TRUE.
!>  Stellaratir Asymmetric boundary flag.
      LOGICAL, PARAMETER :: rho_bndy_asym = .FALSE.

!*******************************************************************************
!  Derived-type declarations
!  1) rho_bndy_class
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Class containing coeffients to convert between rho and r,z.
!-------------------------------------------------------------------------------
      TYPE :: rho_bndy_class
!>  t1 poloidal angle coefficent.
         REAL (rprec), DIMENSION(:), POINTER :: t1 => null()
!>  t2 poloidal angle coefficent.
         REAL (rprec), DIMENSION(:), POINTER :: t2 => null()
      CONTAINS
         PROCEDURE                           :: to_rz => rho_bndy_to_rz
         PROCEDURE                           ::                                &
     &      to_rho => rho_bndy_to_rho
         FINAL                               :: rho_bndy_destruct
      END TYPE

!*******************************************************************************
!>  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for rho_bndy constructor.
!-------------------------------------------------------------------------------
      INTERFACE rho_bndy_class
         MODULE PROCEDURE rho_bndy_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for checking results of the unit tests.
!-------------------------------------------------------------------------------
      INTERFACE check
         MODULE PROCEDURE check_real,                                          &
     &                    check_fourier_1d,                                    &
     &                    check_fourier_2d
      END INTERFACE

      PRIVATE :: check, check_real, check_fourier_1d, check_fourier_2d

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Setup t1 and t2 coefficents.
!>
!>  @param[in]  mpol Number of poloidal modes.
!>  @param[out] t1m  t1 angle coefficent.
!>  @param[out] t2m  t2 angle coefficent.
!>  @returns A pointer to a constructed @ref rho_bndy_class object.
!-------------------------------------------------------------------------------
      FUNCTION rho_bndy_construct(mpol)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (rho_bndy_class), POINTER :: rho_bndy_construct
      INTEGER, INTENT(in)             :: mpol

!  Local Variables
      INTEGER                         :: m

!  Local Parameters
      INTEGER, PARAMETER              :: mexp = 4

!  Start of executable code
      ALLOCATE(rho_bndy_construct)
      ALLOCATE(rho_bndy_construct%t1(0:mpol+1))
      ALLOCATE(rho_bndy_construct%t2(0:mpol+1))

      rho_bndy_construct%t2(0) = 0.5**mexp
      rho_bndy_construct%t1(0) = 0.0
      rho_bndy_construct%t1(1) = 1.0
!  Better fit if axis is fixed. Minics polar m=1 coupling.
      rho_bndy_construct%t1(2) = rho_bndy_construct%t2(0)

      DO m = 1, mpol + 1
         rho_bndy_construct%t2(m) = (REAL(m + 1)/m)**mexp
      END DO
      DO m = 3, mpol + 1
         rho_bndy_construct%t1(m) = (REAL(m - 1)/m)**mexp
      END DO

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Destruct a @ref rho_bndy_class object.
!>
!>  Deallocates memory and uninitializes a @ref rho_bndy_class object.
!>
!>  @param[inout] this A @ref rho_bndy_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE rho_bndy_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (rho_bndy_class), INTENT(inout) :: this

!  Start of executable code.
      IF (ASSOCIATED(this%t1)) THEN
         DEALLOCATE(this%t1)
         this%t1 => null()
      END IF

      IF (ASSOCIATED(this%t2)) THEN
         DEALLOCATE(this%t2)
         this%t2 => null()
      END IF

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Convert from rho to r,z
!>
!>  This routine uses the rho Fourier components to compute the r,z fourier
!>  components. Assuming a fixed center point,
!>
!>  @param[in]  this   A @ref rho_bndy_class instance.
!>  @param[in]  mpol   Number of poloidal modes.
!>  @param[in]  ntor   Number of toroidal modes.
!>  @param[in]  rho    Rho boundary coefficient.
!>  @param[in]  r_ctr  Center R bounday cofficient.
!>  @param[in]  z_ctr  Center R bounday cofficient.
!>  @param[out] r_bndy R boundary coefficient.
!>  @param[out] z_bndy Z boundary coefficient.
!>  @param[in]  parity Stellarator symmetric or asymmetric parity.
!-------------------------------------------------------------------------------
      SUBROUTINE rho_bndy_to_rz(this, mpol, ntor, rho, r_ctr, z_ctr,           &
     &                          r_bndy, z_bndy, parity)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (rho_bndy_class), INTENT(in)                        :: this
      INTEGER, INTENT(in)                                       :: ntor
      INTEGER, INTENT(in)                                       :: mpol
      REAL (rprec), DIMENSION(-ntor:ntor,0:mpol-1), INTENT(in)  :: rho
      REAL (rprec), DIMENSION(0:ntor), INTENT(in)               :: r_ctr
      REAL (rprec), DIMENSION(0:ntor), INTENT(in)               :: z_ctr
      REAL (rprec), DIMENSION(-ntor:ntor,0:mpol-1), INTENT(out) ::             &
     &   r_bndy
      REAL (rprec), DIMENSION(-ntor:ntor,0:mpol-1), INTENT(out) ::             &
     &   z_bndy
      LOGICAL, INTENT(in)                                       ::             &
     &   parity

!  Local Variables
      INTEGER                                                   :: m
      INTEGER                                                   :: n
      REAL (rprec)                                              :: tnorm
      REAL (rprec)                                              :: x
      INTEGER                                                   :: sign

!  Start of executable code.
      r_bndy = 0
      z_bndy = 0

!  m = 0
      IF (parity .eqv. rho_bndy_sym) THEN
         sign = 1
         r_bndy(0,0) = this%t2(0)*rho(0,1) + r_ctr(0)
      ELSE
         sign = -1
         z_bndy(0,0) = this%t2(0)*rho(0,1) + z_ctr(0)
      END IF

      DO n = 1, ntor
         r_bndy(n,0) =       this%t2(0)*rho( n,1)                              &
     &               +  sign*this%t2(0)*rho(-n,1)                              &
     &               + r_ctr(n)
         z_bndy(n,0) = -sign*this%t2(0)*rho( n,1)                              &
     &               +       this%t2(0)*rho(-n,1)                              &
     &               + z_ctr(n)
      END DO

!  m > 0
      DO n = -ntor, ntor
!  There needs to be at least two m modes (m=0, m=1) for an equilibrium without
!  zero volume. So there is no need to check the mpol value here.
         IF (n .ge. 0) THEN
            r_bndy(n,1) =      this%t1(1)*rho(n,0)
            z_bndy(n,1) = sign*this%t1(1)*rho(n,0)
         ELSE
            r_bndy(n,1) = sign*this%t1(1)*rho(-n,0)
            z_bndy(n,1) =      this%t1(1)*rho(-n,0)
         END IF
!  Circular crosssection cases lack the m=2 modes. In those case, the rho(m=2)
!  modes should be zero. Note mpol uses the convention from the wout file and is
!  equal to total number of m modes not the maximum m mode.
         IF (mpol - 1 .ge. 2) THEN
            r_bndy(n,1) = r_bndy(n,1) +      this%t2(1)*rho(n,2)
            z_bndy(n,1) = z_bndy(n,1) - sign*this%t2(1)*rho(n,2)
         END IF

         DO m = 2, mpol - 2
            IF (rho(0,m) .ne. 0) THEN
               r_bndy(n,m) = this%t1(m)*rho(n,m-1)                             &
     &                     + this%t2(m)*rho(n,m+1)
               z_bndy(n,m) = sign*(this%t1(m)*rho(n,m-1) -                     &
     &                             this%t2(m)*rho(n,m+1))
            END IF

!  FIXME: This seems to fail.
!            IF (m .lt. mpol - 2) THEN
!               x = rho(0,m+2)
!            ELSE
!               x = 0.0
!            END IF
!            IF (x .eq. 0) THEN
!               tnorm = this%t1(m+2)*this%t1(m+2)/this%t2(m)
!               r_bndy(n,m) = r_bndy(n,m) +      tnorm*rho(n,m+1)
!               z_bndy(n,m) = z_bndy(n,m) - sign*tnorm*rho(n,m+1)
!            END IF
         END DO

         IF (mpol - 1 .ge. 2) THEN
            m = mpol - 1
            r_bndy(n,m) =      this%t1(m)*rho(n,m-1)
            z_bndy(n,m) = sign*this%t1(m)*rho(n,m-1)
         END IF
      END DO

      END SUBROUTINE
      
!-------------------------------------------------------------------------------
!>  @brief Convert from r,z to rho
!>
!>  This routine uses the r,z Fourier components to compute the rho fourier
!>  components. Assuming a fixed center point,
!>
!>  @param[in]  this   A @ref rho_bndy_class instance.
!>  @param[in]  mpol   Number of poloidal modes.
!>  @param[in]  ntor   Number of toroidal modes.
!>  @param[in]  r_bndy R boundary coefficient.
!>  @param[in]  z_bndy Z boundary coefficient.
!>  @param[out] rho    Rho boundary coefficient.
!>  @param[out] r_ctr  Center R bounday cofficient.
!>  @param[out] z_ctr  Center R bounday cofficient.
!>  @param[in]  parity Stellarator symmetric or asymmetric parity.
!-------------------------------------------------------------------------------
      SUBROUTINE rho_bndy_to_rho(this, mpol, ntor, r_bndy, z_bndy, rho,        &
     &                           r_ctr, z_ctr, parity)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (rho_bndy_class), INTENT(in)                        :: this
      INTEGER, INTENT(in)                                       :: ntor
      INTEGER, INTENT(in)                                       :: mpol
      REAL (rprec), DIMENSION(-ntor:ntor,0:mpol-1), INTENT(in)  ::             &
     &   r_bndy
      REAL (rprec), DIMENSION(-ntor:ntor,0:mpol-1), INTENT(in)  ::             &
     &   z_bndy
      REAL (rprec), DIMENSION(-ntor:ntor,0:mpol-1), INTENT(out) :: rho
      REAL (rprec), DIMENSION(0:ntor), INTENT(out)              :: r_ctr
      REAL (rprec), DIMENSION(0:ntor), INTENT(out)              :: z_ctr
      LOGICAL, INTENT(in)                                       ::             &
     &   parity

!  Local Variables
      INTEGER                                                   :: m
      INTEGER                                                   :: n
      REAL (rprec)                                              :: tnorm
      INTEGER                                                   :: sign

!  Start of executable code.
      rho = 0.0
      r_ctr = 0.0
      z_ctr = 0.0

!  Convert center.
      IF (parity .eqv. rho_bndy_sym) THEN
         sign = 1

         r_ctr(0) = r_bndy(0,0)
         z_ctr(0) = 0.0
         IF (mpol - 1 .ge. 2) THEN
            r_ctr(0) = r_ctr(0) - 0.5*(r_bndy(0,2) + z_bndy(0,2))
         END IF
      ELSE
         sign = -1

         r_ctr(0) = 0.0
         z_ctr(0) = z_bndy(0,0)
         IF (mpol - 1 .ge. 2) THEN
            z_ctr(0) = z_ctr(0) - 0.5*(r_bndy(0,2) - z_bndy(0,2))
         END IF
      END IF

      DO n = 1, ntor
         r_ctr(n) = r_bndy(n,0)
         z_ctr(n) = z_bndy(n,0)

         IF (mpol - 1 .ge. 2) THEN
            r_ctr(n) = r_ctr(n) - 0.5*(r_bndy(n,2) +                           &
     &                                 sign*z_bndy(n,2) +                      &
     &                                 sign*r_bndy(-n,2) +                     &
     &                                 z_bndy(-n,2))
            z_ctr(n) = z_ctr(n) - 0.5*(-sign*r_bndy(n,2) -                     &
     &                                 z_bndy(n,2) +                           &
     &                                 r_bndy(-n,2) +                          &
     &                                 sign*z_bndy(-n,2))
         END IF
      END DO

!  Convert boundary.
      DO n = -ntor, ntor
!  m = 0
         IF (n .ge. 0) THEN
            rho(n,0) = 0.5*(r_bndy(n,1) + sign*z_bndy(n,1))/this%t1(1)
         END IF

!  m = 1
         IF (mpol - 1 .ge. 2) THEN
            rho(n,1) = 0.5*(r_bndy(n,2) + sign*z_bndy(n,2))/this%t1(2)
         END IF

!  m > 1
         DO m = 2, mpol - 2
            tnorm = 0.5/(this%t1(m+1)**2 + this%t2(m-1)**2)
            rho(n,m) = tnorm*((r_bndy(n,m+1) +                                 &
     &                         sign*z_bndy(n,m+1))*this%t1(m+1) +              &
     &                        (r_bndy(n,m-1) -                                 &
     &                         sign*z_bndy(n,m-1))*this%t2(m-1))
         END DO

!  Last m mode.
         IF (mpol - 1 .ge. 2) THEN
            m = mpol - 1
            tnorm = 0.5/(this%t1(m+1)**2 + this%t2(m-1)**2)
            rho(n,m) = tnorm*(     r_bndy(n,m-1) -                             &
     &                        sign*z_bndy(n,m-1))*this%t2(m-1)
         END IF
      END DO

      IF (parity .eqv. rho_bndy_asym) THEN
         rho(0,0) = 0.0
      END IF

      END SUBROUTINE

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Rho boundary unit test function.
!>
!>  This runs the associated unit tests and returns the result. Rho conversion
!>  boundaries only work with fouier compressed spectra. So tests will run
!>  rho -> rz -> rho.
!>
!>  @returns True if the tests pass and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION rho_bndy_test()

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                   :: rho_bndy_test

!  Local Variables
      CLASS (rho_bndy_class), POINTER           :: bndy => null()
      INTEGER                                   :: mpol
      INTEGER                                   :: ntor
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: rho_c
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: rho_s
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: r_ctr_c
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: r_ctr_s
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: z_ctr_c
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: z_ctr_s
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: t_rho_c
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: t_rho_s
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: t_r_ctr_c
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: t_r_ctr_s
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: t_z_ctr_c
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: t_z_ctr_s
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: rbc
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: rbs
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: zbc
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: zbs
      INTEGER                                   :: test_num

      INTEGER                                   :: m
      INTEGER                                   :: n

!  Start of executable code.
!-------------------------------------------------------------------------------
!  Test a circular boundary.
!-------------------------------------------------------------------------------
      mpol = 2
      ntor = 0

      bndy => rho_bndy_class(mpol)

      ALLOCATE(rho_c(-ntor:ntor,0:mpol - 1))
      ALLOCATE(r_ctr_c(0:ntor))
      ALLOCATE(z_ctr_s(0:ntor))

      rho_c(0,0) = 0.25
      rho_c(0,1) = 0.0
      r_ctr_c(0) = 1.0
      z_ctr_s(0) = 0.0

      ALLOCATE(rbc(-ntor:ntor,0:mpol - 1))
      ALLOCATE(zbs(-ntor:ntor,0:mpol - 1))

      CALL bndy%to_rz(mpol, ntor, rho_c, r_ctr_c, z_ctr_s, rbc, zbs,           &
     &                rho_bndy_sym)

      rho_bndy_test = check(1.0_rprec, rbc(0,0), 1,                            &
     &                      'circular_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, rbc(0,1), 2,                           &
     &                      'circular_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbs(0,0), 3,                            &
     &                      'circular_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, zbs(0,1), 4,                           &
     &                      'circular_boundary')
      IF (.not.rho_bndy_test) RETURN

!  Make some duplate rho buffers to check that transformed match the original.
      ALLOCATE(t_rho_c(-ntor:ntor,0:mpol - 1))
      ALLOCATE(t_r_ctr_c(0:ntor))
      ALLOCATE(t_z_ctr_s(0:ntor))

      CALL bndy%to_rho(mpol, ntor, rbc, zbs, t_rho_c, t_r_ctr_c,               &
     &                 t_z_ctr_s, rho_bndy_sym)

      test_num = 4

      rho_bndy_test = check(rho_c, t_rho_c, test_num,                          &
     &                      'circular_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(r_ctr_c, t_r_ctr_c, test_num,                      &
     &                      'circular_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(z_ctr_s, t_z_ctr_s, test_num,                      &
     &                      'circular_boundary')
      IF (.not.rho_bndy_test) RETURN

      DEALLOCATE(t_rho_c)
      DEALLOCATE(t_r_ctr_c)
      DEALLOCATE(t_z_ctr_s)

      DEALLOCATE(rbc)
      DEALLOCATE(zbs)

      DEALLOCATE(rho_c)
      DEALLOCATE(r_ctr_c)
      DEALLOCATE(z_ctr_s)

!-------------------------------------------------------------------------------
!  Test circular boundary with an asymmetric shift.
!-------------------------------------------------------------------------------
      ALLOCATE(rho_s(-ntor:ntor,0:mpol - 1))
      ALLOCATE(r_ctr_s(0:ntor))
      ALLOCATE(z_ctr_c(0:ntor))

      rho_s(0,0) = 0.0
      rho_s(0,1) = 0.0
      r_ctr_s(0) = 0.0
      z_ctr_c(0) = 1.0

      ALLOCATE(rbs(-ntor:ntor,0:mpol - 1))
      ALLOCATE(zbc(-ntor:ntor,0:mpol - 1))

      CALL bndy%to_rz(mpol, ntor, rho_s, r_ctr_s, z_ctr_c, rbs, zbc,           &
     &                rho_bndy_asym)

      rho_bndy_test = check(0.0_rprec, rbs(0,0), 1,                            &
     &                      'circular_asym_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, rbs(0,1), 2,                            &
     &                      'circular_asym_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(1.0_rprec, zbc(0,0), 3,                            &
     &                      'circular_asym_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbc(0,1), 4,                            &
     &                      'circular_asym_boundary')
      IF (.not.rho_bndy_test) RETURN

!  Make some duplate rho buffers to check that transformed match the original.
      ALLOCATE(t_rho_s(-ntor:ntor,0:mpol - 1))
      ALLOCATE(t_r_ctr_s(0:ntor))
      ALLOCATE(t_z_ctr_c(0:ntor))

      CALL bndy%to_rho(mpol, ntor, rbs, zbc, t_rho_s, t_r_ctr_s,               &
     &                 t_z_ctr_c, rho_bndy_asym)

      test_num = 4

      rho_bndy_test = check(rho_s, t_rho_s, test_num,                          &
     &                      'circular_asym_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(r_ctr_s, t_r_ctr_s, test_num,                      &
     &                      'circular_asym_boundary')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(z_ctr_c, t_z_ctr_c, test_num,                      &
     &                      'circular_asym_boundary')
      IF (.not.rho_bndy_test) RETURN

      DEALLOCATE(t_rho_s)
      DEALLOCATE(t_r_ctr_s)
      DEALLOCATE(t_z_ctr_c)

      DEALLOCATE(rbs)
      DEALLOCATE(zbc)

      DEALLOCATE(rho_s)
      DEALLOCATE(r_ctr_s)
      DEALLOCATE(z_ctr_c)

!-------------------------------------------------------------------------------
!  Test a non-axisymmetric circular boundary.
!-------------------------------------------------------------------------------
      mpol = 2
      ntor = 1

      ALLOCATE(rho_c(-ntor:ntor,0:mpol - 1))
      ALLOCATE(r_ctr_c(0:ntor))
      ALLOCATE(z_ctr_s(0:ntor))

      rho_c(0,0) = 0.25
      rho_c(0,1) = 0.0
      rho_c(1,0) = 0.0
      rho_c(-1,0) = 0.0
      rho_c(1,1) = 0.0
      rho_c(-1,1) = 0.0
      r_ctr_c(0) = 1.0
      r_ctr_c(1) = 0.25
      z_ctr_s(0) = 0.0
      z_ctr_s(1) = 0.25

      ALLOCATE(rbc(-ntor:ntor,0:mpol - 1))
      ALLOCATE(zbs(-ntor:ntor,0:mpol - 1))

      CALL bndy%to_rz(mpol, ntor, rho_c, r_ctr_c, z_ctr_s, rbc, zbs,           &
     &                rho_bndy_sym)

      rho_bndy_test = check(1.0_rprec, rbc(0,0), 1,                            &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, rbc(0,1), 2,                           &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, rbc(1,0), 3,                           &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, rbc(-1,0), 4,                           &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, rbc(1,1), 5,                            &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, rbc(-1,1), 6,                           &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbs(0,0), 7,                            &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, zbs(0,1), 8,                           &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, zbs(1,0), 9,                           &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbs(-1,0), 10,                          &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbs(1,1), 11,                           &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbs(-1,1), 12,                          &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

!  Make some duplate rho buffers to check that transformed match the original.
      ALLOCATE(t_rho_c(-ntor:ntor,0:mpol - 1))
      ALLOCATE(t_r_ctr_c(0:ntor))
      ALLOCATE(t_z_ctr_s(0:ntor))

      CALL bndy%to_rho(mpol, ntor, rbc, zbs, t_rho_c, t_r_ctr_c,               &
     &                 t_z_ctr_s, rho_bndy_sym)

      test_num = 12

      rho_bndy_test = check(rho_c, t_rho_c, test_num,                          &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(r_ctr_c, t_r_ctr_c, test_num,                      &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(z_ctr_s, t_z_ctr_s, test_num,                      &
     &                      'circular_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      DEALLOCATE(t_rho_c)
      DEALLOCATE(t_r_ctr_c)
      DEALLOCATE(t_z_ctr_s)

      DEALLOCATE(rbc)
      DEALLOCATE(zbs)

      DEALLOCATE(rho_c)
      DEALLOCATE(r_ctr_c)
      DEALLOCATE(z_ctr_s)

!-------------------------------------------------------------------------------
!  Test a non-axisymmetric circular boundary with an asymmetric shift.
!-------------------------------------------------------------------------------
      ALLOCATE(rho_s(-ntor:ntor,0:mpol - 1))
      ALLOCATE(r_ctr_s(0:ntor))
      ALLOCATE(z_ctr_c(0:ntor))

      rho_s(0,0) = 0.0
      rho_s(0,1) = 0.0
      rho_s(1,0) = 0.0
      rho_s(-1,0) = 0.0
      rho_s(1,1) = 0.0
      rho_s(-1,1) = 0.0
      r_ctr_s(0) = 0.0
      r_ctr_s(1) = 0.25
      z_ctr_c(0) = 0.25
      z_ctr_c(1) = 0.25

      ALLOCATE(rbs(-ntor:ntor,0:mpol - 1))
      ALLOCATE(zbc(-ntor:ntor,0:mpol - 1))

      CALL bndy%to_rz(mpol, ntor, rho_s, r_ctr_s, z_ctr_c, rbs, zbc,           &
     &                rho_bndy_asym)

      rho_bndy_test = check(0.0_rprec, rbs(0,0), 1,                            &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, rbs(0,1), 2,                           &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, rbs(1,0), 3,                           &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, rbs(-1,0), 4,                           &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, rbs(1,1), 5,                            &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, rbs(-1,1), 6,                           &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, zbc(0,0), 7,                           &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbc(0,1), 8,                            &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.25_rprec, zbc(1,0), 9,                           &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbc(-1,0), 10,                          &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbc(1,1), 11,                           &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(0.0_rprec, zbc(-1,1), 12,                          &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

!  Make some duplate rho buffers to check that transformed match the original.
      ALLOCATE(t_rho_s(-ntor:ntor,0:mpol - 1))
      ALLOCATE(t_r_ctr_s(0:ntor))
      ALLOCATE(t_z_ctr_c(0:ntor))

      CALL bndy%to_rho(mpol, ntor, rbs, zbc, t_rho_s, t_r_ctr_s,               &
     &                 t_z_ctr_c, rho_bndy_asym)

      test_num = 12

      rho_bndy_test = check(rho_s, t_rho_s, test_num,                          &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(r_ctr_s, t_r_ctr_s, test_num,                      &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(z_ctr_c, t_z_ctr_c, test_num,                      &
     &                      'circular_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      DEALLOCATE(t_rho_s)
      DEALLOCATE(t_r_ctr_s)
      DEALLOCATE(t_z_ctr_c)

      DEALLOCATE(rbs)
      DEALLOCATE(zbc)

      DEALLOCATE(rho_s)
      DEALLOCATE(r_ctr_s)
      DEALLOCATE(z_ctr_c)

      DEALLOCATE(bndy)

!-------------------------------------------------------------------------------
!  Test a non-axisymmetric non-circular boundary.
!-------------------------------------------------------------------------------
      mpol = 3
      ntor = 2

      bndy => rho_bndy_class(mpol)

      ALLOCATE(rho_c(-ntor:ntor,0:mpol - 1))
      ALLOCATE(r_ctr_c(0:ntor))
      ALLOCATE(z_ctr_s(0:ntor))

      rho_c(2,0) = 0.02
      rho_c(1,0) = 0.2
      rho_c(0,0) = 0.25
      rho_c(-1,0) = 0.0
      rho_c(-2,0) = 0.0
      rho_c(2,1) = -0.05
      rho_c(1,1) = -0.15
      rho_c(0,1) = 0.1
      rho_c(-1,1) = 0.3
      rho_c(-2,1) = 0.03
      rho_c(2,2) = 0.0
      rho_c(1,2) = 0.0
      rho_c(0,2) = 0.0
      rho_c(-1,2) = 0.0
      rho_c(-2,2) = 0.0
      r_ctr_c(0) = 1.0
      r_ctr_c(1) = 0.25
      r_ctr_c(2) = 0.05
      z_ctr_s(0) = 0.0
      z_ctr_s(1) = 0.25
      z_ctr_s(2) = 0.01

      ALLOCATE(rbc(-ntor:ntor,0:mpol - 1))
      ALLOCATE(zbs(-ntor:ntor,0:mpol - 1))

      CALL bndy%to_rz(mpol, ntor, rho_c, r_ctr_c, z_ctr_s, rbc, zbs,           &
     &                rho_bndy_sym)

!  Make some duplate rho buffers to check that transformed match the original.
      ALLOCATE(t_rho_c(-ntor:ntor,0:mpol - 1))
      ALLOCATE(t_r_ctr_c(0:ntor))
      ALLOCATE(t_z_ctr_s(0:ntor))

      CALL bndy%to_rho(mpol, ntor, rbc, zbs, t_rho_c, t_r_ctr_c,               &
     &                 t_z_ctr_s, rho_bndy_sym)

      test_num = 0

      rho_bndy_test = check(rho_c, t_rho_c, test_num,                          &
     &                      'm3_n2_sym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(r_ctr_c, t_r_ctr_c, test_num,                      &
     &                      'm3_n2_sym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(z_ctr_s, t_z_ctr_s, test_num,                      &
     &                      'm3_n2_sym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      DEALLOCATE(t_rho_c)
      DEALLOCATE(t_r_ctr_c)
      DEALLOCATE(t_z_ctr_s)

      DEALLOCATE(rbc)
      DEALLOCATE(zbs)

      DEALLOCATE(rho_c)
      DEALLOCATE(r_ctr_c)
      DEALLOCATE(z_ctr_s)

!-------------------------------------------------------------------------------
!  Test a non-axisymmetric non-circular boundary.
!-------------------------------------------------------------------------------
      ALLOCATE(rho_s(-ntor:ntor,0:mpol - 1))
      ALLOCATE(r_ctr_s(0:ntor))
      ALLOCATE(z_ctr_c(0:ntor))

      rho_s(2,0) = 0.02
      rho_s(1,0) = 0.2
      rho_s(0,0) = 0.0
      rho_s(-1,0) = 0.0
      rho_s(-2,0) = 0.0
      rho_s(2,1) = -0.05
      rho_s(1,1) = -0.15
      rho_s(0,1) = 0.1
      rho_s(-1,1) = 0.3
      rho_s(-2,1) = 0.03
      rho_s(2,2) = 0.0
      rho_s(1,2) = 0.0
      rho_s(0,2) = 0.0
      rho_s(-1,2) = 0.0
      rho_s(-2,2) = 0.0
      r_ctr_s(0) = 0.0
      r_ctr_s(1) = 0.25
      r_ctr_s(2) = 0.05
      z_ctr_c(0) = 1.0
      z_ctr_c(1) = 0.25
      z_ctr_c(2) = 0.01

      ALLOCATE(rbs(-ntor:ntor,0:mpol - 1))
      ALLOCATE(zbc(-ntor:ntor,0:mpol - 1))

      CALL bndy%to_rz(mpol, ntor, rho_s, r_ctr_s, z_ctr_c, rbs, zbc,           &
     &                rho_bndy_asym)

!  Make some duplate rho buffers to check that transformed match the original.
      ALLOCATE(t_rho_s(-ntor:ntor,0:mpol - 1))
      ALLOCATE(t_r_ctr_s(0:ntor))
      ALLOCATE(t_z_ctr_c(0:ntor))

      CALL bndy%to_rho(mpol, ntor, rbs, zbc, t_rho_s, t_r_ctr_s,               &
     &                 t_z_ctr_c, rho_bndy_asym)

      test_num = 0

      rho_bndy_test = check(rho_s, t_rho_s, test_num,                          &
     &                      'm3_n2_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(r_ctr_s, t_r_ctr_s, test_num,                      &
     &                      'm3_n2_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(z_ctr_c, t_z_ctr_c, test_num,                      &
     &                      'm3_n2_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      DEALLOCATE(t_rho_s)
      DEALLOCATE(t_r_ctr_s)
      DEALLOCATE(t_z_ctr_c)

      DEALLOCATE(rbs)
      DEALLOCATE(zbc)

      DEALLOCATE(rho_s)
      DEALLOCATE(r_ctr_s)
      DEALLOCATE(z_ctr_c)

      DEALLOCATE(bndy)

!-------------------------------------------------------------------------------
!  Test a non-axisymmetric non-circular boundary.
!-------------------------------------------------------------------------------
      mpol = 5
      ntor = 3

      bndy => rho_bndy_class(mpol)

      ALLOCATE(rho_c(-ntor:ntor,0:mpol - 1))
      ALLOCATE(r_ctr_c(0:ntor))
      ALLOCATE(z_ctr_s(0:ntor))

      rho_c(3,0) = -0.001
      rho_c(2,0) = 0.02
      rho_c(1,0) = 0.2
      rho_c(0,0) = 0.25
      rho_c(-1,0) = 0.0
      rho_c(-2,0) = 0.0
      rho_c(-3,0) = 0.0
      rho_c(3,1) = -0.002
      rho_c(2,1) = -0.05
      rho_c(1,1) = -0.15
      rho_c(0,1) = 0.1
      rho_c(-1,1) = 0.3
      rho_c(-2,1) = 0.03
      rho_c(-3,1) = -0.003
      rho_c(3,2) = -0.0002
      rho_c(2,2) = -0.005
      rho_c(1,2) = -0.015
      rho_c(0,2) = 0.01
      rho_c(-1,2) = 0.03
      rho_c(-2,2) = 0.003
      rho_c(-3,2) = -0.0003
      rho_c(3,3) = -0.00002
      rho_c(2,3) = -0.0005
      rho_c(1,3) = -0.0015
      rho_c(0,3) = 0.001
      rho_c(-1,3) = 0.003
      rho_c(-2,3) = 0.0003
      rho_c(-3,3) = -0.00003
      rho_c(3,4) = 0.0
      rho_c(2,4) = 0.0
      rho_c(1,4) = 0.0
      rho_c(0,4) = 0.0
      rho_c(-1,4) = 0.0
      rho_c(-2,4) = 0.0
      rho_c(-3,4) = 0.0
      r_ctr_c(0) = 1.0
      r_ctr_c(1) = 0.25
      r_ctr_c(2) = 0.05
      r_ctr_c(3) = -0.01
      z_ctr_s(0) = 0.0
      z_ctr_s(1) = 0.25
      z_ctr_s(2) = 0.01
      z_ctr_s(3) = 0.005

      ALLOCATE(rbc(-ntor:ntor,0:mpol - 1))
      ALLOCATE(zbs(-ntor:ntor,0:mpol - 1))

      CALL bndy%to_rz(mpol, ntor, rho_c, r_ctr_c, z_ctr_s, rbc, zbs,           &
     &                rho_bndy_sym)

!  Make some duplate rho buffers to check that transformed match the original.
      ALLOCATE(t_rho_c(-ntor:ntor,0:mpol - 1))
      ALLOCATE(t_r_ctr_c(0:ntor))
      ALLOCATE(t_z_ctr_s(0:ntor))

      CALL bndy%to_rho(mpol, ntor, rbc, zbs, t_rho_c, t_r_ctr_c,               &
     &                 t_z_ctr_s, rho_bndy_sym)

      test_num = 0

      rho_bndy_test = check(rho_c, t_rho_c, test_num,                          &
     &                      'm5_n3_sym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(r_ctr_c, t_r_ctr_c, test_num,                      &
     &                      'm5_n3_sym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(z_ctr_s, t_z_ctr_s, test_num,                      &
     &                      'm5_n3_sym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      DEALLOCATE(t_rho_c)
      DEALLOCATE(t_r_ctr_c)
      DEALLOCATE(t_z_ctr_s)

      DEALLOCATE(rbc)
      DEALLOCATE(zbs)

      DEALLOCATE(rho_c)
      DEALLOCATE(r_ctr_c)
      DEALLOCATE(z_ctr_s)

!-------------------------------------------------------------------------------
!  Test a non-axisymmetric non-circular boundary with an asymmetric shift.
!-------------------------------------------------------------------------------
      ALLOCATE(rho_s(-ntor:ntor,0:mpol - 1))
      ALLOCATE(r_ctr_s(0:ntor))
      ALLOCATE(z_ctr_c(0:ntor))

      rho_s(3,0) = -0.001
      rho_s(2,0) = 0.02
      rho_s(1,0) = 0.2
      rho_s(0,0) = 0.0
      rho_s(-1,0) = 0.0
      rho_s(-2,0) = 0.0
      rho_s(-3,0) = 0.0
      rho_s(3,1) = -0.002
      rho_s(2,1) = -0.05
      rho_s(1,1) = -0.15
      rho_s(0,1) = 0.1
      rho_s(-1,1) = 0.3
      rho_s(-2,1) = 0.03
      rho_s(-3,1) = -0.003
      rho_s(3,2) = -0.0002
      rho_s(2,2) = -0.005
      rho_s(1,2) = -0.015
      rho_s(0,2) = 0.01
      rho_s(-1,2) = 0.03
      rho_s(-2,2) = 0.003
      rho_s(-3,2) = -0.0003
      rho_s(3,3) = -0.00002
      rho_s(2,3) = -0.0005
      rho_s(1,3) = -0.0015
      rho_s(0,3) = 0.001
      rho_s(-1,3) = 0.003
      rho_s(-2,3) = 0.0003
      rho_s(-3,3) = -0.00003
      rho_s(3,4) = 0.0
      rho_s(2,4) = 0.0
      rho_s(1,4) = 0.0
      rho_s(0,4) = 0.0
      rho_s(-1,4) = 0.0
      rho_s(-2,4) = 0.0
      rho_s(-3,4) = 0.0
      r_ctr_s(0) = 0.0
      r_ctr_s(1) = 0.25
      r_ctr_s(2) = 0.05
      r_ctr_s(3) = -0.01
      z_ctr_c(0) = 1.0
      z_ctr_c(1) = 0.25
      z_ctr_c(2) = 0.01
      z_ctr_c(3) = 0.005

      ALLOCATE(rbs(-ntor:ntor,0:mpol - 1))
      ALLOCATE(zbc(-ntor:ntor,0:mpol - 1))

      CALL bndy%to_rz(mpol, ntor, rho_s, r_ctr_s, z_ctr_c, rbs, zbc,           &
     &                rho_bndy_asym)

!  Make some duplate rho buffers to check that transformed match the original.
      ALLOCATE(t_rho_s(-ntor:ntor,0:mpol - 1))
      ALLOCATE(t_r_ctr_s(0:ntor))
      ALLOCATE(t_z_ctr_c(0:ntor))

      CALL bndy%to_rho(mpol, ntor, rbs, zbc, t_rho_s, t_r_ctr_s,               &
     &                 t_z_ctr_c, rho_bndy_asym)

      test_num = 0

      rho_bndy_test = check(rho_s, t_rho_s, test_num,                          &
     &                      'm5_n3_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(r_ctr_s, t_r_ctr_s, test_num,                      &
     &                      'm5_n3_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      rho_bndy_test = check(z_ctr_c, t_z_ctr_c, test_num,                      &
     &                      'm5_n3_asym_boundary_3D')
      IF (.not.rho_bndy_test) RETURN

      DEALLOCATE(t_rho_s)
      DEALLOCATE(t_r_ctr_s)
      DEALLOCATE(t_z_ctr_c)

      DEALLOCATE(rbs)
      DEALLOCATE(zbc)

      DEALLOCATE(rho_s)
      DEALLOCATE(r_ctr_s)
      DEALLOCATE(z_ctr_c)

      DEALLOCATE(bndy)

      END FUNCTION

!*******************************************************************************
!  CHECK FUNCTIONS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Check two real values to match within a tolarance.
!>
!>  @param[in] expected Expected value for the test.
!>  @param[in] received Recieve value for the test.
!>  @param[in] test_num Test number.
!>  @param[in] name     Name of the test.
!-------------------------------------------------------------------------------
      FUNCTION check_real(expected, received, test_num, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_real
      REAL (rprec), INTENT(in)      :: expected
      REAL (rprec), INTENT(in)      :: received
      INTEGER, INTENT(in)           :: test_num
      CHARACTER (len=*), INTENT(in) :: name

!  Local Parameters
      REAL (rprec), PARAMETER       :: eps = 1.0E-8

!  Start of executable code.
      check_real = (ABS(expected - received) .lt. eps)
      IF (.not.check_real) THEN
         WRITE (*,1000) TRIM(name), test_num, expected, received
      END IF

1000  FORMAT('rho_bndy_input.f: ',a,' test ',i3,' failed.',/,
     &       'Expected ',e12.3,' Recieved ',e12.3)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check 2D Fourier harmonics for all modes.
!>
!>  @param[in]    expected Expected value for the test.
!>  @param[in]    received Recieve value for the test.
!>  @param[inout] test_num Test number.
!>  @param[in]    name     Name of the test.
!>  @param[in]    mpol     Number of poloidal modes.
!>  @param[in]    ntor     Number of toroidal modes.
!-------------------------------------------------------------------------------
      FUNCTION check_fourier_1d(expected, received, test_num, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                :: check_fourier_1d
      REAL (rprec), DIMENSION(:), INTENT(in) :: expected
      REAL (rprec), DIMENSION(:), INTENT(in) :: received
      INTEGER, INTENT(inout)                 :: test_num
      CHARACTER (len=*), INTENT(in)          :: name

!  Local Variables
      INTEGER                                :: n
      CHARACTER (len=200)                    :: temp_name

!  Start of executable code.
      DO n = LBOUND(expected, 1), UBOUND(expected, 1)
         test_num = test_num + 1
         WRITE (temp_name, 1000) TRIM(name), n - LBOUND(expected, 1)
         check_fourier_1d = check(expected(n), received(n),                    &
     &                            test_num, temp_name)
         IF (.not.check_fourier_1d) RETURN
      END DO

1000  FORMAT(a,' n = ',i3)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check 2D Fourier harmonics for all modes.
!>
!>  @param[in]    expected Expected value for the test.
!>  @param[in]    received Recieve value for the test.
!>  @param[inout] test_num Test number.
!>  @param[in]    name     Name of the test.
!>  @param[in]    mpol     Number of poloidal modes.
!>  @param[in]    ntor     Number of toroidal modes.
!-------------------------------------------------------------------------------
      FUNCTION check_fourier_2d(expected, received, test_num, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                  :: check_fourier_2d
      REAL (rprec), DIMENSION(:,:), INTENT(in) :: expected
      REAL (rprec), DIMENSION(:,:), INTENT(in) :: received
      INTEGER, INTENT(inout)                   :: test_num
      CHARACTER (len=*), INTENT(in)            :: name

!  Local Variables
      INTEGER                                  :: m
      INTEGER                                  :: n
      CHARACTER (len=200)                      :: temp_name

!  Start of executable code.
      DO m = LBOUND(expected, 2), UBOUND(expected, 2)
         DO n = LBOUND(expected, 1), UBOUND(expected, 1)
            test_num = test_num + 1
            WRITE (temp_name, 1000) TRIM(name),                                &
     &                              m - LBOUND(expected, 2),                   &
     &                              n - LBOUND(expected, 1)
            check_fourier_2d = check(expected(n,m), received(n,m),             &
     &                               test_num, temp_name)
            IF (.not.check_fourier_2d) RETURN
         END DO
      END DO

1000  FORMAT(a,' m = ',i3,' n = ',i3)

      END FUNCTION

      END MODULE
