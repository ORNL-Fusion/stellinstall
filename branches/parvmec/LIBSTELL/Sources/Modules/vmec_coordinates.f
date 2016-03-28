!*******************************************************************************
!>  @file equilibrium.f
!>  @brief Contains module @ref vmec_coordinates
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref vmec_coordinates_class. This class
!>  holds information for working in the VMEC coordinate system.
!*******************************************************************************

      MODULE vmec_coordinates
      USE stel_kinds
      USE stel_constants

      IMPLICIT NONE

!*******************************************************************************
!  vmec_coordinates module parameters
!*******************************************************************************
!>  S position outside last closed flux surface.
      INTEGER, PARAMETER :: vmec_coordinates_bad_s_flag = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) vmec_coordinates base class
!  1) vmec_coordinates base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a point in the vmec coordinate system.
!>  @ref vmec_coordinates_class
!-------------------------------------------------------------------------------
      TYPE vmec_coordinates_class
!>  Real space position in cylindical coordinates.
         REAL (rprec), DIMENSION(3)          :: cyl
!>  Flux space positon.
         REAL (rprec), DIMENSION(3)          :: flux

!>  esub
         REAL (rprec), DIMENSION(3,3)        :: esub
!>  esup
         REAL (rprec), DIMENSION(3,3)        :: esup

!>  bsup
         REAL (rprec), DIMENSION(3)          :: bsup
!>  bsub
         REAL (rprec), DIMENSION(3)          :: bsub

!>  asup
         REAL (rprec), DIMENSION(3)          :: asup
!>  asub
         REAL (rprec), DIMENSION(3)          :: asub

!>  Status of coordinates.
         INTEGER                             :: status
         REAL (rprec)                        :: tol
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      PRIVATE cross_product

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref vmec_coordinates_class from a real position.
!>
!>  Allocates memory and initializes a @ref vmec_coordinates_class object.
!>
!>  @param[in] cyl Real space position in cylindical coordinates.
!>  @returns A pointer to a constructed @ref vmec_coordinates_class object.
!-------------------------------------------------------------------------------
      FUNCTION vmec_coordinates_real(cyl)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_coordinates_class), POINTER :: vmec_coordinates_real
      REAL (rprec), DIMENSION(3)             :: cyl

!  Start of executable code
      ALLOCATE (vmec_coordinates_real)

      CALL vmec_coordinates_set_flux(vmec_coordinates_real, cyl)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref vmec_coordinates_class from a flux position.
!>
!>  Allocates memory and initializes a @ref vmec_coordinates_class object.
!>
!>  @param[in] flux Flux space position.
!>  @returns A pointer to a constructed @ref vmec_coordinates_class object.
!-------------------------------------------------------------------------------
      FUNCTION vmec_coordinates_flux(flux)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_coordinates_class), POINTER :: vmec_coordinates_flux
      REAL (rprec), DIMENSION(3)             :: flux

!  Start of executable code
      ALLOCATE (vmec_coordinates_flux)

      CALL vmec_coordinates_set_flux(vmec_coordinates_flux, flux)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE vmec_coordinates_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_coordinates_class), POINTER :: this

!  Start of executable code
      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sets the real space position.
!>
!>  Sets the real space position and updates the current vmec coordinates.
!>
!>  @param[inout] this A @ref vmec_coordinates_class instance.
!>  @param[in]    cyl  Real space position in cylindical coordinates.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_coordinates_set_real(this, cyl)
      USE read_wout_mod
      USE cyl_flux

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_coordinates_class), INTENT(inout) :: this
      REAL (rprec), DIMENSION(3), INTENT(in)       :: cyl

!  local variables
      INTEGER                                      :: status
      INTEGER                                      :: nfe
      REAL (rprec)                                 :: fmin
      REAL (rprec), DIMENSION(3)                   :: cyl_local

!  Start of executable code
      this%cyl = cyl
      this%status = 0

      cyl_local = cyl
      cyl_local(2) = nfp*cyl(2)

      CALL cyl2flx(rzl_local, cyl_local, this%flux, ns, ntor, mpol,            &
     &             ntmax, lthreed, lasym, status, nfe, this%tol)
      this%flux(3) = this%flux(3)/nfp

      SELECT CASE (status)

         CASE DEFAULT
            IF (this%flux(1) .gt. 1.0) THEN
               this%status = IBSET(this%status,                                &
     &                             vmec_coordinates_bad_s_flag)
            END IF

         CASE (-3)
            this%status = IBSET(this%status,                                   &
     &                          vmec_coordinates_bad_s_flag)

         CASE (-1)
            IF (this%flux(1) .gt. 1.0) THEN
               this%status = IBSET(this%status,                                &
     &                             vmec_coordinates_bad_s_flag)
            END IF

      END SELECT

      IF (this%status .eq. 0) THEN
         CALL vmec_coordinates_set_flux(this, this%flux)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Sets the flux space position.
!>
!>  Sets the flux space position and updates the current vmec coordinates.
!>
!>  @param[inout] this A @ref vmec_coordinates_class instance.
!>  @param[in]    flux Flux space position in cylindical coordinates.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_coordinates_set_flux(this, flux)
      USE read_wout_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_coordinates_class), INTENT(inout) :: this
      REAL (rprec), DIMENSION(3)                   :: flux

!  local variables
      REAL (rprec)                                 :: ds
      REAL (rprec)                                 :: phi_temp
      REAL (rprec)                                 :: chi_temp
      REAL (rprec)                                 :: jacobian

      INTEGER                                      :: i_low_full
      INTEGER                                      :: i_high_full
      INTEGER                                      :: i_low_half
      INTEGER                                      :: i_high_half
      REAL (rprec)                                 :: w_low_full
      REAL (rprec)                                 :: w_high_full
      REAL (rprec)                                 :: w_low_half
      REAL (rprec)                                 :: w_high_half

      REAL (rprec), DIMENSION(:), ALLOCATABLE      :: cosz
      REAL (rprec), DIMENSION(:), ALLOCATABLE      :: sinz

      REAL (rprec), DIMENSION(3)                   :: a_cyl

!  Start of executable code
      this%flux = flux

      this%status = 0
      IF (flux(1) .gt. 1.0) THEN
         this%status = IBSET(this%status, vmec_coordinates_bad_s_flag)
         RETURN
      END IF

      ds = 1.0/(ns - 1)

!  Find the interpolation weights and indicies for full grid.
      i_low_full = 1 + INT(this%flux(1)*(ns - 1))
      IF (i_low_full .eq. ns) THEN
         i_low_full = ns - 1
      END IF
      i_high_full = i_low_full + 1

      w_low_full = i_low_full - this%flux(1)*(ns - 1)
      w_high_full = 1.0 - w_low_full

!  Find the interpolation weights and indicies for half grid.
      i_low_half = INT(1.5 + this%flux(1)*(ns - 1))
      IF (i_low_half .eq. 1) THEN
         i_low_half = 2
      ELSE IF (i_low_half .eq. ns) THEN
         i_low_half = ns - 1
      END IF
      i_high_half = i_low_half + 1

      w_low_half = i_low_half - 0.5 - this%flux(1)*(ns - 1.0)
      w_high_half = 1.0 - w_low_half

      ALLOCATE(cosz(mnmax))
      ALLOCATE(sinz(mnmax))

      cosz = COS(xm*this%flux(2) - xn*this%flux(3))
      sinz = SIN(xm*this%flux(2) - xn*this%flux(3))

      this%cyl(2) = this%flux(3)
      this%esub = 0.0
      this%cyl(1) =                                                            &
     &   SUM((w_low_full*rmnc(:,i_low_full) +                                  &
     &        w_high_full*rmnc(:,i_high_full))*cosz)

       this%cyl(3) =                                                           &
     &    SUM((w_low_full*zmns(:,i_low_full) +                                 &
     &         w_high_full*zmns(:,i_high_full))*sinz)

      this%esub(1,1) = SUM((rmnc(:,i_high_full) - rmnc(:,i_low_full))          &
     &                     /ds*cosz)
      this%esub(3,1) = SUM((zmns(:,i_high_full) - zmns(:,i_low_full))          &
     &                     /ds*sinz)

      this%esub(1,2) = SUM(-xm*(w_low_full*rmnc(:,i_low_full) +                &
     &                          w_high_full*rmnc(:,i_high_full))*sinz)
      this%esub(3,2) = SUM(xm*(w_low_full*zmns(:,i_low_full) +                 &
     &                         w_high_full*zmns(:,i_high_full))*cosz)

      this%esub(1,3) = SUM(xn*(w_low_full*rmnc(:,i_low_full) +                 &
     &                         w_high_full*rmnc(:,i_high_full))*sinz)
      this%esub(3,3) = SUM(-xn*(w_low_full*zmns(:,i_low_full) +                &
     &                          w_high_full*zmns(:,i_high_full))*cosz)

      this%asub(1) = SUM((lmns(:,i_high_full) - lmns(:,i_low_full))            &
     &                   /ds*sinz)
      this%asub(2) = SUM(xm*(w_low_full*lmns(:,i_low_full) +                   &
     &                       w_high_full*lmns(:,i_high_full))*cosz)
      this%asub(3) = SUM(-xn*(w_low_full*lmns(:,i_low_full) +                  &
     &                        w_high_full*lmns(:,i_high_full))*cosz)

      IF (lasym) THEN
         this%cyl(1) = this%cyl(1) +                                           &
     &      SUM((w_low_full*rmns(:,i_low_full) +                               &
     &           w_high_full*rmns(:,i_high_full))*sinz)

         this%cyl(3) = this%cyl(3) +                                           &
     &      SUM((w_low_full*zmnc(:,i_low_full) +                               &
     &           w_high_full*zmnc(:,i_high_full))*cosz)

         this%esub(1,1) = this%esub(1,1) +                                     &
     &      SUM((rmns(:,i_high_full) - rmns(:,i_low_full))/ds*sinz)
         this%esub(3,1) = this%esub(3,1) +                                     &
     &      SUM((zmnc(:,i_high_full) - zmnc(:,i_low_full))/ds*cosz)

         this%esub(1,2) = this%esub(1,2) +                                     &
     &      SUM(xm*(w_low_full*rmns(:,i_low_full) +                            &
     &              w_high_full*rmns(:,i_high_full))*cosz)
         this%esub(3,2) = this%esub(3,2) +                                     &
     &      SUM(-xm*(w_low_full*zmnc(:,i_low_full) +                           &
     &               w_high_full*zmnc(:,i_high_full))*sinz)

         this%esub(1,3) = this%esub(1,3) +                                     &
     &      SUM(-xn*(w_low_full*rmns(:,i_low_full) +                           &
     &               w_high_full*rmns(:,i_high_full))*cosz)
         this%esub(3,3) = this%esub(3,3) +                                     &
     &      SUM(xn*(w_low_full*zmns(:,i_low_full) +                            &
     &              w_high_full*zmns(:,i_high_full))*sinz)

         this%asub(1) = this%asub(1) +                                         &
     &      SUM((lmnc(:,i_high_full) - lmnc(:,i_low_full))/ds*cosz)
         this%asub(2) = this%asub(2) +                                         &
     &      SUM(-xm*(w_low_full*lmnc(:,i_low_full) +                           &
     &               w_high_full*lmnc(:,i_high_full))*sinz)
         this%asub(3) = this%asub(3) +                                         &
     &      SUM(xn*(w_low_full*lmnc(:,i_low_full) +                            &
     &              w_high_full*lmnc(:,i_high_full))*sinz)
      END IF
      this%esub(2,3) = this%cyl(1)

      phi_temp = w_low_full*phi(i_low_full)                                    &
     &         + w_high_full*phi(i_high_full)
      chi_temp = w_low_full*chi(i_low_full)                                    &
     &         + w_high_full*chi(i_high_full)

      this%asub(1) = phi_temp*this%asub(1)
      this%asub(2) = phi_temp*(1.0 + this%asub(2))
      this%asub(3) = -chi_temp + phi_temp*this%asub(3)

      jacobian = DOT_PRODUCT(this%esub(:,1),                                   &
     &                       cross_product(this%esub(:,2),                     &
     &                                     this%esub(:,3)))

      this%esup(:,1) = cross_product(this%esub(:,2), this%esub(:,3))           &
     &               / jacobian
      this%esup(:,2) = cross_product(this%esub(:,3), this%esub(:,1))           &
     &               / jacobian
      this%esup(:,3) = cross_product(this%esub(:,1), this%esub(:,2))           &
     &               / jacobian

      a_cyl = vmec_coordinates_get_A_cyl(this)
      this%asup(1) = DOT_PRODUCT(a_cyl, this%esup(:,1))
      this%asup(2) = DOT_PRODUCT(a_cyl, this%esup(:,2))
      this%asup(3) = DOT_PRODUCT(a_cyl, this%esup(:,3))

      DEALLOCATE(cosz)
      DEALLOCATE(sinz)

      ALLOCATE(cosz(mnmax_nyq))
      ALLOCATE(sinz(mnmax_nyq))

      cosz = COS(xm_nyq*this%flux(2) - xn_nyq*this%flux(3))
      sinz = SIN(xm_nyq*this%flux(2) - xn_nyq*this%flux(3))

      this%bsup(1) = 0.0
      this%bsup(2) =                                                           &
     &   SUM((w_low_half*bsupumnc(:,i_low_half) +                              &
     &        w_high_half*bsupumnc(:,i_high_half))*cosz)
      this%bsup(3) =                                                           &
     &   SUM((w_low_half*bsupvmnc(:,i_low_half) +                              &
     &        w_high_half*bsupvmnc(:,i_high_half))*cosz)

      this%bsub(1) =                                                           &
     &   SUM((w_low_full*bsubsmns(:,i_low_full) +                              &
     &        w_high_full*bsubsmns(:,i_high_full))*sinz)
      this%bsub(2) =                                                           &
     &   SUM((w_low_half*bsubumnc(:,i_low_half) +                              &
     &        w_high_half*bsubumnc(:,i_high_half))*cosz)
      this%bsub(3) =                                                           &
     &   SUM((w_low_half*bsubvmnc(:,i_low_half) +                              &
     &        w_high_half*bsubvmnc(:,i_high_half))*cosz)

      IF (lasym) THEN
         this%bsup(2) = this%bsup(2) +                                         &
     &      SUM((w_low_half*bsupumns(:,i_low_half) +                           &
     &           w_high_half*bsupumns(:,i_high_half))*sinz)
         this%bsup(3) = this%bsup(3) +                                         &
     &      SUM((w_low_half*bsupvmns(:,i_low_half) +                           &
     &           w_high_half*bsupvmns(:,i_high_half))*sinz)

         this%bsub(1) = this%bsub(1) +                                         &
     &      SUM((w_low_full*bsubsmnc(:,i_low_full) +                           &
     &           w_high_full*bsubsmnc(:,i_high_full))*cosz)
         this%bsub(2) = this%bsub(2) +                                         &
     &      SUM((w_low_half*bsubumns(:,i_low_half) +                           &
     &           w_high_half*bsubumns(:,i_high_half))*sinz)
         this%bsub(3) = this%bsub(3) +                                         &
     &      SUM((w_low_half*bsubvmns(:,i_low_half) +                           &
     &           w_high_half*bsubvmns(:,i_high_half))*sinz)
      END IF

      DEALLOCATE(cosz)
      DEALLOCATE(sinz)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Get the vector potential in cylindrical coordinates.
!>
!>  Get the vector potential in cylindrical coordinates.
!>
!>  @param[in] this A @ref equilibrium_class instance.
!-------------------------------------------------------------------------------
      FUNCTION vmec_coordinates_get_A_cyl(this)
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3) :: vmec_coordinates_get_A_cyl
      TYPE (vmec_coordinates_class), INTENT(in) :: this

!  Start of executable code
      vmec_coordinates_get_A_cyl(1) =                                          &
     &   DOT_PRODUCT(this%asub, this%esup(1,:))
      vmec_coordinates_get_A_cyl(2) =                                          &
     &   DOT_PRODUCT(this%asub, this%esup(2,:))
      vmec_coordinates_get_A_cyl(3) =                                          &
     &   DOT_PRODUCT(this%asub, this%esup(3,:))

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the magnetic field vector in cylindrical coordinates.
!>
!>  Get the magnetic field vector in cylindrical coordinates.
!>
!>  @param[in] this A @ref equilibrium_class instance.
!-------------------------------------------------------------------------------
      FUNCTION vmec_coordinates_get_B_cyl(this)
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3) :: vmec_coordinates_get_B_cyl
      TYPE (vmec_coordinates_class), INTENT(in) :: this

!  Start of executable code
      vmec_coordinates_get_B_cyl(1) =                                          &
     &   DOT_PRODUCT(this%bsup, this%esub(1,:))
      vmec_coordinates_get_B_cyl(2) =                                          &
     &   DOT_PRODUCT(this%bsup, this%esub(2,:))
      vmec_coordinates_get_B_cyl(3) =                                          &
     &   DOT_PRODUCT(this%bsup, this%esub(3,:))

      END FUNCTION

!*******************************************************************************
!  PRIVATE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief a X b
!>
!>  There is a crossproduct subroutine in math_utilities but a function will
!>  work better.
!>
!>  @param[in] a Vector A
!>  @param[in] b Vector B
!>  @returns A X B
!-------------------------------------------------------------------------------
      FUNCTION cross_product(a, b)
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: cross_product
      REAL (rprec), DIMENSION(3), INTENT(in) :: A
      REAL (rprec), DIMENSION(3), INTENT(in) :: B

!  Start of executable code
      cross_product(1) = a(2)*b(3) - a(3)*b(2)
      cross_product(2) = a(3)*b(1) - a(1)*b(3)
      cross_product(3) = a(1)*b(2) - a(2)*b(1)

      END FUNCTION

      END MODULE
