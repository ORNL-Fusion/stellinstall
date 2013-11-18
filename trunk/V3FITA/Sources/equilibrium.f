!*******************************************************************************
!>  @file equilibrium.f
!>  @brief Contains module @ref equilibrium
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref equilibrium_class. Equilibrium is
!>  an abstract interface to the equilibrium solver. Every method is meant to be
!>  overwritten by a subclass method with a few exceptions. For instance a
!>  vacuum equilibrium could use the default return values for
!>  @ref equilibrium_get_ne, equilibrium_get_te and equilibrium_get_sxrem
!>  methods.
!>  @par Sub Classes:
!>  @ref vmec_equilibrium
!*******************************************************************************

      MODULE equilibrium
      USE vmec_equilibrium
      USE vacuum_equilibrium

      IMPLICIT NONE

!*******************************************************************************
!  equilibrium module parameters
!*******************************************************************************
!>  Type descriptor for the @ref equilibrium_class no subclass.
      INTEGER, PARAMETER :: equilibrium_no_type     = -1
!>  Type descriptor for the @ref equilibrium_class subclass
!>  @ref vmec_equilibrium.
      INTEGER, PARAMETER :: equilibrium_vmec_type   = 0
!>  Type descriptor for the @ref equilibrium_class subclass
!>  @ref vacuum_equilibrium.
      INTEGER, PARAMETER :: equilibrium_vacuum_type = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) equilibrium base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a equilibrium.
!>  @par Sub Classes:
!>  @ref vmec_equilibrium
!-------------------------------------------------------------------------------
      TYPE equilibrium_class
!>  Type descirptor of the equilibrium subclass type.
!>  @par Possible values are:
!>  * @ref equilibrium_no_type
!>  * @ref equilibrium_vmec_type
         INTEGER                      :: type = equilibrium_no_type
!>  Controls if an equilibrium needs to be resolved or not.
         LOGICAL                      :: not_converged
!>  Controls if an equilibrium is forces to be resolved or not.
         LOGICAL                      :: force_solve
!>  An instance of a @ref vmec_equilibrium object.
         TYPE (vmec_class), POINTER   :: vmec => null()
!>  An instance of a @ref vacuum_equilibrium object.
         TYPE (vacuum_class), POINTER :: vacuum => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the constriction of @ref equilibrium_class types using
!>  @ref equilibrium_construct_vmec
!-------------------------------------------------------------------------------
      INTERFACE equilibrium_construct
         MODULE PROCEDURE equilibrium_construct_vmec,                          &
     &                    equilibrium_construct_vacuum
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref equilibrium_class containing a @ref vmec_equilibrium
!>  object.
!>
!>  Allocates memory and initializes a @ref equilibrium_class object.
!>  @param[in] vmec_object An instance of a @ref vmec_equilibrium subclass.
!>  @param[in] force_solve If true, forces the equilbirum to resolve every time.
!>  @returns A pointer to a constructed @ref equilibrium_class object.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_construct_vmec(vmec_object, force_solve)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), POINTER :: equilibrium_construct_vmec
      TYPE (vmec_class), POINTER, INTENT(in) :: vmec_object
      LOGICAL, INTENT(in)                    :: force_solve

!  Start of executable code
      ALLOCATE(equilibrium_construct_vmec)

      equilibrium_construct_vmec%type = equilibrium_vmec_type
      equilibrium_construct_vmec%vmec => vmec_object
      equilibrium_construct_vmec%not_converged = .true.
      equilibrium_construct_vmec%force_solve = force_solve

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref equilibrium_class containing a
!>  @ref vacuum_equilibrium object.
!>
!>  Allocates memory and initializes a @ref equilibrium_class object.
!>
!>  @param[in] vacuum_object An instance of a @ref vacuum_equilibrium subclass.
!>  @param[in] force_solve   If true, forces the equilbirum to resolve every
!>                           time.
!>  @returns A pointer to a constructed @ref equilibrium_class object.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_construct_vacuum(vacuum_object, force_solve)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), POINTER :: equilibrium_construct_vacuum
      TYPE (vacuum_class), POINTER, INTENT(in) :: vacuum_object
      LOGICAL, INTENT(in)                      :: force_solve

!  Start of executable code
      ALLOCATE(equilibrium_construct_vacuum)

      equilibrium_construct_vacuum%type = equilibrium_vacuum_type
      equilibrium_construct_vacuum%vacuum => vacuum_object
      equilibrium_construct_vacuum%not_converged = .false.
      equilibrium_construct_vacuum%force_solve = force_solve

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref equilibrium_class object.
!>
!>  Deallocates memory and uninitializes a @ref equilibrium_class object.
!>
!>  @param[inout] this A @ref equilibrium_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%vmec)) THEN
         CALL vmec_destruct(this%vmec)
         this%vmec => null()
      END IF

      IF (ASSOCIATED(this%vacuum)) THEN
         CALL vacuum_destruct(this%vacuum)
         this%vacuum => null()
      END IF

      this%type = equilibrium_no_type
      this%not_converged = .true.
      this%force_solve = .false.

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sets the value of a reconstruction equilibrium parameter.
!>
!>  This method is virtual. The actual setting of the parameter should be
!>  handled by a subclass method. The subclass is responcible for informing the
!>  equilibirum if the equilibrium needs to be reconverged/resolved.
!>  @see vmec_equilibrium::vmec_set_param
!>  @see vacuum_equilibrium::vacuum_set_param
!>
!>  @param[inout] this    A @ref equilibrium_class instance.
!>  @param[in]    id      ID of the parameter.
!>  @param[in]    i_index The ith index of the parameter.
!>  @param[in]    j_index The jth index of the parameter.
!>  @param[in]    value   The value of the parameter.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_set_param(this, id, i_index, j_index,             &
     &                                 value)

      IMPLICIT NONE

!  Declare Sruments
      TYPE (equilibrium_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: id
      INTEGER, INTENT(in)                     :: i_index
      INTEGER, INTENT(in)                     :: j_index
      REAL (rprec), INTENT(in)                :: value

!  Start of executable code
      SELECT CASE (this%type)

!  Set param needs to come first so that function gets executed. Otherwise, when
!  this%not_converged was true, the expression will evaluate to false all the
!  time and not bother calling vmec_set_param.
         CASE (equilibrium_vmec_type)
            this%not_converged = vmec_set_param(this%vmec, id,                 &
     &                                          i_index, j_index,              &
     &                                          value)                         &
     &                         .or. this%not_converged

         CASE (equilibrium_vacuum_type)
            this%not_converged = vacuum_set_param(this%vacuum, id,             &
     &                                            i_index, j_index,            &
     &                                            value)                       &
     &                         .or. this%not_converged

      END SELECT

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Sets the magnetic cache of the equilibrium.
!>
!>  This method is virtual. The actual setting of the magnetic cache should be
!>  handled a subclass method. The magnetic cache is any pre compuation that
!>  needs to be performed so that magnetic signals may compute a modeled signal.
!>  If the equilibrium does not require any information to be cached, this
!>  method does not need to be overridden. This method is meant to only be
!>  called once to tell the equilibrium it needs to use a magnetic cache.
!>  @see vmec_equilibrium::vmec_set_magnetic_cache
!>
!>  @param[inout] this             A @ref equilibrium_class instance.
!>  @param[in]    total_tor_planes Total number of toroidal planes.
!>  @param[in]    tor_planes_used  Number of stellarator symmetric toroidal
!>                                 planes.
!>  @note Stellarator symmetry may be exploited speed up computation and reduce
!>  memory usage. If stellarator symmetry is not used then tor_planes_used
!>  should match total_tor_planes.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_set_magnetic_cache(this, total_tor_planes,        &
     &                                          tor_planes_used)

      IMPLICIT NONE

!  Declare Sruments
      TYPE (equilibrium_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: total_tor_planes
      INTEGER, INTENT(in)                     :: tor_planes_used

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_set_magnetic_cache(this%vmec, total_tor_planes,          &
     &                                   tor_planes_used)

      END SELECT

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Get the id for a reconstruction parameter.
!>
!>  This method is virtual. The actual getting of the reconstruction parameter
!>  id should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_param_id
!>  @see vacuum_equilibrium::vacuum_get_param_id
!>
!>  @param[in] this       A @ref equilibrium_class instance.
!>  @param[in] param_name Name of a reconstruction parameter.
!>  @returns The id for a reconstruction parameter.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_param_id(this, param_name)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: equilibrium_get_param_id
      TYPE (equilibrium_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)        :: param_name

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_param_id =                                         &
     &         vmec_get_param_id(this%vmec, param_name)

         CASE (equilibrium_vacuum_type)
            equilibrium_get_param_id =                                         &
     &         vacuum_get_param_id(this%vacuum, param_name)

         CASE DEFAULT
            equilibrium_get_param_id = data_no_id

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the value of a reconstruction equilibrium parameter.
!>
!>  This method is virtual. The actual getting of the reconstruction parameter
!>  value should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_param_value
!>  @see vacuum_equilibrium::vacuum_get_param_value
!>
!>  @param[in] this    A @ref equilibrium_class instance.
!>  @param[in] id      ID of the parameter.
!>  @param[in] i_index The ith index of the parameter.
!>  @param[in] j_index The jth index of the parameter.
!>  @returns The value of the parameter.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_param_value(this, id, i_index, j_index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_param_value
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id
      INTEGER, INTENT(in)                  :: i_index
      INTEGER, INTENT(in)                  :: j_index

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_param_value =                                      &
     &         vmec_get_param_value(this%vmec, id, i_index, j_index)

         CASE (equilibrium_vacuum_type)
            equilibrium_get_param_value =                                      &
     &         vacuum_get_param_value(this%vacuum, id, i_index, j_index)

         CASE DEFAULT
            equilibrium_get_param_value = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the name of a reconstruction equilibrium parameter.
!>
!>  This method is virtual. The actual getting of the reconstruction parameter
!>  name should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_param_name
!>  @see vacuum_equilibrium::vacuum_get_param_name
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns The name of the parameter.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_param_name(this, id)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=data_name_length) :: equilibrium_get_param_name
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_param_name =                                       &
     &         vmec_get_param_name(this%vmec, id)

         CASE (equilibrium_vacuum_type)
            equilibrium_get_param_name =                                       &
     &         vacuum_get_param_name(this%vacuum, id)

         CASE DEFAULT
            equilibrium_get_param_name = ''

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron density at a position.
!>
!>  This method is virtual. The actual getting of the electron density should be
!>  handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_ne
!>
!>  @param[in] this   A @ref equilibrium_class instance.
!>  @param[in] x_cart Cartesian position to get the electron density at.
!>  @returns The electron density at x_cart.
!>  @note Note the electron density is not relevant to all equilibrium types.
!>  For instance, vacuum equilibrium should not override this method.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_ne(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_ne
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_ne = vmec_get_ne(this%vmec, x_cart)

         CASE DEFAULT
            equilibrium_get_ne = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the electron temperature at a position.
!>
!>  This method is virtual. The actual getting of the electron temperature
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_te
!>
!>  @param[in] this   A @ref equilibrium_class instance.
!>  @param[in] x_cart Cartesian position to get the electron temperature at.
!>  @returns The electron temperature at x_cart.
!>  @note Note the electron temperature is not relevant to all equilibrium
!>  types. For instance, vacuum equilibrium should not override this method.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_te(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_te
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_te = vmec_get_te(this%vmec, x_cart)

         CASE DEFAULT
            equilibrium_get_te = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the soft x-ray emissivity at a position.
!>
!>  This method is virtual. The actual getting of the electron temperature
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_sxrem
!>
!>  @param[in] this   A @ref equilibrium_class instance.
!>  @param[in] x_cart Cartesian position to get the soft x-ray emissivity at.
!>  @returns The soft x-ray emissivity at x_cart.
!>  @note Note the soft x-ray emissivity is not relevant to all equilibrium
!>  types. For instance, vacuum equilibrium should not override this method.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_sxrem(this, x_cart, index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_sxrem
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      INTEGER, INTENT(in)                    :: index

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_sxrem = vmec_get_sxrem(this%vmec, x_cart,          &
     &                                             index)

         CASE DEFAULT
            equilibrium_get_sxrem = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the plasma pressure at a position.
!>
!>  This method is virtual. The actual getting of the plasma pressure should be
!>  handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_p
!>
!>  @param[in] this   A @ref equilibrium_class instance.
!>  @param[in] x_cart Cartesian position to get the plasma pressure at.
!>  @returns The plasma pressure at x_cart.
!>  @note Note the plasma pressure is not relevant to all equilibrium types. For
!>  instance, vacuum equilibrium should not override this method.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_p(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_p
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_p = vmec_get_p(this%vmec, x_cart)

         CASE DEFAULT
            equilibrium_get_p = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the magnetic field vector at a position.
!>
!>  This method is virtual. The actual getting of the magnetic field vector
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_B_vec
!>  @see vacuum_equilibrium::vacuum_get_B_vec
!>
!>  @param[in] this   A @ref equilibrium_class instance.
!>  @param[in] x_cart Cartesian position to get the magnetic field vector at.
!>  @returns The magnetic field vector at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_B_vec(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3) :: equilibrium_get_B_vec
      TYPE (equilibrium_class), INTENT(in)   :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_B_vec = vmec_get_B_vec(this%vmec, x_cart)

         CASE (equilibrium_vacuum_type)
            equilibrium_get_B_vec = vacuum_get_B_vec(this%vacuum,              &
     &                                               x_cart)

         CASE DEFAULT
            equilibrium_get_B_vec = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the loop integrated magnetic field at a position.
!>
!>  This method is virtual. The actual getting of the loop integrated magnetic
!>  field should be handled by a subclass method. Returns the closed loop
!>  integration of the magnetic field to determine the enclosed current.
!>  @see vmec_equilibrium::vmec_get_Int_B_dphi
!>  @see vacuum_equilibrium::vacuum_get_Int_B_dphi
!>
!>  @param[in] this  A @ref equilibrium_class instance.
!>  @param[in] r     Radial position to integrate about.
!>  @param[in] theta Poloidal angle to integrate about.
!>  @returns The loop integrated magnetic field at x_cart.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_Int_B_dphi(this, r, theta)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_Int_B_dphi
      TYPE (equilibrium_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)             :: r
      REAL (rprec), INTENT(in)             :: theta

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_Int_B_dphi =                                       &
     &         vmec_get_Int_B_dphi(this%vmec, r, theta)

         CASE (equilibrium_vacuum_type)
            equilibrium_get_Int_B_dphi =                                       &
     &         vacuum_get_Int_B_dphi(this%vacuum, r, theta)

         CASE DEFAULT
            equilibrium_get_Int_B_dphi = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the r and z positions of the outer surface at a toroidal angle.
!>
!>  This method is virtual. The actual getting of the magnetic field vector
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_plasma_edge
!>
!>  @param[in]  this A @ref equilibrium_class instance.
!>  @param[in]  phi  Toroidal angle to determine the outer surface at.
!>  @param[out] r    The radial postions of the other surface in a single
!>                   toroidal angle.
!>  @param[out] z    The Z postions of the other surface in a single toroidal
!>                   angle.
!>  @returns The number of elements in the r and z arrays.
!>
!>  @todo An "edge" could be defined as a the last field line that intersects a
!>  limiter for the vacuum case. By defineing a r_last parameter, a feild line
!>  can be followed around the machine. When that field line intersecs the
!>  desired phi plane, record the r and z positions. Field line following will
!>  also need to be performed for siesta.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_plasma_edge(this, phi, r, z)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: equilibrium_get_plasma_edge
      TYPE (equilibrium_class), INTENT(in)             :: this
      REAL (rprec), INTENT (in)                        :: phi
      REAL (rprec), DIMENSION(:), POINTER, INTENT(out) :: r
      REAL (rprec), DIMENSION(:), POINTER, INTENT(out) :: z

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_plasma_edge =                                      &
     &         vmec_get_plasma_edge(this%vmec, phi, r, z)

         CASE DEFAULT
            equilibrium_get_plasma_edge = 0
            r => null()
            z => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration radial grid points.
!>
!>  This method is virtual. The actual getting of the radial grid points
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_magnetic_volume_rgrid
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @returns The radial grid points. 
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_rgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_rgrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_rgrid =>                           &
     &         vmec_get_magnetic_volume_rgrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_rgrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration z grid points.
!>
!>  This method is virtual. The actual getting of the z grid points should be
!>  handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_magnetic_volume_zgrid
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @returns The z grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_zgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_zgrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_zgrid =>                           &
     &         vmec_get_magnetic_volume_zgrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_zgrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration jr grid points.
!>
!>  This method is virtual. The actual getting of the jr grid points should be
!>  handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_magnetic_volume_jrgrid
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @returns The jr grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_jrgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_jrgrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_jrgrid =>                          &
     &         vmec_get_magnetic_volume_jrgrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_jrgrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration jphi grid points.
!>
!>  This method is virtual. The actual getting of the jphi grid points should be
!>  handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_magnetic_volume_jphigrid
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @returns The jphi grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_jphigrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_jphigrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_jphigrid =>                        &
     &         vmec_get_magnetic_volume_jphigrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_jphigrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume magnetic volume integration jz grid points.
!>
!>  This method is virtual. The actual getting of the jz grid points should be
!>  handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_magnetic_volume_jzgrid
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @returns The jz grid points.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_magnetic_volume_jzgrid(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:,:), POINTER ::                               &
     &   equilibrium_get_magnetic_volume_jzgrid
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_magnetic_volume_jzgrid =>                          &
     &         vmec_get_magnetic_volume_jzgrid(this%vmec)

         CASE DEFAULT
            equilibrium_get_magnetic_volume_jzgrid => null()

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get volume integration element.
!>
!>  This method is virtual. The actual getting of the volume integration element
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_get_magnetic_volume_int_element
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @returns The volume integration element.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_volume_int_element(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: equilibrium_get_volume_int_element
      TYPE (equilibrium_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)

        CASE (equilibrium_vmec_type)
            equilibrium_get_volume_int_element =                               &
     &         vmec_get_volume_int_element(this%vmec)

         CASE DEFAULT
            equilibrium_get_volume_int_element = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get external current.
!>
!>  This method is virtual. The actual getting of the external current should be
!>  handled by a subclass method. To void extra memory operations, arrays are
!>  returned as pointers. The memory management is handled by the equilibrium
!>  subclass instance. If the equilibrium has no external currents return a
!>  null() pointer or don't override this method. MAKEGRID has a mode where
!>  currents maybe scaled. Inform the signal using this to scale the currents
!>  if need be.
!>  @see vmec_equilibrium::vmec_get_ext_currents
!>  @see vacuum_equilibrium::vacuum_get_ext_currents
!>
!>  @param[in]  this           A @ref equilibrium_class instance.
!>  @param[out] scale_currents Informs the caller that currents need to be
!>                             scaled.
!>  @returns The external currents.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_get_ext_currents(this, scale_currents)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:), POINTER  ::                                  &
     &   equilibrium_get_ext_currents
      TYPE (equilibrium_class), INTENT(in) :: this
      LOGICAL, INTENT(out)                 :: scale_currents

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_get_ext_currents =>                                    &
     &         vmec_get_ext_currents(this%vmec, scale_currents)

         CASE (equilibrium_vacuum_type)
            equilibrium_get_ext_currents =>                                    &
     &         vacuum_get_ext_currents(this%vacuum, scale_currents)

         CASE DEFAULT
            scale_currents = .false.
            equilibrium_get_ext_currents => null()

      END SELECT

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Checks if a parameter id is a scaler value.
!>
!>  This method is virtual. The actual check should be handled by a subclass
!>  method.
!>  @see vmec_equilibrium::vmec_is_scaler_value
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the parameter is a scaler and false if otherwise.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_scaler_value(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_scaler_value
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_scaler_value =                                      &
     &         vmec_is_scaler_value(this%vmec, id)

         CASE DEFAULT
            equilibrium_is_scaler_value = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if a parameter id is a 1d array.
!>
!>  This method is virtual. The actual check should be handled by a subclass
!>  method.
!>  @see vmec_equilibrium::vmec_is_1d_array
!>  @see vacuum_equilibrium::vacuum_is_1d_array
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the parameter is a 1d array and false if otherwise.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_1d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_1d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_1d_array = vmec_is_1d_array(this%vmec, id)

         CASE (equilibrium_vacuum_type)
            equilibrium_is_1d_array = vacuum_is_1d_array(this%vacuum,          &
     &                                                   id)

         CASE DEFAULT
            equilibrium_is_1d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if a parameter id is a 2d array.
!>
!>  This method is virtual. The actual check should be handled by a subclass
!>  method.
!>  @see vmec_equilibrium::vmec_is_2d_array
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the parameter is a 2d array and false if otherwise.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_2d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_2d_array
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_2d_array = vmec_is_2d_array(this%vmec, id)

         CASE DEFAULT
            equilibrium_is_2d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if a parameter id is a reconstruction parameter.
!>
!>  This method is virtual. The actual check should be handled by a subclass
!>  method.
!>  @see vmec_equilibrium::vmec_is_recon_param
!>  @see vacuum_equilibrium::vacuum_is_recon_param
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @param[in] id   ID of the parameter.
!>  @returns True if the parameter is a reconstruction parameter and false if
!>  otherwise.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_is_recon_param(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: equilibrium_is_recon_param
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            equilibrium_is_recon_param =                                       &
     &         vmec_is_recon_param(this%vmec, id)

         CASE (equilibrium_vacuum_type)
            equilibrium_is_recon_param =                                       &
     &         vacuum_is_recon_param(this%vacuum, id)

         CASE DEFAULT
            equilibrium_is_recon_param = .false.

      END SELECT

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Solves the equilibrium.
!>
!>  This method is virtual. The actual solving of the equilibrium should be
!>  handled by a subclass method. The equilibrium is only reconverged if it is
!>  in an unconverged state. If the equilibrium was converged, this sets the
!>  @ref equilibrium_class::not_converged to false to inform future calls that
!>  the equilibirum is already in a converged state.
!>  @see vmec_equilibrium::vmec_converge
!>
!>  @param[inout] this     A @ref equilibrium_class instance.
!>  @param[inout] num_iter Counter to track the number of iterations.
!>  @param[in]    iou      Input/output unit of the file to write logs to.
!>  @returns True if the convergece was sucessful and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION equilibrium_converge(this, num_iter, iou)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                 :: equilibrium_converge
      TYPE (equilibrium_class), INTENT(inout) :: this
      INTEGER, INTENT(inout)                  :: num_iter
      INTEGER, INTENT(in)                     :: iou

!  local variables
      INTEGER                                 :: init_num_iter

!  Start of executable code
      init_num_iter = num_iter
      IF (this%not_converged .or. this%force_solve) THEN
         SELECT CASE (this%type)

            CASE (equilibrium_vmec_type)
               equilibrium_converge = vmec_converge(this%vmec, num_iter,       &
     &                                              iou)
               this%not_converged = .not.equilibrium_converge

            CASE (equilibrium_vacuum_type)
               equilibrium_converge = .true.
               this%not_converged = .false.

            CASE DEFAULT
               equilibrium_converge = .false.
               this%not_converged = .true.

         END SELECT
      ELSE
         equilibrium_converge = .true.
      ENDIF

      WRITE (*,1000) equilibrium_converge, num_iter,                           &
     &               num_iter - init_num_iter
      WRITE (iou,1000) equilibrium_converge, num_iter,                         &
     &                 num_iter - init_num_iter
1000  FORMAT('Equilibrium converged ',l,2(2x,i7))

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Save the internal state of the equilibrium.
!>
!>  This method is virtual. The actual saving of the equilibrium state should be
!>  handled by a subclass method. 
!>  @see vmec_equilibrium::vmec_converge
!>
!>  @param[inout] this A @ref equilibrium_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_save_state(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(inout) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_save_state(this%vmec)

      END SELECT

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Reset the internal state of the equilibrium.
!>
!>  This method is virtual. The actual reseting of the equilibrium state should
!>  be handled by a subclass method. Sub classes need to make sure that the
!>  equilibrium is reset to a converged state.
!>  @see vmec_equilibrium::vmec_converge
!>
!>  @param[inout] this A @ref equilibrium_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_reset_state(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(inout) :: this

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_reset_state(this%vmec)

      END SELECT

!  Reseting should set the equilibrium to a converged state.
      this%not_converged = .false.

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out the equilibrium to an output file.
!>
!>  This method is virtual. The actual writing of the equilibrium should be
!>  handled by a subclass method.
!>  @see vmec_equilibrium::vmec_write
!>  @see vacuum_equilibrium::vacuum_write
!>
!>  @param[in] this A @ref equilibrium_class instance.
!>  @param[in] iou  Input/output unit of the output file.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_write(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: iou

!  Start of executable code
      SELECT CASE (this%type)

      CASE (equilibrium_vmec_type)
         CALL vmec_write(this%vmec, iou)

      CASE (equilibrium_vacuum_type)
         CALL vacuum_write(this%vacuum, iou)

      END SELECT

      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Define NetCDF variables for the result file
!>
!>  This method is virtual. The actual defining of the equilibrium variables
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_def_result
!>
!>  @param[in] this             A @ref equilibrium_class instance.
!>  @param[in] result_ncid      NetCDF file id of the result file.
!>  @param[in] maxnsetps_dim_id NetCDF dimension id of the number of steps
!>                              dimension.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_def_result(this, result_ncid,                     &
     &                                  maxnsetps_dim_id)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: result_ncid
      INTEGER, INTENT(in)                  :: maxnsetps_dim_id

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_def_result(this%vmec, result_ncid,                       &
     &                           maxnsetps_dim_id)

      END SELECT

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write inital data to NetCDF result file
!>
!>  This method is virtual. The actual writing of the equilibrium variables
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_write_init_data
!>
!>  @param[in] this        A @ref equilibrium_class instance.
!>  @param[in] result_ncid NetCDF file id of the result file.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_write_init_data(this, result_ncid)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: result_ncid

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_write_init_data(this%vmec, result_ncid)

      END SELECT

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write step data to NetCDF result file
!>
!>  This method is virtual. The actual writing of the equilibrium variables
!>  should be handled by a subclass method.
!>  @see vmec_equilibrium::vmec_write_step_data
!>
!>  @param[in] this         A @ref equilibrium_class instance.
!>  @param[in] iou          NetCDF file id of the result file.
!>  @param[in] current_step Step index to write variables to.
!-------------------------------------------------------------------------------
      SUBROUTINE equilibrium_write_step_data(this, result_ncid,                &
     &                                       current_step)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (equilibrium_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: result_ncid
      INTEGER, INTENT(in)                  :: current_step

!  Start of executable code
      SELECT CASE (this%type)

         CASE (equilibrium_vmec_type)
            CALL vmec_write_step_data(this%vmec, result_ncid,                  &
     &                                current_step)

      END SELECT

      END SUBROUTINE

      END MODULE
