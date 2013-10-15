!*******************************************************************************
!>  @file wout_file_cache.f
!>  @brief Contains module @ref wout_file_cache.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref wout_file_cache_class. This caches
!>  the data in from the wout file.
!>  @note The @ref read_wout_mod cannot be used in the full file context here to
!>  avoid variable name clashes with the phi variable in @ref magnetic.
!*******************************************************************************

      MODULE wout_file_cache
      USE stel_kinds

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) wout_file_cache class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a wout_file_cache. This contains a copy of every
!>  variable that is needed from the wout file to compute a modeled signal.
!>  Fixed parameters do not need to be copied.
!-------------------------------------------------------------------------------
      TYPE wout_file_cache_class
!>  cosmn component of cylindrical R, full mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: rmnc => null()
!>  sinmn component of cylindrical R, full mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: rmns => null()
!>  cosmn component of cylindrical Z, full mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: zmnc => null()
!>  sinmn component of cylindrical Z, full mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: zmns => null()
!>  sinmn component of lambda, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: lmns => null()
!>  cosmn component of lambda, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: lmnc => null()
!>  cosmn component of jacobian, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: gmnc => null()
!>  cosmn component of jacobian, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: gmns => null()
!>  sinmn covariant s-component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsubsmns => null()
!>  cosmn covariant s-component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsubsmnc => null()
!>  cosmn covariant u-component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsubumnc => null()
!>  sinmn covariant u-component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsubumns => null()
!>  cosmn covariant v-component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsubvmnc => null()
!>  sinmn covariant v-component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsubvmns => null()
!>  cosmn covariant u^component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsupumnc => null()
!>  sinmn covariant u^component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsupumns => null()
!>  cosmn covariant v^component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsupvmnc => null()
!>  sinmn covariant v^component of B, half mesh.
         REAL (rprec), DIMENSION(:,:), POINTER :: bsupvmns => null()
!>  Toroidal flux on full mesh.
         REAL (rprec), DIMENSION(:), POINTER   :: phi => null()
!>  Rotational transform on full mesh.
         REAL (rprec), DIMENSION(:), POINTER   :: iotaf => null()
!>  Rotational transform on half mesh.
         REAL (rprec), DIMENSION(:), POINTER   :: iotas => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the constriction of @ref wout_file_cache_class types using
!>  @ref wout_file_cache_construct_new or @ref wout_file_cache_construct_copy
!-------------------------------------------------------------------------------
      INTERFACE wout_file_cache_construct
         MODULE PROCEDURE wout_file_cache_construct_new,                       &
     &                    wout_file_cache_construct_copy
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a new @ref wout_file_cache_class object.
!>
!>  Allocates memory and initializes a @ref wout_file_cache_class object.
!>  @returns A pointer to a constructed @ref wout_file_cache_class object.
!-------------------------------------------------------------------------------
      FUNCTION wout_file_cache_construct_new()
      USE read_wout_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (wout_file_cache_class), POINTER ::                                 &
     &   wout_file_cache_construct_new

!  Start of executable code
      ALLOCATE(wout_file_cache_construct_new)

      ALLOCATE(wout_file_cache_construct_new%rmnc(mnmax,ns))
      ALLOCATE(wout_file_cache_construct_new%zmns(mnmax,ns))
      ALLOCATE(wout_file_cache_construct_new%lmns(mnmax,ns))
      ALLOCATE(wout_file_cache_construct_new%gmnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_new%bsubumnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_new%bsubvmnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_new%bsubsmns(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_new%bsupumnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_new%bsupvmnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_new%phi(ns))
      ALLOCATE(wout_file_cache_construct_new%iotaf(ns))
      ALLOCATE(wout_file_cache_construct_new%iotas(ns))

      wout_file_cache_construct_new%rmnc = rmnc
      wout_file_cache_construct_new%zmns = zmns
      wout_file_cache_construct_new%lmns = lmns
      wout_file_cache_construct_new%gmnc = gmnc
      wout_file_cache_construct_new%bsubumnc = bsubumnc
      wout_file_cache_construct_new%bsubvmnc = bsubvmnc
      wout_file_cache_construct_new%bsubsmns = bsubsmns
      wout_file_cache_construct_new%bsupumnc = bsupumnc
      wout_file_cache_construct_new%bsupvmnc = bsupvmnc
      wout_file_cache_construct_new%phi = phi
      wout_file_cache_construct_new%iotaf = iotaf
      wout_file_cache_construct_new%iotas = iotas

      IF (lasym) THEN
         ALLOCATE(wout_file_cache_construct_new%rmns(mnmax,ns))
         ALLOCATE(wout_file_cache_construct_new%zmnc(mnmax,ns))
         ALLOCATE(wout_file_cache_construct_new%lmnc(mnmax,ns))
         ALLOCATE(wout_file_cache_construct_new%gmns(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_new%bsubumns(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_new%bsubvmns(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_new%bsubsmnc(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_new%bsupumns(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_new%bsupvmns(mnmax_nyq,ns))

         wout_file_cache_construct_new%rmns = rmns
         wout_file_cache_construct_new%zmnc = zmnc
         wout_file_cache_construct_new%lmnc = lmnc
         wout_file_cache_construct_new%gmns = gmns
         wout_file_cache_construct_new%bsubumns = bsubumns
         wout_file_cache_construct_new%bsubvmns = bsubvmns
         wout_file_cache_construct_new%bsubsmnc = bsubsmnc
         wout_file_cache_construct_new%bsupumns = bsupumns
         wout_file_cache_construct_new%bsupvmns = bsupvmns
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a new @ref wout_file_cache_class object by copying an old
!>  one.
!>
!>  Allocates memory and initializes a @ref wout_file_cache_class object.
!>  @returns A pointer to a constructed @ref wout_file_cache_class object.
!-------------------------------------------------------------------------------
      FUNCTION wout_file_cache_construct_copy(cache_object)
      USE read_wout_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (wout_file_cache_class), POINTER ::                                 &
     &   wout_file_cache_construct_copy
      TYPE (wout_file_cache_class), INTENT(in) :: cache_object

!  Start of executable code
      ALLOCATE(wout_file_cache_construct_copy)

      ALLOCATE(wout_file_cache_construct_copy%rmnc(mnmax,ns))
      ALLOCATE(wout_file_cache_construct_copy%zmns(mnmax,ns))
      ALLOCATE(wout_file_cache_construct_copy%lmns(mnmax,ns))
      ALLOCATE(wout_file_cache_construct_copy%gmnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_copy%bsubumnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_copy%bsubvmnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_copy%bsubsmns(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_copy%bsupumnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_copy%bsupvmnc(mnmax_nyq,ns))
      ALLOCATE(wout_file_cache_construct_copy%phi(ns))
      ALLOCATE(wout_file_cache_construct_copy%iotaf(ns))
      ALLOCATE(wout_file_cache_construct_copy%iotas(ns))

      wout_file_cache_construct_copy%rmnc = cache_object%rmnc
      wout_file_cache_construct_copy%zmns = cache_object%zmns
      wout_file_cache_construct_copy%lmns = cache_object%lmns
      wout_file_cache_construct_copy%gmnc = cache_object%gmnc
      wout_file_cache_construct_copy%bsubumnc = cache_object%bsubumnc
      wout_file_cache_construct_copy%bsubvmnc = cache_object%bsubvmnc
      wout_file_cache_construct_copy%bsubsmns = cache_object%bsubsmns
      wout_file_cache_construct_copy%bsupumnc = cache_object%bsupumnc
      wout_file_cache_construct_copy%bsupvmnc = cache_object%bsupvmnc
      wout_file_cache_construct_copy%phi = cache_object%phi
      wout_file_cache_construct_copy%iotaf = cache_object%iotaf
      wout_file_cache_construct_copy%iotas = cache_object%iotas

      IF (lasym) THEN
         ALLOCATE(wout_file_cache_construct_copy%rmns(mnmax,ns))
         ALLOCATE(wout_file_cache_construct_copy%zmnc(mnmax,ns))
         ALLOCATE(wout_file_cache_construct_copy%lmnc(mnmax,ns))
         ALLOCATE(wout_file_cache_construct_copy%gmns(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_copy%bsubumns(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_copy%bsubvmns(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_copy%bsubsmnc(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_copy%bsupumns(mnmax_nyq,ns))
         ALLOCATE(wout_file_cache_construct_copy%bsupvmns(mnmax_nyq,ns))

         wout_file_cache_construct_copy%rmns = cache_object%rmns
         wout_file_cache_construct_copy%zmnc = cache_object%zmnc
         wout_file_cache_construct_copy%lmnc = cache_object%lmnc
         wout_file_cache_construct_copy%gmns = cache_object%gmns
         wout_file_cache_construct_copy%bsubumns = cache_object%bsubumns
         wout_file_cache_construct_copy%bsubvmns = cache_object%bsubvmns
         wout_file_cache_construct_copy%bsubsmnc = cache_object%bsubsmnc
         wout_file_cache_construct_copy%bsupumns = cache_object%bsupumns
         wout_file_cache_construct_copy%bsupvmns = cache_object%bsupvmns
      END IF

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref wout_file_cache_class object.
!>
!>  Deallocates memory and uninitializes a @ref wout_file_cache_class object.
!>
!>  @param[inout] this A @ref wout_file_cache_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE wout_file_cache_destruct(this)
      USE read_wout_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (wout_file_cache_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%rmnc)) THEN
         DEALLOCATE(this%rmnc)
         this%rmnc => null()
      END IF

      IF (ASSOCIATED(this%rmns)) THEN
         DEALLOCATE(this%rmns)
         this%rmns => null()
      END IF

      IF (ASSOCIATED(this%zmnc)) THEN
         DEALLOCATE(this%zmnc)
         this%zmnc => null()
      END IF

      IF (ASSOCIATED(this%zmns)) THEN
         DEALLOCATE(this%zmns)
         this%zmns => null()
      END IF

      IF (ASSOCIATED(this%lmns)) THEN
         DEALLOCATE(this%lmns)
         this%lmns => null()
      END IF

      IF (ASSOCIATED(this%lmnc)) THEN
         DEALLOCATE(this%lmnc)
         this%lmnc => null()
      END IF

      IF (ASSOCIATED(this%gmnc)) THEN
         DEALLOCATE(this%gmnc)
         this%gmnc => null()
      END IF

      IF (ASSOCIATED(this%gmns)) THEN
         DEALLOCATE(this%gmns)
         this%gmns => null()
      END IF

      IF (ASSOCIATED(this%bsubumnc)) THEN
         DEALLOCATE(this%bsubumnc)
         this%bsubumnc => null()
      END IF

      IF (ASSOCIATED(this%bsubvmnc)) THEN
         DEALLOCATE(this%bsubvmnc)
         this%bsubvmnc => null()
      END IF

      IF (ASSOCIATED(this%bsubsmns)) THEN
         DEALLOCATE(this%bsubsmns)
         this%bsubsmns => null()
      END IF

      IF (ASSOCIATED(this%bsupumnc)) THEN
         DEALLOCATE(this%bsupumnc)
         this%bsupumnc => null()
      END IF

      IF (ASSOCIATED(this%bsupvmnc)) THEN
         DEALLOCATE(this%bsupvmnc)
         this%bsupvmnc => null()
      END IF

      IF (ASSOCIATED(this%bsubumns)) THEN
         DEALLOCATE(this%bsubumns)
         this%bsubumns => null()
      END IF

      IF (ASSOCIATED(this%bsubvmns)) THEN
         DEALLOCATE(this%bsubvmns)
         this%bsubvmns => null()
      END IF

      IF (ASSOCIATED(this%bsubsmnc)) THEN
         DEALLOCATE(this%bsubsmnc)
         this%bsubsmnc => null()
      END IF

      IF (ASSOCIATED(this%bsupumns)) THEN
         DEALLOCATE(this%bsupumns)
         this%bsupumns => null()
      END IF

      IF (ASSOCIATED(this%bsupvmns)) THEN
         DEALLOCATE(this%bsupvmns)
         this%bsupvmns => null()
      END IF

      IF (ASSOCIATED(this%phi)) THEN
         DEALLOCATE(this%phi)
         this%phi => null()
      END IF

      IF (ASSOCIATED(this%iotaf)) THEN
         DEALLOCATE(this%iotaf)
         this%iotaf => null()
      END IF

      IF (ASSOCIATED(this%iotas)) THEN
         DEALLOCATE(this%iotas)
         this%iotas => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Save the internal state of the wout file.
!>
!>  Saves the wout file variables.
!>
!>  @param[inout] this A @ref wout_file_cache_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE wout_file_cache_save_state(this)
      USE read_wout_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (wout_file_cache_class), INTENT(inout) :: this

!  Start of executable code
      this%rmnc = rmnc
      this%zmns = zmns
      this%lmns = lmns
      this%gmnc = gmnc
      this%bsubumnc = bsubumnc
      this%bsubvmnc = bsubvmnc
      this%bsubsmns = bsubsmns
      this%bsupumnc = bsupumnc
      this%bsupvmnc = bsupvmnc
      this%phi = phi
      this%iotaf = iotaf
      this%iotas = iotas

      IF (lasym) THEN
         this%rmns = rmns
         this%zmnc = zmnc
         this%lmnc = lmnc
         this%gmns = gmns
         this%bsubumns = bsubumns
         this%bsubvmns = bsubvmns
         this%bsubsmnc = bsubsmnc
         this%bsupumns = bsupumns
         this%bsupvmns = bsupvmns
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Reset the internal state of the wout file.
!>
!>  Copies the saved wout file variables back to @ref read_wout_mod.
!>
!>  @param[inout] this A @ref wout_file_cache_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE wout_file_cache_reset_state(this)
      USE read_wout_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (wout_file_cache_class), INTENT(in) :: this

!  local variables
      INTEGER                                  :: status

!  Start of executable code
      rmnc = this%rmnc
      zmns = this%zmns
      lmns = this%lmns
      gmnc = this%gmnc
      bsubumnc = this%bsubumnc
      bsubvmnc = this%bsubvmnc
      bsubsmns = this%bsubsmns
      bsupumnc = this%bsupumnc
      bsupvmnc = this%bsupvmnc
      phi = this%phi
      iotaf = this%iotaf
      iotas = this%iotas

      IF (lasym) THEN
         rmns = this%rmns
         zmnc = this%zmnc
         lmnc = this%lmnc
         gmns = this%gmns
         bsubumns = this%bsubumns
         bsubvmns = this%bsubvmns
         bsubsmnc = this%bsubsmnc
         bsupumns = this%bsupumns
         bsupvmns = this%bsupvmns
      END IF

!  Restore currents and RZL arrays. Compute_Currents allocates the arrays
!  without checking if they are already allocated.
      IF (ALLOCATED(currumnc)) THEN
         DEALLOCATE(currumnc)
      END IF
      IF (ALLOCATED(currvmnc)) THEN
         DEALLOCATE(currvmnc)
      END IF
      IF (ALLOCATED(currumns)) THEN
         DEALLOCATE(currumns)
      END IF
      IF (ALLOCATED(currvmns)) THEN
         DEALLOCATE(currvmns)
      END IF
      CALL Compute_Currents(status)
!  The rzl_local array will not be updated if the array is allocated.
      DEALLOCATE(rzl_local)
      CALL LoadRZL

      END SUBROUTINE

      END MODULE
