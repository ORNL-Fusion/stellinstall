!*******************************************************************************
!>  @file threed1_file_cache.f
!>  @brief Contains module @ref threed1_file_cache.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref threed1_file_cache. This caches
!>  the data from the threed1 file.
!*******************************************************************************

      MODULE threed1_file_cache
      USE stel_kinds
      USE v3f_vmec_comm

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) threed1_file_cache class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a threed1_file_cache. This contains a copy of every
!>  variable that is needed from the threed1 file to compute a modeled signal.
!-------------------------------------------------------------------------------
      TYPE threed1_file_cache_class
!>  Internal inductance.
         REAL (rprec) :: vvc_smaleli
!>  Mean elongation.
         REAL (rprec) :: vvc_kappa_p
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref threed1_file_cache_class types using
!>  @ref threed1_file_cache_construct_new or                                   &
!>  @ref threed1_file_cache_construct_copy
!-------------------------------------------------------------------------------
      INTERFACE threed1_file_cache_construct
         MODULE PROCEDURE threed1_file_cache_construct_new,                    &
     &                    threed1_file_cache_construct_copy
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a new @ref threed1_file_cache_class object.
!>
!>  Allocates memory and initializes a @ref threed1_file_cache_class object.
!>
!>  @returns A pointer to a constructed @ref threed1_file_cache_class object.
!-------------------------------------------------------------------------------
      FUNCTION threed1_file_cache_construct_new()

      IMPLICIT NONE

!  Declare Arguments
      TYPE (threed1_file_cache_class), POINTER ::                              &
     &   threed1_file_cache_construct_new

!  Start of executable code
      ALLOCATE(threed1_file_cache_construct_new)

      threed1_file_cache_construct_new%vvc_smaleli = vvc_smaleli
      threed1_file_cache_construct_new%vvc_kappa_p = vvc_kappa_p

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a new @ref threed1_file_cache_class object by copying an
!>  old one.
!>
!>  Allocates memory and initializes a @ref threed1_file_cache_class object.
!>  This newly created object is a copy of the old one.
!>
!>  @param[in] cache_object @ref threed1_file_cache_class object to copy.
!>  @returns A pointer to a constructed @ref threed1_file_cache_class object.
!-------------------------------------------------------------------------------
      FUNCTION threed1_file_cache_construct_copy(cache_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (threed1_file_cache_class), POINTER ::                              &
     &   threed1_file_cache_construct_copy
      TYPE (threed1_file_cache_class), INTENT(in) :: cache_object

!  Start of executable code
      ALLOCATE(threed1_file_cache_construct_copy)

      threed1_file_cache_construct_copy%vvc_smaleli =                          &
     &   cache_object%vvc_smaleli
      threed1_file_cache_construct_copy%vvc_kappa_p =                          &
     &   cache_object%vvc_kappa_p

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref threed1_file_cache_class object.
!>
!>  Deallocates memory and uninitializes a @ref threed1_file_cache_class object.
!>
!>  @param[inout] this A @ref threed1_file_cache_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE threed1_file_cache_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (threed1_file_cache_class), POINTER :: this

!  Start of executable code
      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Save the internal state of the threed1 file.
!>
!>  Saves the threed1 file variables.
!>
!>  @param[inout] this A @ref threed1_file_cache_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE threed1_file_cache_save_state(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (threed1_file_cache_class), INTENT(inout) :: this

!  Start of executable code
      this%vvc_smaleli = vvc_smaleli
      this%vvc_kappa_p = vvc_kappa_p

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Reset the internal state of the wout file.
!>
!>  Copies the saved wout file variables back to @ref v3f_vmec_comm.
!>
!>  @param[inout] this A @ref threed1_file_cache_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE threed1_file_cache_reset_state(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (threed1_file_cache_class), INTENT(inout) :: this

!  Start of executable code
      vvc_smaleli = this%vvc_smaleli
      vvc_kappa_p = this%vvc_kappa_p

      END SUBROUTINE

      END MODULE
