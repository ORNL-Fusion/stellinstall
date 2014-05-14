!*******************************************************************************
!  File sxrch_mc.f
!  Contains module sxrch_mc
!           SUBROUTINE sxrch_mc_model_compute
!           FUNCTION emissivity_model
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine sxrch_mc_model_compute
!  Both an sxrch_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE sxrch_mc
!    (sxrch - Soft X-Ray Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
! SECTION IV.    _OLD (SUPERCEDED) CODE
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE sxrch_mc

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations, constants, utilities
!-------------------------------------------------------------------------------
      USE stel_kinds
      USE v3f_global
      USE vmec_utils
      USE read_wout_mod
      USE coordinate_utilities

!-------------------------------------------------------------------------------
!  sxrch Derived Type
!-------------------------------------------------------------------------------
      USE sxrch_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!------------------------------------------------------------------------------
      IMPLICIT NONE
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************

      PRIVATE :: sxr_function

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Compute an sxrch signal
!
!    Information comes from the sxrch_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = sxrch
!        Soft X-Ray Chord Diagnostic
!        signal_model_compute_sxrch 
!
!-------------------------------------------------------------------------------

      SUBROUTINE sxrch_mc_model_compute(a_sxrch,a_model,mod_signal,            &
     &                                  mod_sigma)
      IMPLICIT NONE
!-------------------------------------------------------------------------------
! ARGUMENTS
! a_sxrch        - type sxrch_desc - holds sxr chord info
! a_model       - type model - holds eq_state and density model information
! mod_signal    - output of the generated signal
! mod_sigma     - output sigma      
!------------------------------------------------------------------------------- 
      TYPE (sxrch_desc), INTENT (inout)    :: a_sxrch
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:)   :: mod_signal, mod_sigma

!------------------------------------------------------------------------------- 
! local variables
!------------------------------------------------------------------------------- 
      CHARACTER(len=*), PARAMETER   :: sub_name =                              &
     &  'sxrch_mc_model_compute'
      INTEGER                       :: istat1
      CHARACTER(len=1), ALLOCATABLE :: context(:)
      INTEGER                       :: context_length

!-------------------------------------------------------------------------------
!  START OF EXECUTABLE CODE
!-------------------------------------------------------------------------------         
!  Allocate data and sigma
      IF (ASSOCIATED(mod_signal)) THEN
         DEALLOCATE(mod_signal,mod_sigma, stat=istat1)
         CALL assert_eq(0,istat1,sub_name //                                   &
     &                     'mod_signal, sigma dealloc')
      ENDIF
      ALLOCATE(mod_signal(1),mod_sigma(1),stat=istat1)
      CALL assert_eq(0,istat1,sub_name // 'mod_signal, sigma alloc')

!  Load up the wout file
      IF (.not. a_model%eqstate%l_def_aux1) THEN
         CALL eq_aux_1_define(a_model%eqstate)
      END IF

!  Compute the line soft xray emissivity
      context_length = SIZE(TRANSFER(a_model, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(a_model, context)
      mod_signal(1) = path_integrate(a_sxrch%chordPath,                        &
     &                               sxr_function, context)
      DEALLOCATE(context)

! make up a sigma         
!         mod_sigma(1)=0.5   Changed to zero JDH 2011-10-17
      mod_sigma(1)=0.0
         
      END SUBROUTINE sxrch_mc_model_compute

!*******************************************************************************
!  SXRem Integration Function
!
!  lineAve = Int[sxrem*dl]
!
!*******************************************************************************
      FUNCTION sxr_function(context, xcart, dxcart, dx)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: sxr_function

!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      TYPE(model) :: a_model
      REAL(rprec), DIMENSION(3) :: r_cyl

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      a_model = TRANSFER(context, a_model)
      r_cyl = cart_to_cyl(xcart)
      sxr_function = model_get_sxrem_xcyl(a_model, r_cyl)*dx

      END FUNCTION

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2009-03-03. First version of module. Based on coding from previous
!    module signal_model
!
!  GJH 2010-01-20. SXR version of coding based on module mddc_mc
!
!  JDH 2011-07-28
!    Refactoring - move emissivity function to the model. Test.
!    Moved to _old - sxrch_mc_model_compute_old and emissivity_model_old
!
!  JDH 2011-08-01
!    Refactor sxrc -> sxrch
!
!  JDH 2011-09-08
!    Remove old code. Use eq_get_flux instead of cyl2flx.
!
!  JDH 2011-10-17
!     Change to cartesian coordinates in sxrch
!
!  JDH 2011-10-25
!     Some refactoring, changed s_array to sxrem_array.

      END MODULE sxrch_mc
