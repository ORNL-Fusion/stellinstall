!*******************************************************************************
!  File ipch_mc.f
!  Contains module ipch_mc
!           SUBROUTINE ipch_mc_model_compute
!           FUNCTION emissivity_model
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine ipch_mc_model_compute
!  Both an ipch_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE ipch_mc
!    (ipch - Interferometry-Polarimetry Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
! SECTION IV.    _OLD (SUPERCEDED) CODE
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE ipch_mc

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
!  ipch Derived Type
!-------------------------------------------------------------------------------
      USE ipch_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE pol_context
         TYPE(model) :: a_model
         REAL(rprec) :: wavelength
      END TYPE

!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************

      PRIVATE :: int_function, pol_function, pol_context

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Compute an ipch signal
!
!    Information comes from the ipch_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = ipch
!        Interferometry-Polarimetry Diagnostic
!        signal_model_compute_ipch 
!
!-------------------------------------------------------------------------------

      SUBROUTINE ipch_mc_model_compute(a_ipch,a_model,mod_signal,              &
     &   mod_sigma)
      USE stel_constants
      USE read_wout_mod
      IMPLICIT NONE

!-------------------------------------------------------------------------------
! ARGUMENTS
! a_ipch        - type ipch_desc - holds ip chord info
! a_model       - type model - holds eq_state and density model information
! mod_signal    - output of the generated signal
! mod_sigma     - output sigma      
!------------------------------------------------------------------------------- 
      TYPE (ipch_desc), INTENT (inout)     :: a_ipch
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:) :: mod_signal, mod_sigma

!------------------------------------------------------------------------------- 
! local variables
!------------------------------------------------------------------------------- 
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipch_mc_model_compute'
      INTEGER, PARAMETER :: ipoints = 200 ! number of integration points
      INTEGER            :: ierr
      CHARACTER(len=1), ALLOCATABLE :: context(:)
      TYPE(pol_context)  :: p_context
      INTEGER            :: context_length

!-------------------------------------------------------------------------------
! cyl2flx declarations  - used for the subroutine cyl2flx call
!------------------------------------------------------------------------------- 
      REAL(rprec), DIMENSION(3) :: r_cyl
      REAL(rprec), DIMENSION(3) :: xcart, dxcart
      REAL(rprec), DIMENSION(3) :: bcart = 0

!------------------------------------------------------------------------------- 
! other declarations
!
!-------------------------------------------------------------------------------        
      INTEGER       :: istat1
      REAL(rprec)   :: dx               ! integration step size
      INTEGER       :: istep

!-------------------------------------------------------------------------------
!  START OF EXECUTABLE CODE
!-------------------------------------------------------------------------------         
!  Allocate data and sigma
      IF (ASSOCIATED(mod_signal)) THEN
         DEALLOCATE(mod_signal,mod_sigma, stat=istat1)
         CALL assert_eq(0,istat1,sub_name //                                   &
     &                  'mod_signal, sigma dealloc')
      ENDIF
      ALLOCATE(mod_signal(1),mod_sigma(1),stat=istat1)
      CALL assert_eq(0,istat1,sub_name // 'mod_signal, sigma alloc')

!  Load up the wout file
      IF (.not. a_model%eqstate%l_def_aux1) THEN
         CALL eq_aux_1_define(a_model%eqstate)
      END IF

      SELECT CASE(a_ipch % ip_type)
            
         CASE DEFAULT
!  Compute the line average density
            context_length = SIZE(TRANSFER(a_model, context))
            ALLOCATE(context(context_length))
            context = TRANSFER(a_model, context)
            mod_signal(1) = path_integrate(a_ipch%chordPath,                   &
     &                                     int_function,                       &
     &                                     context)
            DEALLOCATE(context)

         CASE ('p')
!  Compute the faraday rotation angle in radians.
            p_context%a_model = a_model
            p_context%wavelength = a_ipch%wavelength
            context_length = SIZE(TRANSFER(p_context, context))
            ALLOCATE(context(context_length))
            context = TRANSFER(p_context, context)
            mod_signal(1) = path_integrate(a_ipch%chordPath,                   &
     &                                     pol_function,                       &
     &                                     context)
            DEALLOCATE(context)

      END SELECT

      IF (a_ipch%ip_type .eq. 'p' .and. a_ipch%inDegrees) THEN
         mod_signal(1) = mod_signal(1)/degree
      END IF

! make up a sigma         
!         mod_sigma(1)=0.5   Changed to zero JDH 2011-10-17
      mod_sigma(1)=0.0
         
      END SUBROUTINE ipch_mc_model_compute

!*******************************************************************************
!  Interferometer Integration Function
!
!  lineAve = Int[ne*dl]
!
!*******************************************************************************
      FUNCTION int_function(context, xcart, dxcart, dx)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: int_function

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
      int_function = model_get_ne_xcyl(a_model, r_cyl)*dx

      END FUNCTION

!*******************************************************************************
!  Polarimetry Integration Function
!
!  Phi = 2.62E-13*lamba^2*Int[ne*B.dl]
!
!*******************************************************************************
      FUNCTION pol_function(context, xcart, dxcart, dx)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: pol_function

!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      TYPE(pol_context) :: p_context
      REAL(rprec), DIMENSION(3) :: r_cyl, bcart

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      p_context = TRANSFER(context, p_context)
      r_cyl = cart_to_cyl(xcart)
      CALL eq_get_B_field(p_context%a_model%eqstate, r_cyl, bcart)
      pol_function = 2.62E-13*(p_context%wavelength**2)                        &
     &             * model_get_ne_xcyl(p_context%a_model, r_cyl)               &
     &             * DOT_PRODUCT(bcart, dxcart)

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
!
!  JDH 2012-03-15
!    First version for ip - based on sxr

      END MODULE ipch_mc
