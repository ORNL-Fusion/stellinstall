!*******************************************************************************
!  File boundary_mc.f
!  Contains module boundary_mc
!           SUBROUTINE boundary_model_compute
!           FUNCTION boundary_model
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine boundary_model_compute
!  Both an boundary_desc and a model are needed for this computation.
!*******************************************************************************
!  MODULE boundary_mc
!    (boundary - Outer boundary Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
! SECTION IV.    _OLD (SUPERCEDED) CODE
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE boundary_mc

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations, constants, utilities
!-------------------------------------------------------------------------------
      USE v3f_global

!-------------------------------------------------------------------------------
!  boundary Derived Type
!-------------------------------------------------------------------------------
      USE boundary_T

!------------------------------------------------------------------------------
      IMPLICIT NONE

!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS

!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Compute an boundary signal
!
!    Information comes from the boundary_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = boundary
!        Outer boundary geometric
!        signal_model_compute_boundary
!
!-------------------------------------------------------------------------------

      SUBROUTINE boundary_model_compute(a_boundary, a_model, mod_signal,       &
     &                                  mod_sigma)
      IMPLICIT NONE
!-------------------------------------------------------------------------------
! ARGUMENTS
! a_boundary    - type boundary_desc - holds boundary info
! a_model       - type model - holds eq_state and rbx and zbx coefficients
! mod_signal    - output of the generated signal
! mod_sigma     - output sigma
!-------------------------------------------------------------------------------
      TYPE (boundary_desc), INTENT (inout) :: a_boundary
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:)   :: mod_signal, mod_sigma

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------
      CHARACTER(len=*), PARAMETER   :: sub_name =                              &
     &  'boundary_model_compute'
      INTEGER                       :: istat1
      INTEGER                       :: m, n, s

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

!  Update eqstate if needed.
      IF (.not. a_model%eqstate%l_def_aux1) THEN
         CALL eq_aux_1_define(a_model%eqstate)
      END IF

      s = a_model%eqstate%fixp%ns_array(1)
      m = a_boundary%m_index
      n = a_boundary%n_index

      SELECT CASE (a_boundary%coefficientName)
         CASE ('rbc')
            mod_signal(1) = model_a%eqstate%aux1%rmnc(m*n + m + 1, s)

         CASE ('zbs')
            mod_signal(1) = model_a%eqstate%aux1%zmns(m*n + m + 1, s)

         CASE ('rbs')
            mod_signal(1) = model_a%eqstate%aux1%rmns(m*n + m + 1, s)

         CASE ('zbc')
            mod_signal(1) = model_a%eqstate%aux1%zmnc(m*n + m + 1, s)

         CASE DEFAULT
            WRITE (*,*) 'Unknown coeffient type '                              &
     &                  // a_boundary%coefficientName
            STOP
      END SELECT

! make up a sigma
!         mod_sigma(1)=0.5   Changed to zero JDH 2011-10-17
      mod_sigma(1)=0.0

      END SUBROUTINE boundary_model_compute

      END MODULE boundary_mc
