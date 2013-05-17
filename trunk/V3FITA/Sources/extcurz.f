!*******************************************************************************
!  File magnetic.f
!  Contains module extcurz
!  Defines the base class of the type extcurz:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE extcurz

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) extcurz base class
!
!*******************************************************************************

      TYPE extcurz_class
         REAL(rprec)                        :: s0
         REAL(rprec), DIMENSION(:), POINTER :: svec
      END TYPE extcurz_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION extcurz_construct(s0, eqstate)
      USE eq_T

      IMPLICIT NONE

!  Declare Arguments
      TYPE (extcurz_class)           :: extcurz_construct
      TYPE (eq_state), INTENT(inout) :: eqstate
      REAL(rprec), INTENT(in)        :: s0

!  local variables
      INTEGER                        :: ns, status, i
      REAL(rprec)                    :: ds

!  Start of executable code
      extcurz_construct%s0 = s0

      IF (.not. eqstate%l_def_aux1) THEN
         CALL eq_aux_1_define(eqstate)
      END IF

      ns = eqstate%aux1%ns - 1
      ds = 1.0/ns

      ALLOCATE(extcurz_construct%svec(ns), STAT=status)

      DO i = 1, ns
         extcurz_construct%svec(i) = i*ds + ds/2.0
      END DO

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE extcurz_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (extcurz_class), INTENT(inout) :: this

!  Start of executable code
      this%s0 = 0.0

      DEALLOCATE(this%svec)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION extcurz_get_modeled_signal(this, a_model)
      USE model_T
      USE line_segment, only: line_seg
      USE stel_constants, only: twopi, mu0

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)       :: extcurz_get_modeled_signal
      TYPE (extcurz_class), INTENT(in) :: this
      TYPE (model), INTENT(inout)      :: a_model

!  local variables
      REAL(rprec)                      :: bsubv00c
      INTEGER                          :: u_bound

!  Start of executable code
      IF (.not. a_model%eqstate%l_def_aux1) THEN
         CALL eq_aux_1_define(a_model%eqstate)
      END IF

      u_bound = UBOUND(this%svec,1)
      CALL line_seg(this%s0, bsubv00c, this%svec,                              &
     &              a_model%eqstate%aux1%bsubvmnc(1,2:u_bound + 1),            &
     &              u_bound)

      extcurz_get_modeled_signal(1) = twopi*bsubv00c/mu0

      END FUNCTION

      END MODULE
