!*******************************************************************************
!  File v3fit_context.f
!  Contains module v3fit_context
!  Defines the base class of the type v3fit_context_class:
!
!  Super Class: None
!  Sub Classes: None
!
!*******************************************************************************

      MODULE v3fit_context
      USE v3fit_input
      USE model
      USE commandline_parser

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) v3fit context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  This object holds all the memory for running v3fit.
!-------------------------------------------------------------------------------
      TYPE v3fit_context_class
         TYPE (commandline_parser_class), POINTER :: cl_parser => null()
         TYPE (model_class), POINTER              :: model => null()
         TYPE (signal_pointer), DIMENSION(:), POINTER ::                       &
     &      signals => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION v3fit_context_construct(cl_parser)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), POINTER :: v3fit_context_construct
      TYPE (commandline_parser_class), POINTER :: cl_parser

!  Start of executable code
      ALLOCATE(v3fit_context_construct)

      v3fit_context_construct%cl_parser => cl_parser

!  Allocate the number of of signals needed. This is allocated to the size
!  na_s_desc. na_s_desc will default to the maximum number of signals needed.
      ALLOCATE(v3fit_context_construct%signals(na_s_desc))

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE v3fit_context_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), POINTER, INTENT(inout) :: this

!  local variables
      INTEGER                                            :: i

!  Start of executable code
      IF (ASSOCIATED(this%cl_parser)) THEN
         CALL commandline_parser_destruct(this%cl_parser)
         this%cl_parser => null()
      ENDIF

      IF (ASSOCIATED(this%model)) THEN
         CALL model_destruct(this%model)
         this%model => null()
      ENDIF

!  Deconstruct and deallocate all the signals.
      IF (ASSOCIATED(this%signals)) THEN
         DO i = 1, SIZE(this%signals)
            IF (ASSOCIATED(this%signals(i)%p)) THEN
               CALL signal_destruct(this%signals(i)%p)
               this%signals(i)%p => null()
            ENDIF
         END DO
         DEALLOCATE(this%signals)
         this%signals => null()
      ENDIF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Reduce memory useage by resizing the arrays to the smallest necessary size.
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_resize(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout)    :: this

!  local variables
      INTEGER                                      :: i
      INTEGER                                      :: minsize = 0
      TYPE (signal_pointer), DIMENSION(:), POINTER :: temp_signal

!  Start of executable code
      WRITE (*,*) ' *** Resizing context arrays'

!  Search an array for the first instance of a null pointer.
      IF (ASSOCIATED(this%signals)) THEN
         DO i = 1, SIZE(this%signals)
            IF (.not.ASSOCIATED(this%signals(i)%p)) THEN
               minsize = i - 1
               EXIT
            ENDIF
         END DO

!  Reallocate the array if it is smaller
         IF (minsize .lt. SIZE(this%signals)) THEN
            ALLOCATE(temp_signal(minsize))
            DO i = 1, minsize
                temp_signal(i)%p => this%signals(i)%p
            END DO

            DEALLOCATE(this%signals)
            this%signals => temp_signal
         ENDIF
      ENDIF

      END SUBROUTINE

      END MODULE
