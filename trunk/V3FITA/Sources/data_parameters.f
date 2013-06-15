!*******************************************************************************
!  File data_parameters.f
!
!  This modules contains parameters used by equilibrium models.
!*******************************************************************************
      MODULE data_parameters

      IMPLICIT NONE

!*******************************************************************************
!  signal module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: data_name_length = 80
      INTEGER, PARAMETER :: data_filename_length = 300

      INTEGER, PARAMETER :: data_no_id = -1

      CONTAINS
!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
      FUNCTION filename_base(filename)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_filename_length) :: filename_base
      CHARACTER (len=*), INTENT(in)        :: filename

!  Start of executable code
      IF ((filename(1:4) .eq. 'v3fi') .and.                                    &
     &    (filename(6:6) .eq. '.')) THEN
         filename_base = TRIM(filename(7:LEN_TRIM(filename)))
      ELSE
         filename_base = TRIM(filename)
      END IF

      END FUNCTION

      END MODULE
