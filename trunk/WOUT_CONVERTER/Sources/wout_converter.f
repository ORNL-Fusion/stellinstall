!*******************************************************************************
!>  @file wout_converter.f
!>  @brief Utility to convert a NetCDF based wout file to a text based wout 
!>  file.
!>  @ref vmec_equilibrium
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief wout_converter main program.
!>
!>  Converts a text based wout file to a netcdf file.
!-------------------------------------------------------------------------------
      PROGRAM wout_converter
      USE read_wout_mod
      USE v3_utilities

      IMPLICIT NONE

!  local variables
      CHARACTER (len=100) :: filename
      INTEGER             :: num_cl_args
      INTEGER             :: i, status

!  Start of executable code
      num_cl_args = iargc()

      IF (num_cl_args .lt. 1) THEN
         CALL err_fatal('Expected wout filename for argument.')
      END IF

      CALL getarg(1, filename)

      CALL read_wout_file(filename, status)
      CALL assert_eq(0, status, 'Failed to load wout file ' // filename)

      IF (filename(1:5) .eq. 'wout_') THEN
         filename = TRIM(filename(6:LEN_TRIM(filename)))
      END IF

      i = INDEX(filename, '.nc', BACK=.true.)

      IF (i .gt. 0) THEN
         filename = filename(1:i - 1)
      END IF

      CALL write_wout_text(filename, status)
      CALL assert_eq(0, status, 'Failed to write wout file')

      END PROGRAM


