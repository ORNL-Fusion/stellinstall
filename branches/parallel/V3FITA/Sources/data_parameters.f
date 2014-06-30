!*******************************************************************************
!>  @file data_parameters.f
!>  @brief Contains module @ref data_parameters.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This modules contains parameters used by equilibrium models.
!*******************************************************************************
      MODULE data_parameters

      IMPLICIT NONE

!*******************************************************************************
!  signal module parameters
!*******************************************************************************
!>  Maximum length of short data strings.
      INTEGER, PARAMETER           :: data_short_name_length = 30
!>  Maximum length of data strings.
      INTEGER, PARAMETER           :: data_name_length = 80
!>  Maximum length of file name. This needs to be large to accommodate the
!>  path to the file as well.
      INTEGER, PARAMETER           :: data_filename_length = 300

!>  Default parameter id specifiying no id.
!>  @see model
!>  @see equilibrium
      INTEGER, PARAMETER           :: data_no_id = -1

!>  Max number of parameter indicies.
      INTEGER, PARAMETER           :: data_max_indicies = 2

#if defined(MPI_OPT)
!>  Directory prefix name.
      CHARACTER (len=*), PARAMETER :: dir_prefix = 'process_dir'
!>  Directory string length.
      INTEGER, PARAMETER           :: dir_prefix_len = 14
#endif

      CONTAINS
!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Generates a base filename.
!>
!>  Parses the filename to generate a base name. Looks for a v3fit. and trims
!>  off the remainder.
!>
!>  @param[in] namelist_file The file name of the input file.
!-------------------------------------------------------------------------------
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
