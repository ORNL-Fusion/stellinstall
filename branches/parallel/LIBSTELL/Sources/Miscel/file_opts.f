!*******************************************************************************
!>  @file file_opts.f
!>  @brief Contains cross platform routines for manipulating files on the file
!>  system.
!>
!>  Defines a functions to move, copy and delete a file that is cross platform.
!>  While, vmec_system contains a crossplatform way to call a command line
!>  commands, the command used will vary by platform.
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Moves the source file to the destination.
!>
!>  Moves a file by calling the appropriate command for the platform. If the
!>  destination is an existing file, that file will be overwritten.
!>
!>  @param[in]  file_source The file to move.
!>  @param[in]  file_dest   The destination of the moved file.
!>  @param[out] error       The result of the move.
!-------------------------------------------------------------------------------
      SUBROUTINE move_file(file_source, file_dest, error)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: file_source
      CHARACTER (len=*), INTENT(in) :: file_dest
      INTEGER, INTENT(out)          :: error

!  local parameters
!  To avoid any possible interactive command lines, make sure the command
!  overwrites the file if it exists.
#if defined(WIN32)
      CHARACTER (len=*), PARAMETER :: cmd = 'move /y '
#else
      CHARACTER (len=*), PARAMETER :: cmd = 'mv -f '
#endif

!  Start of executable code

      CALL vmec_system(cmd // TRIM(file_source) // ' ' // file_dest,           &
     &                 error)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Copies the source file to the destination.
!>
!>  Copies a file by calling the appropriate command for the platform. If the
!>  destination is an existing file, that file will be overwritten.
!>
!>  @param[in]  file_source The file to copy.
!>  @param[in]  file_dest   The destination of the copied file.
!>  @param[out] error       The result of the copy.
!-------------------------------------------------------------------------------
      SUBROUTINE copy_file(file_source, file_dest, error)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: file_source
      CHARACTER (len=*), INTENT(in) :: file_dest
      INTEGER, INTENT(out)          :: error

!  local parameters
!  To avoid any possible interactive command lines, make sure the command
!  overwrites the file if it exists.
#if defined(WIN32)
      CHARACTER (len=*), PARAMETER :: cmd = 'copy /y '
#else
      CHARACTER (len=*), PARAMETER :: cmd = 'cp -f '
#endif

!  Start of executable code

      CALL vmec_system(cmd // TRIM(file_source) // ' ' // file_dest,           &
     &                 error)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deletes a file
!>
!>  Deletes a file by calling the appropriate command for the platform. Ignores
!>  errors if the file doesn't exisit.
!>
!>  @param[in]  file_source The file to delete.
!>  @param[out] error       The result of the copy.
!-------------------------------------------------------------------------------
      SUBROUTINE delete_file(file_source, error)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: file_source
      INTEGER, INTENT(out)          :: error

!  local parameters
!  To avoid any possible interactive command lines, make sure the command
!  forces execution.
#if defined(WIN32)
      CHARACTER (len=*), PARAMETER :: cmd = 'del /f '
#else
      CHARACTER (len=*), PARAMETER :: cmd = 'rm -f '
#endif

!  Start of executable code

      CALL vmec_system(cmd // TRIM(file_source), error)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Makes a directory
!>
!>  Creates new directory by calling the appropriate command for the plaform.
!>
!>  @param[in]  file_source The file to delete.
!>  @param[out] error       The result of the copy.
!-------------------------------------------------------------------------------
      SUBROUTINE create_directory(directory_source, error)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: directory_source
      INTEGER, INTENT(out)          :: error

!  local parameters
#if defined(WIN32)
      CHARACTER (len=*), PARAMETER :: cmd = 'mkdir '
#else
      CHARACTER (len=*), PARAMETER :: cmd = 'mkdir -p '
#endif

!  Start of executable code

      CALL vmec_system(cmd // TRIM(directory_source), error)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deletes a directory
!>
!>  Deletes directory by calling the appropriate command for the plaform. This
!>  deletes all files inside the directory as well.
!>
!>  @param[in]  file_source The file to delete.
!>  @param[out] error       The result of the copy.
!-------------------------------------------------------------------------------
      SUBROUTINE delete_directory(directory_source, error)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: directory_source
      INTEGER, INTENT(out)          :: error

!  local parameters
#if defined(WIN32)
      CHARACTER (len=*), PARAMETER :: cmd = 'rmdir /Q /S '
#else
      CHARACTER (len=*), PARAMETER :: cmd = 'rm -rf '
#endif

!  Start of executable code

      CALL vmec_system(cmd // TRIM(directory_source), error)

      END SUBROUTINE
