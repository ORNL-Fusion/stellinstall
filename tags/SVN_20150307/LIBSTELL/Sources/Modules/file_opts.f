!*******************************************************************************
!>  @file file_opts.f
!>   @brief Contains module @ref file_opts
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Contains cross platform routines for manipulating files on the file system.
!>  Defines a functions to move, copy and delete a file that is cross platform.
!>  While, vmec_system contains a crossplatform way to call a command line
!>  commands, the command used will vary by platform.
!*******************************************************************************

      MODULE file_opts

      IMPLICIT NONE

      CONTAINS

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

!-------------------------------------------------------------------------------
!>  @brief Queries if the a path is absoulte.
!>
!>  Determines if the root of a path exists or if the path is relative. If the
!>  path is absoulte return true. If the path is relative return false. An empty
!>  path is treated as absoulte.
!>
!>  @param[in] path Path to test.
!>  @returns True if the path is absolute.
!-------------------------------------------------------------------------------
      FUNCTION is_absolute_path(path)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: is_absolute_path
      CHARACTER (len=*), INTENT(in) :: path

!  Start of executable code
!  Treat empty paths as absolute.
      IF (path .eq. '') THEN
         is_absolute_path = .true.
         RETURN
      END IF

#if defined(WIN32)
!  On windows systems, root paths start with a variey of cases. Taken from
!
!  http://msdn.microsoft.com/en-us/library/windows/desktop/aa365247(v=vs.85).aspx#paths
!
!    * A UNC name of any format, which always start with two backslash characters
!      ("\\"). For more information, see the next section.
!    * A disk designator with a backslash, for example "C:\" or "d:\".
!    * A single backslash, for example, "\directory" or "\file.txt". This is
!      also referred to as an absolute path.
!
      is_absolute_path = (path(1:2)         .eq. '\\') .or.                    &
     &                   (INDEX(path, ':\') .eq. 0)    .or.                    &
     &                   (path(1:1)         .eq. '\')
#else
!  On unix systems, an absolute path starts with either / or ~
      is_absolute_path = (path(1:1) .eq. '/') .or.                             &
     &                   (path(1:1) .eq. '~')
#endif

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Returns the directory of a file.
!>
!>  Searches a path for last path separator and returns the directory path.
!>
!>  @param[in] path Path to file.
!>  @returns The path of the file.
!-------------------------------------------------------------------------------
      FUNCTION get_path_of_file(file)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=300)           :: get_path_of_file
      CHARACTER (len=*), INTENT(in) :: file

!  Start of executable code
#if defined(WIN32)
      get_path_of_file =                                                       &
     &   file(1:(INDEX(TRIM(file), '\', back = .true.) - 1))
#else
      get_path_of_file =                                                       &
     &   file(1:(INDEX(TRIM(file), '/', back = .true.) - 1))
#endif

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Builds a path.
!>
!>  Builds a path by concatenating a string onto a path.
!>
!>  @param[in] path Existing path.
!>  @param[in] name String to add to path.
!>  @returns The new path.
!-------------------------------------------------------------------------------
      FUNCTION build_path(path, name)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=300)           :: build_path
      CHARACTER (len=*), INTENT(in) :: path
      CHARACTER (len=*), INTENT(in) :: name

!  Start of executable code
      IF (TRIM(path) .eq. '') THEN
         build_path = name
         RETURN
      END IF

#if defined(WIN32)
      build_path = TRIM(path) // '\' // name
#else
      build_path = TRIM(path) // '/' // name
#endif

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Change working directory
!>
!>  Changes the current working directory of the running process. CHDIR is an
!>  extension on many compilers. Check for compatability.
!>
!>  @param[in] path Path of the working directory to change to.
!-------------------------------------------------------------------------------
      SUBROUTINE change_directory(path, error)

#if defined(__INTEL_COMPILER)
      USE IFLPORT
#endif

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: path
      INTEGER, INTENT(out)          :: error

!  local parameters
      error = CHDIR(TRIM(path))

      END SUBROUTINE

      END MODULE
