!*******************************************************************************
!  File commandline_parser.f
!  Contains module commandline_parser
!  Defines the base class of the type commandline_parser_class:
!
!  Super Class: None
!  Sub Classes: None
!
!*******************************************************************************

      MODULE commandline_parser
      USE stel_kinds
      USE data_parameters

      IMPLICIT NONE

!*******************************************************************************
!  commandline parser module parameters
!*******************************************************************************

      INTEGER, PARAMETER, PRIVATE :: max_arg_length = 8
      INTEGER, PARAMETER, PRIVATE :: max_length = data_filename_length         &
     &                                          + max_arg_length + 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) commandline parser base class
!
!*******************************************************************************
      TYPE commandline_parser_class
         CHARACTER (len=max_arg_length), DIMENSION(:), POINTER ::              &
     &      arg => null()
         CHARACTER (len=data_filename_length), DIMENSION(:), POINTER ::        &
     &      value => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION commandline_parser_construct()

      IMPLICIT NONE

!  Declare Arguments
      TYPE (commandline_parser_class), POINTER ::                              &
     &   commandline_parser_construct

!  local variables
      CHARACTER (len=max_length) :: temp
      INTEGER                    :: num_args = 0, i, value_index = 0

!  Start of executable code
      temp = ''
      ALLOCATE(commandline_parser_construct)

!  Read the zeroith arg to get the number of arguments.
      CALL getcarg(0, temp, num_args)

      IF (num_args .le. 0) THEN
         CALL commandline_parser_print_help
      END IF

!  Allocate the arrays and
      ALLOCATE(commandline_parser_construct%arg(num_args))
      ALLOCATE(commandline_parser_construct%value(num_args))

!  Loop through the command line arguments, and setup the argument and value
!  arrays
      DO i = 1, num_args
         CALL getcarg(i, temp, num_args)

!  Check for a - as the first character.
         IF (temp(1:1) .eq. '-') THEN
            value_index = INDEX(temp, '=')
            IF (value_index .eq. 0) THEN
!  Check for help command.
               IF (TRIM(temp) .eq. '-h') THEN
                  CALL commandline_parser_print_help
               END IF

               commandline_parser_construct%arg(i) = TRIM(temp)
               commandline_parser_construct%value(i) = ''
            ELSE
               commandline_parser_construct%arg(i) =                           &
     &            temp(1:value_index - 1)
               commandline_parser_construct%value(i) =                         &
     &            temp(value_index + 1:LEN_TRIM(temp))
            END IF
         ELSE
!  Implicity set the argument to -file
            commandline_parser_construct%arg(i) = '-file'
            commandline_parser_construct%value(i) =                            &
     &         temp(1:data_filename_length)
         END IF
      END DO

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE commandline_parser_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (commandline_parser_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%arg)) THEN
         DEALLOCATE(this%arg)
         this%arg => null()
      END IF

      IF (ASSOCIATED(this%value)) THEN
         DEALLOCATE(this%value)
         this%value => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get a string value for an argument
!-------------------------------------------------------------------------------
      FUNCTION commandline_parser_get_string(this, arg)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_filename_length) ::                                   &
     &   commandline_parser_get_string
      TYPE (commandline_parser_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)  :: arg

!  Local arguments
      INTEGER i

!  Start of executable code
!  Loop through the arguments until the correct arg is found.
      commandline_parser_get_string = ''

      IF (ASSOCIATED(this%arg)) THEN
         DO i = 1, SIZE(this%arg)
            IF (TRIM(this%arg(i)) .eq. TRIM(arg)) THEN
               commandline_parser_get_string = this%value(i)
               RETURN
            END IF
         END DO
      END IF

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Query if flag is set.
!-------------------------------------------------------------------------------
      FUNCTION commandline_parser_is_flag_set(this, arg)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: commandline_parser_is_flag_set
      TYPE (commandline_parser_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)  :: arg

!  Local arguments
      INTEGER i

!  Start of executable code
!  Loop through the arguments until the correct arg is found.
      commandline_parser_is_flag_set = .false.

      IF (ASSOCIATED(this%arg)) THEN
         DO i = 1, SIZE(this%arg)
            IF (TRIM(this%arg(i)) .eq. TRIM(arg)) THEN
               commandline_parser_is_flag_set = .true.
               RETURN
            END IF
         END DO
      END IF

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Print out the help message. When help is encountered, exit the program.
!-------------------------------------------------------------------------------
      SUBROUTINE commandline_parser_print_help

      IMPLICIT NONE

!  Start of executable code
!  All command line messages need to fit within this width.
!                '           s              c                          '
      WRITE(*,*) '                                                     '
      WRITE(*,*) '                        V3FIT                        '
      WRITE(*,*) '                                                     '
      WRITE(*,*) 'Usage: xv3fit [-arg][=option] ...                    '
      WRITE(*,*) '                                                     '
      WRITE(*,*) 'Options:                                             '
      WRITE(*,*) 'All options are displayed as [arg][takesoption][text]'
      WRITE(*,*) '  -h       N Display this information                '
      WRITE(*,*) '  -d       N Use default namelist input file         '
      WRITE(*,*) '  -file    Y Specify the v3fit namelist input file   '
      WRITE(*,*) '                                                     '
      WRITE(*,*) 'Notes:                                               '
      WRITE(*,*) '  If no -arg is found -file is implicitly implied.   '
      WRITE(*,*) 'The default input filename is v3fit.in.              '
      WRITE(*,*) '                                                     '

      STOP

      END SUBROUTINE

      END MODULE
