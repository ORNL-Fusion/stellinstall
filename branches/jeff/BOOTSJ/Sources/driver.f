      PROGRAM driver
      USE stel_kinds
      USE parambs, ONLY : lscreen
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat, numargs, iunit_in = 10
      CHARACTER*120 :: extension
      CHARACTER*50 :: arg1, arg2
      REAL(rprec) :: curtor
C-----------------------------------------------
!
!     driver: reads from command line file the wout file extension and surfaces
!     (half-radial) at which the bootstrap current is  required
!     writes the bootstrap current, jdotb to a file, jbbs.extension
!
!     call this as follows:   xbootjs input.boots (t or f)
!
!     where input.boz contains the wout file extension and the jrad values (as a
!     blank-delimited list, not necessarily on a single line):
!
!     FILE_EXTENSION
!     2  3   5   10  12
!
!     The surface numbers are relative to half-mesh vmec quantities.
!     thus, 2 is the first half mesh point.
!
!     The optional (T) or (F) argument allows (default) or suppresses output to screen.
!
      lscreen = .true.

      CALL getcarg(1, arg1, numargs)
      IF (numargs .gt. 1) CALL getcarg(2, arg2, numargs)

      IF (numargs .lt. 1 .or.
     1   (arg1 .eq. '-h' .or. arg1 .eq. '/h')) THEN
         PRINT *,
     1   ' ENTER INPUT FILE NAME ON COMMAND LINE'
         PRINT *,' For example: xbootsj in_bootsj.tftr'
         PRINT *
         PRINT *,
     1   ' Optional command line argument to suppress screen output:'
         PRINT *,' xbootsj input_file <(T or F)>'
         PRINT *
         PRINT *,' WHERE F will suppress screen output'
         STOP
      ELSE IF (numargs .gt. 1) THEN
         IF (arg2(1:1).eq.'f' .or. arg2(1:1).eq.'F') lscreen = .false.
      ENDIF


!      INPUT FILE
!      1st line:   extension WOUT file
!      2nd line:   surfaces to be computed


      CALL safe_open(iunit_in, istat, TRIM(arg1), 'old', 'formatted')
      IF (istat .ne. 0) STOP 'error opening input file in bootsj'

      READ (iunit_in, *,iostat=istat) extension

      CALL bootsj(curtor, TRIM(extension), iunit_in)

      END PROGRAM driver
