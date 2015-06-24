!  V3RFUN is under going a complete refactoring. Allow the old code to be
!  compiled by default. The main goal of this cleanup is to allow an extention
!  D3D style plain magnetic probes.
!#undef NEW_CODE
#ifdef NEW_CODE
!-------------------------------------------------------------------------------
!  The @header2, @begin_table, @item3 and @end_table commands are custom defined
!  commands in Doxygen.in. They are defined under ALIASES. For the page created
!  here, the 80 column limit is exceeded. Arguments of aliases are separated by
!  ','. If you intended ',' to be a string you must use an escaped comma '\,'.
!
!>  @page v3rfun_cl_parsing_sec Command Line Arguments
!>
!>  @tableofcontents
!>
!>  @section v3rfun_cl_parsing_intro Introduction
!>  This contains a description of the command line arguments. All arguments
!>  take the form of
!>
!>  @fixed_width{xv3rfun [-arg] filename}
!>
!>  @section v3rfun_cl_parsing_arg_sec Command Line Arguments
!>  @header2{Argument, Takes Value, Discription}
!>  @begin_table
!>     @item3{@fixed_width{-h},     N, Displays the help text and exits the program.}
!>  @end_table
!>
!>  @section cl_pasring_prog_ref_secProgrammers Reference
!>  Reference material for the coding to implement command line parsing is found
!>  in the @ref v3rfun.f
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file v3rfun.f
!>  @brief Contains the main routines for V3RFUN.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  The V3RFUN code is part of a suite of codes, called V3FIT, designed for
!>  equilibrium reconstruction of stellarator plasmas. It takes as input
!>  descriptions of the field coils and diagnostic coils, and generates files
!>  containing the coil and plasma response functions and the plasma response
!>  functions. The output files from the V3RFUN code are designed to be used as
!>  input to the V3FIT code. Physical quantities in the V3FIT project are
!>  measured in SI units, unless otherwise specified.
!>
!>  @author James D. Hanson
!>  @author Mark Cianciosa
!>
!>  Below is a brief discription of the major top level objects of the code. For
!>  discriptions of lower level objects consult the referenced top level
!>  objects.
!*******************************************************************************
!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief V3RFUN main program.
!>
!>  Highest level V3RFUN routine.
!-------------------------------------------------------------------------------
      PROGRAM v3rfun
      USE v3rfun_context
      USE profiler
      USE v3_utilities
!$    USE omp_lib

      IMPLICIT NONE

!  local variables
!  Forward declare the v3rfun_parse_command_line function.
      INTERFACE
         FUNCTION v3rfun_parse_command_line()
         USE v3rfun_input
         IMPLICIT NONE
         CHARACTER (len=v3rfun_file_length) :: v3rfun_parse_command_line
         END FUNCTION
      END INTERFACE

      TYPE (v3rfun_context_class), POINTER :: context => null()
      TYPE (diagnostic_dot_coil), POINTER  :: head
      INTEGER                              :: id_num
#if defined(MPI_OPT)
      INTEGER                              :: error

!  Start of executable code
!  The Intel version of the MPI libraries does not provide a correct value for
!  time until after MPI_INIT is called. Make sure this is the first think called
!  so that correct timing information can be used.
      CALL MPI_INIT(error)
#endif

      CALL profiler_construct

!$    CALL OMP_SET_NUM_THREADS(OMP_GET_NUM_PROCS())
!$OMP PARALLEL
!$    error = OMP_GET_MAX_THREADS()
!$OMP END PARALLEL
!$    WRITE (*,*) error

      context =>                                                               &
     &   v3rfun_context_construct(TRIM(v3rfun_parse_command_line()))

!  Compute the coil response functions.
      head => context%coils
      id_num = 1
      DO WHILE(ASSOCIATED(head))
         SELECT CASE (head%d_type)

            CASE ('flux_loop', 'flux_loop_circular', 'magnetic_probe',         &
     &            'magnetic_probe_tokamak', 'b_rogowski', 'i_rogowski',        &
     &            'f_rogowski')
               CALL v3rfun_context_write_mrf(context, head, id_num)

            CASE ('b_point_probe')
               CALL v3rfun_context_write_point(context, head, id_num)

            CASE DEFAULT
               CALL err_fatal('Unknown magnetic diagnostic type.')

         END SELECT

         head => head%next
         id_num = id_num + 1
      END DO

      CALL v3rfun_context_destruct(context)

      CALL profiler_destruct

#if defined(MPI_OPT)
      CALL MPI_FINALIZE(error)
#endif

      END PROGRAM v3rfun

!-------------------------------------------------------------------------------
!>  @brief Parse command line arguments.
!>
!>  The command line arguments are simple. Expect either the -h flag or read the
!>  namelist input filename. If any other argument is encountered display the
!>  help and terminate the program.
!>
!>  @returns The namelist input filename.
!-------------------------------------------------------------------------------
      FUNCTION v3rfun_parse_command_line()
      USE v3rfun_input

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=v3rfun_file_length) :: v3rfun_parse_command_line

!  local variables
      INTEGER                            :: num_args, i
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Read the zeroith arg to get the number of arguments. This should also be the
!  command name.
      CALL getcarg(0, v3rfun_parse_command_line, num_args)

      IF (num_args .le. 0 .or. num_args .gt. 2) THEN
         CALL v3rfun_print_help
      END IF

      DO i = 1, num_args
         CALL getcarg(i, v3rfun_parse_command_line, num_args)
         IF (TRIM(v3rfun_parse_command_line) .eq. '-h') THEN
            CALL v3rfun_print_help
         END IF
      END DO

      CALL profiler_set_stop_time('v3rfun_parse_command_line',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Print out help text.
!>
!>  Command line help text should take the form of '-flag'. A Y or N indicating
!>  the flag takes a value or not followed by a space. A short message
!>  describing the flag. A reference for all the command line arguments can be
!>  found at @ref FIXME.
!-------------------------------------------------------------------------------
      SUBROUTINE v3rfun_print_help()
      USE v3rfun_context
      USE magnetic_response

      IMPLICIT NONE

!  Start of executable code
!  All command line messages need to fit within this width.
!                '           s              c                          '
      WRITE(*,*) '                                                     '
      WRITE(*,1000) code_name
      WRITE(*,1100) magnetic_response_current
      WRITE(*,*) '                                                     '
      WRITE(*,*) 'Usage: xv3rfun [-arg] filename ...                   '
      WRITE(*,*) '                                                     '
      WRITE(*,*) 'Options:                                             '
      WRITE(*,*) 'All options are displayed as [arg][takesoption][text]'
      WRITE(*,*) '  -h       N Display this information                '
      WRITE(*,*) '                                                     '

      STOP

1000  FORMAT('Code: ',a)
1100  FORMAT('Version: ',a)

      END SUBROUTINE
#else

!*******************************************************************************
!  File v3rfun.f
!*******************************************************************************
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  BUILD RESPONSE FUNCTION TABLES                **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          v3rfun generates various response function tables for   **
!**	     USE with v3fit                                              **
!**                                                                  **
!**     CALLING ARGUMENTS: This is a stand alone code                **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**********************************************************************
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This file uses the following modules:
!       stel_kinds
!       stel_constants
!       bsc
!       bsc_cdf
!       ezcdf
!       coils_dot
!       diagnostic_dot
!       mddc_T
!       v3rfun_data
!
!    The modules coils_dot and diagnostic_dot also depend on
!       safe_open_mod 
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  See Section V, at end of file
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!
!  Executing 'xv3rfun -h' will printout a help message
!
!  There are three input files read in by xv3rfun:
! 1) A namelist input file
!      The first argument on the execute line is the name of this file
!      If there is no argument on the execute line, the default name of the
!         namelist input file is rfun.in
!      For example, 'xv3rfun my_nl_file' reads namelist input from a file named
!         my_nl_file
!
! 2) A coils_dot file
!       The name of the coils_dot file is specified in the namelist input variable
!       name_coils_dot. This file contains the specification of the magnetic field
!       generating coils. The parsing of this file is done by the module coils_dot.
!
! 3) A diagnostic_dot file
!       The name of the diagnostic_dot file is specified in the namelist input variable
!       name_diagnostic_dot. This file contains the specification of the diagnostic
!       coils (also called the sensors). The parsing of this file is done by the
!       module diagnostic_dot.
!
!  There are many output files generated by xv3rfun:
! 1) a LIST file - names of the mdsig files
! 2) a mutual inductance file - mutual inductances of external coils with each
!    diagnostic
! 3) many mdisg, one for each diagnostic coil
!      

!*******************************************************************************
!  CODE V3RFUN
!    
! SECTION I.	MAIN PROGRAM
! SECTION II.	PRIMARY SUBROUTINES (CALLED FROM MAIN)
! SECTION III.	AUXILIARY SUBROUTINES
! SECTION IV.	COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

!*******************************************************************************
! SECTION I. MAIN PROGRAM
!*******************************************************************************

      PROGRAM v3rfun
      USE profiler

      IMPLICIT NONE
      

!----------------------------------------------------------------------
!-- Initialize and get inputs                                        --
!----------------------------------------------------------------------
      CALL profiler_construct

      CALL get_input

      CALL write_mdsigs

      CALL profiler_destruct

      STOP
      END PROGRAM v3rfun

!*******************************************************************************
! SECTION II. PRIMARY SUBROUTINES (CALLED FROM MAIN)
!*******************************************************************************

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE get_input()

      USE stel_kinds
      USE bsc_T
!      USE coils_dot
      USE biotsavart, zero_toss => zero
      USE diagnostic_dot
      USE v3rfun_data
      USE safe_open_mod
      USE v3_utilities
      USE vsvd0, ONLY: nigroup

      IMPLICIT NONE

!----------------------------------------------------------------------
!  Argument declaration                                   
!----------------------------------------------------------------------
 
!----------------------------------------------------------------------
!  Local Variable declaration                                   
!----------------------------------------------------------------------
!  First set of declarations:
      INTEGER :: istat
      
!  Third set of declarations: variables needed for coils_dot file read.
!  Note: n_field_periods is declared in module v3rfun_data
!      INTEGER :: n_field_periods

!  Declarations for library function getarg
!    iargc is an external integer function
      INTEGER :: iargc
      INTEGER :: num_cl_args, i
      CHARACTER(LEN=80), DIMENSION(10) :: command_line_args
      CHARACTER(LEN=80) :: input_filename
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'get_input: '

!  Declarations for variables used for coil group rotations
      INTEGER :: iextcur, icoil
      REAL(rprec), DIMENSION(3) :: cgxcent, cg_rot_xcent_use, mean_r
      REAL(rprec), DIMENSION(3) :: zero_a3 = (/zero,zero,zero/)
      REAL(rprec) :: cur_coil, cur_total
      TYPE (bsc_rs) :: this_bsc_rs      !  Derived type from bsc, for rotation and shift

!----------------------------------------------------------------------
!  Start of Executable Code                                    
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  Parse the execute line. 
!    (Set l_mdsig and input_filename)
!    Note that iargc and getarg are NOT Fortran Intrinsic functions/subroutines.
!    However, they are defined for the following compilers:
!       IBM Xl Fortran
!    Note that the equivalent Fortran 2003 intrinsics will be 
!    command_argument_count and get_command_argument
!----------------------------------------------------------------------
      num_cl_args = iargc()

      DO i = 1,MIN(num_cl_args,SIZE(command_line_args))
         CALL getarg(i,command_line_args(i))
      END DO
      
!  Choose action based on the number of arguments
!  (Will need a more sophisticated parsing algorithm if the 
!  number of command line argments gets much larger)
      SELECT CASE(num_cl_args)
      CASE(0)
         input_filename = 'rfun.in'
      
      CASE(1)
         input_filename = TRIM(ADJUSTL(command_line_args(1)))

         IF (TRIM(ADJUSTL(command_line_args(1))) .eq. '-h') THEN
            CALL help
            STOP ' -h option halts execution'
         ENDIF
         
         IF (TRIM(ADJUSTL(command_line_args(1))) .eq. '-mdsig') THEN
            WRITE(*,*) 'No longer need -mdsig option'
            input_filename = 'rfun.in'
         ENDIF

      CASE(2)
         IF (TRIM(ADJUSTL(command_line_args(1))) .eq. '-mdsig') THEN
            WRITE(*,*) 'No longer need -mdsig option'
            input_filename = TRIM(ADJUSTL(command_line_args(2)))
         ELSE IF 
     &      (TRIM(ADJUSTL(command_line_args(2))) .eq. '-mdsig') THEN
               WRITE(*,*) 'No longer need -mdsig option'
               input_filename = TRIM(ADJUSTL(command_line_args(1)))
         ELSE
            CALL err_fatal(sub_name // ' 2 cl args, expect -mdsig',            &
     &         int=num_cl_args)
         END IF

      CASE DEFAULT
         CALL err_fatal(sub_name // ' Expected 0, 1 or 2 cl arguments',        &
     &     int=num_cl_args)

      END SELECT
      
!----------------------------------------------------------------------
!-- Initialization Values                                          --
!----------------------------------------------------------------------
      cg_shift_1 = zero
      cg_shift_2 = zero
      cg_rot_xcent = zero
      cg_rot_theta = zero
      cg_rot_phi = zero
      cg_rot_angle = zero
      l_rot_coil_center = .true.

!----------------------------------------------------------------------
!-- Namelist Input File                                          --
!----------------------------------------------------------------------
      OPEN (unit=iou_nli, status='old', file=TRIM(input_filename), 
     1      iostat=istat)
      IF (istat .ne. 0) THEN
         WRITE (6, *) ' Error opening input file: ' // input_filename
         STOP 
      END IF
      name_coils_dot = ''
      name_diagnostic_dot = ''

!      READ (iou_nli, nml=v3r_coils, iostat=istat)
!      IF (istat .ne. 0) STOP 'Error reading input namelist'
      
      READ (iou_nli, nml=v3r_coils)

      CLOSE (iou_nli)
      IF (myid .eq. 'none') myid=' '
      WRITE (*,5000) rfvers

!  Check and modify Namelist Input values
       kp_store = kp
      IF (lstell_sym) THEN
         zmax = max(abs(zmax),abs(zmin))
         zmin = - zmax
         IF (mod(kp,2) .ne. 0) THEN
            WRITE(*,*) ' lstell_sym TRUE and kp is odd, = ',kp
            WRITE(*,*) ' Quiting, in subroutine get_input'
            STOP
         ENDIF
         kp_store = kp / 2 + 1
      ENDIF

      kp_shell_store = kp_shell
      IF (lstell_sym) THEN
         kp_shell_store = kp_shell/2 + 1
      ENDIF

!----------------------------------------------------------------------
!-- Text Output File                                       --
!----------------------------------------------------------------------
      OPEN(unit=iou_out,status='unknown',file='rfun.out'//'_'//idrfun)
      WRITE (iou_out,5000) rfvers
      WRITE (iou_out, nml=v3r_coils)

!----------------------------------------------------------------------
!  Coils_dot read
!----------------------------------------------------------------------
!  l_read_coils_dot - logical in namelist input.
      IF (l_read_coils_dot) THEN
         IF (name_coils_dot .eq. '') THEN
            WRITE(*,*) 'l_read_coils_dot is TRUE, but'
            WRITE(*,*) 'name_coils_dot = blank. '
            WRITE(*,*) 'Setting l_read_coils_dot to FALSE'
            l_read_coils_dot = .false.
         END IF
      END IF
      
      IF (l_read_coils_dot) THEN
!  JDH 2011-09-09. Use module biotsavart to parse coils file
!    Module variables from biotsavart:
!    nfp_bs        Integer, number of field periods
!    coil_group    Array of Type bsc_coilcoll - contains the coils
         CALL parse_coils_file(name_coils_dot)
         n_field_periods_cd = nfp_bs
         n_field_cg = SIZE(coil_group)
         WRITE(6,*) 'coils file read, n_field_cg = ', n_field_cg
      ENDIF
!----------------------------------------------------------------------
!  Rotate and shift of coil groups
!----------------------------------------------------------------------
!  Coils are stored in array of bsc_coilcoll named coil_group, 
!  declared in module biotsavart
!  Loop over coil groups

      IF (n_field_cg .gt. nigroup) THEN
         CALL err_fatal('Number of coils groups exceeds max size ' //          &
     &                  'of coil groups. Increase size of nigroup ' //         &
     &                  'in LIBSTELL/Sources/Modules/vsvd0.f')
      END IF

      WRITE(*,*)
      WRITE(*,*) ' Rotate and Shift of the Coil Groups'
      DO iextcur = 1,n_field_cg
         WRITE(*,*)
         WRITE(*,*) ' Coil Group ', iextcur,' with s_name ',
     &      coil_group(iextcur) % s_name

!  Debug/test - spill the first coil in the coil group. Comment out if unneeded.
!         CALL bsc_spill_coil(coil_group(iextcur) % coils(1),                   &
!     &      'coil(1) before rs:')
!         WRITE(*,*)

!    Compute current-averaged center of coil group (cgxcent)
         cgxcent(1:3) = zero
         cur_total = zero
         DO icoil = 1, coil_group(iextcur) % ncoil
            cur_coil = coil_group(iextcur) % coils(icoil) % current
            cur_total = cur_total + cur_coil
            CALL bsc_mean_r(coil_group(iextcur) % coils(icoil),                &
     &         mean_r)
            cgxcent(1:3) = cgxcent(1:3) + mean_r(1:3) * cur_coil
         END DO
         IF (cur_total .ne. 0) cgxcent(1:3) = cgxcent(1:3) / cur_total
         IF (l_rot_coil_center(iextcur)) THEN
            cg_rot_xcent_use = cgxcent
         ELSE
            cg_rot_xcent_use = cg_rot_xcent(iextcur,1:3)
         ENDIF
 
!    Generate bsc_rs for first shift, and apply it
         CALL bsc_construct_rs(this_bsc_rs,zero,zero,zero,zero_a3,             &
     &      cg_shift_1(iextcur,1:3))
         CALL bsc_rot_shift(coil_group(iextcur),this_bsc_rs)

!    Generate bsc_rs for rotation and second shift, and apply it
         CALL bsc_construct_rs(this_bsc_rs,cg_rot_theta(iextcur),              &
     &      cg_rot_phi(iextcur),cg_rot_angle(iextcur),                         &
     &      cg_rot_xcent_use(1:3),cg_shift_2(iextcur,1:3))
         CALL bsc_rot_shift(coil_group(iextcur),this_bsc_rs)


         WRITE(*,1000) '   Current-Averaged center of cg = ',                  &
     &      cgxcent(1:3)
         WRITE(*,1000) '   First Shift = ', cg_shift_1(iextcur,1:3)
         WRITE(*,1000) '   Center of Rotation Used =  ',                       &
     &      cg_rot_xcent_use
         WRITE(*,1000) '   Rotation theta, phi, angle = ',                     &
     &      cg_rot_theta(iextcur), cg_rot_phi(iextcur),                        &
     &      cg_rot_angle(iextcur)
         WRITE(*,1000) '   Second Shift = ', cg_shift_2(iextcur,1:3)
1000  FORMAT(a34,3(2x,es12.5))

!  Debug/test - spill the first coil in the coil group. Comment out if unneeded.
!         CALL bsc_spill_coil(coil_group(iextcur) % coils(1),                   &
!     &      'coil(1) after rs:')

      END DO

!----------------------------------------------------------------------
!  Point field_coils (in v3rfun_data) to the module biotsavart variable
!----------------------------------------------------------------------
      IF (l_read_coils_dot) THEN
         field_coils => coil_group
      ENDIF 

!----------------------------------------------------------------------
!   Logic for number of field periods to use
!----------------------------------------------------------------------
!  JDH 2007-10-06
!   Print out a warning if both  n_field_periods_cd and  n_field_periods_nli
!   are 1. If axisymmetric system, could be a problem.
      IF ((n_field_periods_cd .eq. 1) .and.                                    &
     &    (n_field_periods_nli .eq. 1)) THEN
            WRITE(*,*) 'I noticed that both n_field_periods_cd and'
            WRITE(*,*) 'n_field_periods_nli are 1. If you are modeling'
            WRITE(*,*) 'an axisymmetric system, then you probably'
            WRITE(*,*) 'need to use n_field_periods_nli large, so that'
            WRITE(*,*) 'plasma responses for non-axisymmetric'
            WRITE(*,*) 'diagnostics are correctly averaged'
      ENDIF
         
 !  JDH 2007-08-15. Changed logic so use NLI value when disagree
 !  This allows using many field periods to approximate the
 !  plasma response appropriate for an axisymmetric plasma
      IF (l_read_coils_dot) THEN
         n_field_periods = n_field_periods_cd
         IF (n_field_periods_nli .ne. 0) THEN
            IF (n_field_periods_nli .ne. n_field_periods ) THEN
               WRITE(*,*) 'Problem. # field periods from coils_dot'
               WRITE(*,*) 'is different from # field periods from'
               WRITE(*,*) 'namelist input file. Will use _nli value'
               WRITE(*,*) '  n-_cd = ', n_field_periods_cd
               WRITE(*,*) '  n-_nli = ', n_field_periods_nli
               n_field_periods = n_field_periods_nli
            ENDIF
         ENDIF
      ELSE
         n_field_periods = n_field_periods_nli
      ENDIF
      
      IF (n_field_periods .le. 0) THEN
         WRITE(6,*) 'Problem. # field periods <= 0'
         WRITE(6,*) '  n_field_periods = ', n_field_periods
      END IF

!----------------------------------------------------------------------
!  Allocate space for the field_coils array
!  Second part of coils_dot read
!----------------------------------------------------------------------
          
      IF (l_read_coils_dot) THEN
!  JDH 2011-09-09. Use module biotsavart to parse coils file
!    Both read_coils passes were done in parse_coils_file
!    Here, point to the biotsavart module variable
         field_coils => coil_group
      ENDIF 
      
!  Read In for Diagnostic Coils
      CALL diagnostic_dot_read(name_diagnostic_dot,mddc_desc_array,            &
     &   n_diagn_c)

      WRITE(6,*) 'diagnostic. read, n_diagn_c = ', n_diagn_c

      RETURN
 5000 format (/,10x,'v3rfun Version ',a,/)
      END SUBROUTINE get_input


!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

      SUBROUTINE write_mdsigs
!  Subroutine to compute and write out the response functions
!  into an mdsig (Magnetic Diagnostic SIGnal) file.
!  Uses derived types that will be used in V3FIT
!  JDH 2007-06-11 Convert to mddc_desc
      
      USE stel_kinds
      USE stel_constants
      USE v3rfun_data
      USE ezcdf
      USE bsc_cdf
      USE diagnostic_T
      USE diagnostic_cdf
      USE mddc_T
      USE safe_open_mod
      USE v3_utilities
      USE coordinate_utilities
      USE magnetic_response

!      Variables declared in v3rfun_data used in this subroutine:
!  iou_mdsig_ncdf        i/o unit number for mdsig file
!  iou_list			     i/o unit number for LIST file
!  iou_mi                i/o unit number for the Mutual Inductance file
!  rfvers                version number                            
!  name_coils_dot        name of the 'coils.' file
!  rmin                  Minimum R for plasma grid
!  rmax                  Maximum R for plasma grid
!  zmin                  Minimum z for plasma grid
!  zmax                  Maximum z for plasma grid
!  n_field_periods       number of field periods
!  lstell_sym            Logical - True for stellarator symmetry
!  kp                    number of phi planes per field period in plasma grid
!  use_con_shell         Logical - True for computing the signal responce of a
!                                  conducting shell.
!  major_radius          Major radius for a conducting shell.
!  minor_radius          Minor radius for a conducting shell.

      IMPLICIT NONE

!  local variables
      REAL(rprec), DIMENSION(:,:,:,:,:), ALLOCATABLE :: xcart_grid_e
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE       :: phi_grid_e
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::                  &
     &   pl_response
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::                  &
     &   pl_response_ss
      REAL(rprec), DIMENSION(3)                      :: xcyl
      REAL(rprec), DIMENSION(3)                      :: acart
      INTEGER :: i, j, k, l, idc 
      INTEGER :: iss, jss, kss
      REAL(rprec) :: phi0
      REAL(rprec) :: fperiod, delr, delz, delp, delu, delp_shell
      CHARACTER (len=120) :: temp_filename, name_dd_nopath
      INTEGER :: numc, maxindex(1)
      REAL(rprec), DIMENSION(:), POINTER :: extcur => null()
      REAL(rprec), DIMENSION(:), POINTER :: rdiag_coilg_1 => null()

      TYPE (magnetic_response_class), POINTER :: response => null()
      TYPE(mddc_mrf), POINTER   :: null_mrf => null()
      TYPE(diagnostic_desc)     :: local_diagnostic_desc
      REAL(rprec), DIMENSION(:,:,:), POINTER :: a_r => null()
      REAL(rprec), DIMENSION(:,:,:), POINTER :: a_f => null()
      REAL(rprec), DIMENSION(:,:,:), POINTER :: a_z => null()

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: phi_grid_shell
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: xcart_s_grid_e
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::                    &
     &   s_response, s_response_ss
      REAl(rprec), DIMENSION(:,:), POINTER :: a_s_r => null()
      REAl(rprec), DIMENSION(:,:), POINTER :: a_s_f => null()
      REAl(rprec), DIMENSION(:,:), POINTER :: a_s_z => null()
      REAL(rprec)                          :: u
      INTEGER                              :: n_u

      INTEGER                :: status
      CHARACTER(len=8)       :: my_date
      CHARACTER(len=10)      :: my_time
      CHARACTER(len=5)       :: my_zone
      CHARACTER (len=8+10+5) :: date_run

      INTERFACE
         FUNCTION cart2cyl_v(vcart, xphi)
         USE stel_kinds
         IMPLICIT NONE
         REAL (rprec), DIMENSION(3)             :: cart2cyl_v
         REAL (rprec), DIMENSION(3), INTENT(in) :: vcart
         REAL (rprec), INTENT(in)               :: xphi
         END FUNCTION
      END INTERFACE

!  xcarte_grid_e    array of cartesian coordinates on the extended
!                      computational grid
!  phi_grid_e       array of phi values on the extended computational grid 
!  pl_response      array of plasma response values
!  pl_response_ss   array of plasma response values. On smaller array, when
!                      lstell_sym = .true.
!  xcyl				array to hold a single set of cylindrical coordinates
!  acyl				array to hold a single vector potential in cyl coordinates
!  acart			array to hold a single vector potential in cart coordinates
!  i				index, runs over radial grid
!  j				index, runs over z grid
!  k				index, runs over phi grid
!  iss              index for stellarator symmetric position, radial
!  jss              index for stellarator symmetric position, z
!  kss              index for stellarator symmetric position, phi
!  l				index, runs over field periods
!  idc				index, runs over the diagnostic coils
!  istat			status variable
!  phi0				phi value on computational grid
!  fperiod			extent of a field period, in phi
!  delr				radial spacing of computational mesh
!  delz				z spacing of computational mesh
!  delp				phi spacing of computational mesh
!  mdsig_filename   file name of the mdsig file
!  list_filename    file name of the list-file
!  extcur           array of current multipliers for the field-coil groups
!  rdiag_coilg_1    array of responses of a diagnostic to field-coil groups


!  Start of executable code

!  Preliminaries for the plasma responses --------------------------------------
!  Set the grid differental elements
      fperiod = twopi/n_field_periods
      delr = (rmax - rmin)/(ir - 1)
      delz = (zmax - zmin)/(jz - 1)
      delp = fperiod/kp

      ALLOCATE(xcart_grid_e(ir,jz,kp,n_field_periods,3))
      ALLOCATE(phi_grid_e(kp,n_field_periods))
      ALLOCATE(pl_response(ir,jz,kp,3))
      IF (lstell_sym) THEN
         ALLOCATE(pl_response_ss(ir,jz,kp_store,3))
      ENDIF

!  Allocate space for conducting shell
      IF (use_con_shell) THEN
         n_u = MAX(INT(n_field_periods*kp_shell*minor_radius/                  &
     &                 major_radius), 10)
         delu = twopi/n_u
         delp_shell = fperiod/kp_shell

         ALLOCATE(xcart_s_grid_e(n_u, kp_shell, n_field_periods, 3))
         ALLOCATE(s_response(n_u,kp_shell,3))
         ALLOCATE(phi_grid_shell(kp_shell,n_field_periods))
         IF (lstell_sym) THEN
            ALLOCATE(s_response_ss(n_u,kp_shell_store,3))
         END IF
      END IF

!  Compute the coordinates, store them in xcart_grid_e. Also compute and store
!  (in phi_grid_e) the phi values
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,k,l,xcyl,phi0)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO k = 1,kp
         phi0 = (k - 1)*delp
         DO l = 1, n_field_periods
            xcyl(2) = phi0 + (l - 1)*fperiod
            phi_grid_e(k,l) = xcyl(2)
            DO i = 1,ir
               xcyl(1) = rmin + (i - 1)*delr
               DO j = 1,jz
                  xcyl(3) = zmin + (j - 1)*delz
                  xcart_grid_e(i,j,k,l,1:3) = cyl_to_cart(xcyl)
               END DO ! over j
            END DO ! over i
         END DO ! over l
      END DO ! over k
!$OMP END DO
!$OMP END PARALLEL

!  Compute the coordinates for the conducting shell. The phi values were
!  computed in the previous loop.
      IF (use_con_shell) THEN
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,k,l,xcyl,u,phi0)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO k = 1,kp_shell
            phi0 = (k - 1)*delp_shell
            DO l = 1,n_field_periods
               xcyl(2) = phi0 + (l - 1)*fperiod
               phi_grid_shell(k,l) = xcyl(2)
               DO j = 1,n_u
                  u = (j - 1)*delu
                  xcyl(1) = major_radius + minor_radius*cos(u)
                  xcyl(3) = minor_radius*sin(u)
                  xcart_s_grid_e(j,k,l,1:3) = cyl_to_cart(xcyl)
               END DO
            END DO
         END DO
!$OMP END DO
!$OMP END PARALLEL
      END IF

!  Open mdsig List File --------------------------------------------------------
!  Strip off the leading path elements from name_diagnostic_dot
      name_dd_nopath = name_diagnostic_dot(INDEX(name_diagnostic_dot,          &
     &                                           '/' , .true.) + 1:)
      IF (LEN_TRIM(myid) .eq. 0) THEN
         temp_filename = TRIM(name_dd_nopath) // '_mdsig.LIST'
      ELSE
         temp_filename = TRIM(name_dd_nopath) // '_' // TRIM(myid) //          &
     &                   '_mdsig.LIST'
      END IF
      WRITE(*,*) 'list_filename is ', temp_filename
      CALL safe_open(iou_list, status, TRIM(temp_filename), 'unknown',         &
     &               'formatted')

!  Open Mutual Inducatance File ------------------------------------------------
      IF (LEN_TRIM(myid) .eq. 0) THEN
         temp_filename = TRIM(name_dd_nopath) // '_MI'
      ELSE
         temp_filename = TRIM(name_dd_nopath) // '_' // TRIM(myid) //          &
     &                   '_MI'
      END IF
      WRITE(*,*) 'mi_filename is ', TRIM(temp_filename)
      CALL safe_open(iou_mi, status, TRIM(temp_filename), 'unknown',             &
     &               'formatted', delim_in='n')

!----------------------------------------------------------------------
!  Preliminaries for the coil responses
!----------------------------------------------------------------------
!  Allocate space for rdiag_coilg_1 array
      ALLOCATE(rdiag_coilg_1(n_field_cg))
      ALLOCATE(extcur(n_field_cg))

      DO i = 1, n_field_cg
!  Find extcur array element for computing scaled inductance matrix
!  next 3 lines taken out of makegrid code
         numc = field_coils(i)%ncoil
         maxindex = MAXLOC(ABS(field_coils(i)%coils(1:numc)%current))
         extcur(i) = field_coils(i)%coils(maxindex(1))%current
         IF (extcur(i) == zero) THEN
            extcur(i) = 1
         END IF
      END DO

      CALL date_and_time(my_date, my_time, my_zone)
      date_run = my_date // ' ' // my_time // ' ' // my_zone

!  Write mutual inductance file header.
      WRITE(iou_mi,*) ' MUTUAL INDUCTANCES computed by ', code_name
      WRITE(iou_mi,*) ' code version is ', rfvers
      WRITE(iou_mi,*) ' Date run: ', TRIM(date_run)
      WRITE(iou_mi,*) ' Field coil information from ',                         &
     &                TRIM(name_coils_dot)
      WRITE(iou_mi,*) ' Diagnostic information from ',                         &
     &                TRIM(name_diagnostic_dot)
      WRITE(iou_mi,*) 
      
!  Loop over all diagnostic coils.
      DO idc =  1,n_diagn_c
         WRITE (6,'(5x,a12,i4)') 'Diagnostic #', idc

!  Coil Response ---------------------------------------------------------------
         DO i = 1, n_field_cg    !  Loop over coil groups
            CALL bsc_fluxba(field_coils(i), mddc_desc_array(idc)%mdcoil,       &
     &                      len_integrate_ddc, rdiag_coilg_1(i))
            rdiag_coilg_1(i) = rdiag_coilg_1(i)                                &
     &                       * mddc_desc_array(idc)%flux_factor
         END DO

!  Plasma Response -------------------------------------------------------------
!  Parallelize the outer most loops. Any local variables inside the loop must
!  be private.
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,k,l,acart)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO i = 1,ir
            DO j = 1,jz
               DO k = 1,kp
                  pl_response(i,j,k,1:3) = zero
                  DO l = 1, n_field_periods
                     CALL bsc_a(mddc_desc_array(idc)%mdcoil,                   &
     &                          xcart_grid_e(i,j,k,l,1:3), acart)
                     pl_response(i,j,k,1:3) = pl_response(i,j,k,1:3) +         &
     &                  cart2cyl_v(acart, phi_grid_e(k,l))
                  END DO ! l, field periods
                  pl_response(i,j,k,1:3) = pl_response(i,j,k,1:3) *            &
     &               mddc_desc_array(idc)%flux_factor       ! Conversion Factor
               END DO ! k
            END DO ! j
         END DO ! i
!$OMP END DO
!$OMP END PARALLEL
         
!  When using stellator symmetry, only half of the field period needs to be
!  stored due to the symmetry. A vertical point will have the same current as a
!  the symmetric point below the midplane. The counter point to point j is found
!  from by counting the grid in reverse.
         IF (lstell_sym) THEN
            DO i = 1,ir
               iss = i
               DO j = 1,jz
                  jss = jz + 1 - j
                  DO k = 1,kp_store
                     IF (k .eq. 1) THEN
                        kss = 1
                     ELSE
                        kss = kp + 2 - k
                     ENDIF
                     pl_response_ss(i,j,k,1) = pl_response(i,j,k,1) -          &
     &                  pl_response(iss,jss,kss,1)
                     pl_response_ss(i,j,k,2) = pl_response(i,j,k,2) +          &
     &                  pl_response(iss,jss,kss,2)
                     pl_response_ss(i,j,k,3) = pl_response(i,j,k,3) +          &
     &                  pl_response(iss,jss,kss,3)
                  END DO ! k
               END DO ! j
            END DO ! i
         END IF ! lstell_sym

!  Conducting Shell Response ---------------------------------------------------
         IF (use_con_shell) THEN

!  Parallelize the outer most loops. Any local variables inside the loop must
!  be private.
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,k,l,acart)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
            DO j = 1, n_u
               DO k = 1, kp_shell
                  s_response(j,k,1:3) = zero
                  DO l = 1,n_field_periods
                     CALL bsc_a(mddc_desc_array(idc)%mdcoil,                   &
     &                          xcart_s_grid_e(j,k,l,1:3), acart)
                     s_response(j,k,1:3) = s_response(j,k,1:3) +               &
     &                  cart2cyl_v(acart, phi_grid_shell(k,l))
                  END DO
               END DO
            END DO
!$OMP END DO
!$OMP END PARALLEL
            s_response = s_response*mddc_desc_array(idc)%flux_factor

!  When using stellator symmetry, only half of the field period needs to be
!  stored due to the symmetry. A poloidal point will have the same current as a
!  the symmetric point below the midplane. The j=1 point corresponds to the
!  same u position on the opposite side of the field perior. For j!=0, the
!  symmetry point is found by counting indices backwards.
            IF (lstell_sym) THEN
               DO j = 1, n_u
                  IF (j .eq. 1) THEN
                     jss = 1
                  ELSE
                     jss = n_u + 2 - j
                  END IF
                  DO k = 1, kp_shell_store
                     IF (k .eq. 1) THEN
                        kss = 1
                     ELSE
                        kss = kp_shell + 2 - k
                     END IF

                     s_response_ss(j,k,1) = s_response(j,k,1)                  &
     &                                    - s_response(jss,kss,1)
                     s_response_ss(j,k,2) = s_response(j,k,2)                  &
     &                                    + s_response(jss,kss,2)
                     s_response_ss(j,k,3) = s_response(j,k,3)                  &
     &                                    + s_response(jss,kss,3)
                  END DO
               END DO
            END IF
         END IF

!  Assign the response function grids to pointers.
         IF (lstell_sym) THEN
            a_r => pl_response_ss(:,:,:,1)
            a_f => pl_response_ss(:,:,:,2)
            a_z => pl_response_ss(:,:,:,3)
            IF (use_con_shell) THEN
               a_s_r => s_response_ss(:,:,1)
               a_s_f => s_response_ss(:,:,2)
               a_s_z => s_response_ss(:,:,3)
            END IF
         ELSE
            a_r => pl_response(:,:,:,1)
            a_f => pl_response(:,:,:,2)
            a_z => pl_response(:,:,:,3)
            IF (use_con_shell) THEN
               a_s_r => s_response(:,:,1)
               a_s_f => s_response(:,:,2)
               a_s_z => s_response(:,:,3)
            END IF
         ENDIF

!  Construct a responce object.
         response => magnetic_response_construct(code_name, date_run,          &
     &                                           TRIM(name_coils_dot),         &
     &                                           rdiag_coilg_1, extcur,        &
     &                                           kp, kp_shell,                 &
     &                                           rmin, rmax, zmin, zmax,       &
     &                                           n_field_periods,              &
     &                                           lstell_sym,                   &
     &                                           a_r, a_f, a_z,                &
     &                                           a_s_r, a_s_f, a_s_z)

!  Write Coil Response ---------------------------------------------------------
!  Open up an mdsig file, for netcdf writing
         IF (LEN_TRIM(myid) .eq. 0) THEN
            temp_filename =                                                    &
     &         TRIM(ADJUSTL(mddc_desc_array(idc)%mdcoil%s_name)) //            &
     &         '_mdsig.nc'
         ELSE
            temp_filename =                                                    &
     &         TRIM(ADJUSTL(mddc_desc_array(idc)%mdcoil%s_name)) //            &
     &         '_' // TRIM(myid) // '_mdsig.nc'
         ENDIF
         status = 0
         CALL cdf_open(iou_mdsig_ncdf, TRIM(temp_filename), 'w', status)
         CALL assert_eq(0, status,                                             &
     &                  'write_mdsigs mdsig file open status')

!  Define NetCDF variables.
         CALL cdf_define(iou_mdsig_ncdf, 'diagnostic_desc_d_type',             &
     &                   mddc_desc_array(idc)%mddc_type)
         CALL cdf_define(iou_mdsig_ncdf, 'diagnostic_desc_s_name',             &
     &                   mddc_desc_array(idc)%s_name)
         CALL cdf_define(iou_mdsig_ncdf, 'diagnostic_desc_l_name',             &
     &                   mddc_desc_array(idc)%l_name)
         CALL cdf_define(iou_mdsig_ncdf, 'diagnostic_desc_units',              &
     &                   mddc_desc_array(idc)%units)
         CALL cdf_define(iou_mdsig_ncdf,                                       &
     &                   'diagnostic_desc_sigma_default',                      &
     &                   mddc_desc_array(idc)%sigma_default)
         CALL magnetic_response_define(response, iou_mdsig_ncdf)


!  Write NetCDF variables. The cdf_write subroutine switches from define mode to
!  write mode.
         CALL cdf_write(iou_mdsig_ncdf, 'diagnostic_desc_d_type',              &
     &                  mddc_desc_array(idc)%mddc_type)
         CALL cdf_write(iou_mdsig_ncdf, 'diagnostic_desc_s_name',              &
     &                  mddc_desc_array(idc)%s_name)
         CALL cdf_write(iou_mdsig_ncdf, 'diagnostic_desc_l_name',              &
     &                  mddc_desc_array(idc)%l_name)
         CALL cdf_write(iou_mdsig_ncdf, 'diagnostic_desc_units',               &
     &                  mddc_desc_array(idc)%units)
         CALL cdf_write(iou_mdsig_ncdf, 'diagnostic_desc_sigma_default',       &
     &                  mddc_desc_array(idc)%sigma_default)
         CALL magnetic_response_write(response, iou_mdsig_ncdf)

!  Close NetCDF file.
         CALL cdf_close(iou_mdsig_ncdf)

!  Write entry in list file
          WRITE(iou_list,fmt='(i4.4,x,a)') idc, TRIM(temp_filename)

!  Write to Mutual Inductance file
         WRITE(iou_mi,'(a,t20,i4)') 'Diagnostic Coil: ', idc
         WRITE(iou_mi,'(a,t20,a)') 'Short Name: ',                             &
     &      TRIM(mddc_desc_array(idc)%s_name)
         WRITE(iou_mi,'(a,t20,a)') 'MDDC Type: ',                              &
     &      TRIM(mddc_desc_array(idc)%mddc_type)
         WRITE(iou_mi,'(a,t20,a)') 'Signal units: ',                           &
     &      TRIM(mddc_desc_array(idc)%units)
         WRITE(iou_mi,*) ' i      ID          Inductance'
         DO i = 1, n_field_cg
            WRITE(iou_mi,'(i4,1x,a11,1x,es14.6)') i,                           &
     &         TRIM(field_coils(i) % s_name), rdiag_coilg_1(i)
         END DO
         WRITE(iou_mi,*)

!  Cleanup
         CALL magnetic_response_destruct(response)
      END DO ! over diagnostic coils

!  Deallocate space  ZZZZZZZZ
      DEALLOCATE(xcart_grid_e, phi_grid_e, pl_response)
      IF (lstell_sym) THEN
         DEALLOCATE(pl_response_ss)
      ENDIF
      DEALLOCATE(rdiag_coilg_1, extcur)

      IF (use_con_shell) THEN
         DEALLOCATE(xcart_s_grid_e, s_response, phi_grid_shell)
         IF (lstell_sym) THEN
            DEALLOCATE(s_response_ss)
         END IF
      END IF

!  Close up the list file and the mutual inductance file
      CLOSE (iou_list)
      CLOSE (iou_mi)

      RETURN
      END SUBROUTINE write_mdsigs

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

!*******************************************************************************
! SECTION III.	AUXILIARY SUBROUTINES
!*******************************************************************************
!  Subroutine to convert a tangent vector from Cartesian coordinates
!  to cylindrical coordinates
      FUNCTION cart2cyl_v(vcart, xphi)
      USE stel_kinds

      IMPLICIT NONE
      
!   Declare Arguments
      REAL (rprec), DIMENSION(3)             :: cart2cyl_v
      REAL (rprec), DIMENSION(3), INTENT(in) :: vcart
      REAL (rprec), INTENT(in)               :: xphi

!  local variables
      REAL (rprec)                           :: cphi
      REAL (rprec)                           :: sphi

!  Start of executable code
      cphi = cos(xphi)
      sphi = sin(xphi)
      cart2cyl_v(1) = vcart(1)*cphi + vcart(2)*sphi
      cart2cyl_v(2) = -vcart(1)*sphi + vcart(2)*cphi
      cart2cyl_v(3) = vcart(3)

      END FUNCTION

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE help
!  Subroutine to print out help information
!  Called from get_input
!  subroutine will be called when '-h' is first argument on execute line
      USE v3rfun_data

      WRITE(*,*) 'V3RFUN - Compute Response Functions'
      WRITE(*,*) 'Version (may not be accurate)', rfvers 
      WRITE(*,*)
      WRITE(*,*) ' Zero, or one command line arguments'
      WRITE(*,*) ' One argument is the NLI file name.'
      WRITE(*,*) ' If the NLI file name is absent, then the default'
      WRITE(*,*) ' NLI file name is: rfun.in'
      WRITE(*,*) ' (There is no longer a need for the -mdsig option)'
      WRITE(*,*)
      WRITE(*,*) 'An input namelist specification SAMPLE is :'
      WRITE(*,fmt='(a)')
     &'   &v3r_coils ',
     &'    idrfun  = ''ncsx_multi'' ',
     &'    name_coils_dot = ''coils.ncsx'' ',
     &'    name_diagnostic_dot = ''diagnostics.short'' ',
     &'    rmin=0.45 ',
     &'    rmax=2.35 ',
     &'    zmin=-0.95 ',
     &'    zmax=0.95 ',
     &'    ir=101 ',
     &'    jz=101 ',
     &'    kp=32 ',
     &'    lstell_sym = T ',
     &'  / '
     
      WRITE(*,*) 'END HELP'

      RETURN
      END SUBROUTINE help

!*******************************************************************************
! SECTION IV.	COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!**          22/04/02..........first created                         **
!  JDH 15.08.02 Add in BSC module for field generating coils.
!  JDH 22.10.02 Newer bsc, add bsc coilcoll diagnostic_coils
!  JDH 19.11.2002. Major rewrite. Uses bsc, coils_dot, and diagnostic_dot
!     modules. Pulled data and kind modules into the v3rfun file.
!
!  12.03.2002 JDH
!    Added Auxiliary subroutines cyl2cart_x and cart2cyl_v. First attempt at real 
!    coding for the plasma response.
!
!  12.04.2002 JDH
!    Changed to use the stel_kinds module, external to this file.
!    Also - changes to coils_dot and diagnostic_dot, so they now depend on safe_open_mod
!
!  12.05.2002 JDH
!    Incorporated changes from Ed Lazarus - Write netcdf file

!  12.10.2002 JDH
!     Write out binary file of plasma response.
!
!  JDH 12.19.02 - Should separate off the test subroutine, and write a separate test code

!  12.20.2002
!    Added lots of comments. Moved rfun_signal to test - no longer a module
!    Changed name from rfun_dat to v3rfun_data

!  01.02.2003
!    Moved coding that was in test_rfun into a separate program.
!    Left the subroutine there, but it doesn't do anything.

!  01.02.2003
!    Eliminated some extraneous variables. Revised i/o unit number variable names.
!    Completely eliminated test_rfun.

!  JDH 01.02.2003 - List of i/o stuff
!    iou - integer input/output unit number
!    safe - is the file opened with safe_open_mod or v3_utilities
!    iouname - variable name for the i/o unit number
!    filename - name of the file, * indicates more stuff
!    var - name of the variable that carries the file name
!    types are text, namelist, unformatted, netCDF 
!    subroutine - where the file is used
!  Note: iou_ variables are declared in module v3rfun_data
!
! iou 	safe 	iouname 		filename	var			in/out 	type 		subroutine
! 11	No		iou_nli			arbitrary	input_id	In		Namelist	get_input
! ?     Yes		?				arbitrary	Good		In		text		coils_dot
! ?     Yes		?				arbitrary	Good		In		text		diagnostic_dot
! 10	No		iou_out			rfun.out*	none		out		text		get_input
! 15	No		iou_crf_unf		rfcoil.dat	none		out		unformatted	coil_rfun
! 20	No		iou_crf_ncdf	crfun*.nc	rfun_filnm	out		netCDF		coil_rfun
! 30	No		iou_prf_ncdf	*prfun.nc	cdffil		out		netCDF		plasma_rfun
! 35	Yes		iou_list		*.LIST		none		out		text		plasma_rfun
! 35	Yes		iou_list		*_mdsig.LIST none		out		text		write_mdsigs
! 40	Yes		iou_mdsig_ncdf	*_mdsig.nc	none		out		netCDF		write_mdsigs

!  01.20.2003
!    Changed plasma response function so that it gets multiplied by
!    diagnostic_coils % coils(k) % raux. This does unit conversion for magnetic probes.
!    Add logicals for writing binary coil response file, and for reading a coils_dot file.
!
!  01.23.2003
!    Finished eliminating the unformatted (binary) file write of the plasma response.
!    Added separate variables for n_field_periods from NLI and Coils_dot. Logic
!    for which to use is in subroutine get_input.
!
!  03.25.2003 SPH
!    Added a printout to screen informing which diagnostic is being stored in response table
!    Added ADJUSTL to avoid leading blanks in file name (plasma_rfun)
!
!  04.10.2003 SPH
!    Added extcur and inductance array to coil_rfun
!
!  06.26.2003 JDH
!    Added coding to do the extra summation in the case where the stellarator
!    symmetry can be utilized.
!    See Stellarator Symmetry Powerpoint of 6.10.03
!    Changed from cdf_prfun_define and cdf_prfun_write to cdf_prfun_define_write
!
!  11.07.03 JDH
!    Commented out a print line.
!
!  09-13-2004 JDH
!     Added -mdsig option. (mdsig - Magnetic Diagnostic SIGnal). Write all responses
!     (both coil and plasma) into a single file, one per diagnostic. Uses newer
!     structures from modules daignsotic_T and signal_T, and uses cdf writes.
!     Changed execute line parsing a  bit. Added ddesc_array to v3rfun_data.f
!
!  12-02-2004 JDH
!    Added raux factor to calculation of rdiag_coilg in write_mdsigs. 
!    Fixed date and time.
!
!  12-13-2004 JDH
!    Modified to take into account changes to signal_T, _cdf and
!    diagnostic_T, _cdf.
!
!  03-08-2005 JDH
!    No changes to the v3rfun file, but changed diagnostic_dot, so that units
!    of rogowski diagnostics are T [Integral (B dl) / Integral(dl)]
!    TO DO. Eliminate use of raux. Use flux_factor instead.
!
! 12-14-2005 JMS
!    added modules for implementation of interferometry and polarimetry and added a 
!    logical toggle WANT_F_ROTATION in MAIN() to allow the option to call v3rfun_intpol_io_driver
!    (which enables the writing out of B/I info to output .nc files) 
!
!  2007-06-28 JDH
!    Fixed test for kp odd, only needed when lstell_sym is TRUE.
!

!  2007-10-06 JDH
!    Earlier changes, and today. Fix problem with number of field periods,
!    when trying to model an axisymmetric system.
!
!  2007-10-07 JDH
!   Messed up earlier commit. Also, took out CVS keyword expansion stuff at
!   beginning of file.
!
!  2008-05-16 JDH
!    Eliminated (iprec). Fixed location of LIST file. Changed rfvers.
!    Revised help message.
!
!  2011-011-07 JDH
!    Significant Rewrite - eliminate unused coding
!    Eliminate old separate plasma response and coil response files.
!      Now will ONLY write mdsig files.
!    Eliminate interferometry and polarimetry
!
!  2011-011-08 JDH
!    Add coil group rotation and shift.

#endif
