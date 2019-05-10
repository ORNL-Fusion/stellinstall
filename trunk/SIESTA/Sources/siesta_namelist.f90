!-------------------------------------------------------------------------------
!  The @header, @table_section, @table_subsection, @item and @end_table commands
!  are custom defined commands in Doxygen.in. They are defined under ALIASES.
!  For the page created here, the 80 column limit is exceeded. Arguments of
!  aliases are separated by ','. If you intended ',' to be a string you must use
!  an escaped comma '\,'.
!
!>  @page siesta_namelist_sec Namelist siesta_main_nli definition
!>  
!>  @tableofcontents
!>  @section siesta_namelist_intro_sec Introduction
!>  This page documents the contents of a SIESTA namelist input file. SIESTA 
!>  namelist variables are defined in the @fixed_width{siesta_info} common 
!>  block.
!>
!>  @section siesta_namelist_var_sec Namelist Variables
!>  @header{Input variable, Description, Code Reference}
!>
!>  @table_section{siesta_flag_sec, Control Flags}
!>     @item{ladd_pert,        Use helical perturbation.}
!>     @item{lresistive,       Use resistive perturbation.}
!>     @item{lrestart,         Restart SIESTA from pervious run.}
!>     @item{l_tracing,        Produce output file for fieldline tracing.}
!>     @item{lcolscale,        Use column scaling.}
!>     @item{l_silo_output,    Produce silo output.}
!>     @item{l_silo3D,         Produce silo 3D output.}
!>     @item{l_output_alliter, Write output files on all iterations instead of
!>                             only iterations that lower the MHD energy and
!>                             force residuals.}
!>     @item{l_VMEC_Uniform,   }
!>     @item{lasym,            Use non stellarator symmetric terms.}
!>     @item{l_vessel,         Use free-boundary grid extension to vessel wall.}
!>     @item{lrecon,           Add additional output to the restart file when
!>                             used in a reconstruction context.}
!>  @end_table
!>
!>  @table_section{siesta_algrothim_sec, Algrothim Control Variables}
!>     @item{niter,         Maximum number of iterations after diagonal prec.}
!>     @item{ftol,          Minimum force tolarance for converged solution.}
!>     @item{mupar,         }
!>     @item{levmarq_param, Inital value of Levenberg-Marquardt parameter.}
!>     @item{eta_factor,    Resistivity value.}
!>     @item{nprecon,       Skip diagonal preconditioner if greater than zero.}
!>     @item{ngmres_type,   }
!>     @item{iortho,        }
!>  @end_table
!>
!>  @table_section{siesta_island_sec, Island Parameters}
!>     @item{mres,    M numbers of island resonances.}
!>     @item{helpert, Sizes of the helical perturbation.}
!>     @table_subsection{siesta_island_dep_sec, Deprecated Parameters}
!>        @item{mres1,    First M number of island resonance.}
!>        @item{nres1,    First N number of island resonance.}
!>        @item{rad1,     First radius of island resonance.}
!>        @item{helpert1, First size of the helical perturbation.}
!>        @item{mres2,    Second M number of island resonance.}
!>        @item{nres2,    Second N number of island resonance.}
!>        @item{rad2,     Second radius of island resonance.}
!>        @item{helpert2, Second size of the helical perturbation.}
!>  @end_table
!>
!>  @table_section{siesta_grid_size_sec, Grid Sizes}
!>     @item{nsin,     Size of plasma radial grid.}
!>     @item{nsin_ext, Size of extended radial grid.}
!>     @item{mpolin,   Number of poloidal modes.}
!>     @item{ntorin,   Number of toroidal modes.}
!>     @item{nfpin,    Number of field periods to use. Setting this to a
!                      anything less than one will use the value form the wout
!                      file.}
!>     @table_subsection{siesta_grid_size_out_sec, Output Grid Sizes}
!>        @item{nphis, Number of cylindrical phi planes.}
!>        @item{nrs,   Number of radial grid points.}
!>        @item{nzs,   Number of vertical grid points.}
!>        @item{nvs,   Number of flux space toroidal points.}
!>        @item{nus,   Number of flux space poloidal points.}
!>        @item{nss,   Number of flux space radial points.}
!>  @end_table
!>
!>  @table_section{siesta_file_name_sec, File Names}
!>      @item{wout_file,   Filename of the VMEC woutfile.}
!>      @item{restart_ext, Name of the restart file extension.}
!>      @item{mgrid_file,  Filename of the MGRID file.}
!>      @item{vessel_file, Filename of the extended surfaces.}
!>  @end_table
!>
!>  @table_section{siesta_test_sec, Test Controls}
!>      @item{hesspass_test, }
!>      @item{mupar_test,    }
!>  @end_table
!>
!>  @section siesta_namelist_prog_ref_sec Programmers Reference
!>  Reference material for the coding to implement this namelist is found in the
!>  @ref siesta_namelist module.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file siesta_namelist.f90
!>  @brief Contains module @ref siesta_namelist.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains all the variables and maximum sizes of the inputs for a
!>  SIESTA namelist input file. The module contained within does not represent 
!>  an object instance. Instead all variables are contained in a global context.
!>  This is required due to limitations of FORTRAN 95 and namelist inputs.
!>
!>  @ref siesta_namelist_sec "Namelist siesta_info definition"
!>
!>  @note Some of the references are missing here. This is due to a bug in
!>  Doxygen when variable decalarations span multiple lines.
!*******************************************************************************
      MODULE siesta_namelist
      USE shared_data, ONLY: ngmres_type, iortho, lcolscale, lasym,            &
                             hesspass_test, mupar_test, lrecon
      USE Hessian, ONLY: levmarq_param, mupar
      USE stel_kinds

      IMPLICIT NONE

!*******************************************************************************
!  siesta_namelist input module parameters
!*******************************************************************************
!>  Input string length.
      INTEGER, PARAMETER ::  siesta_namelist_name_length = 256

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) siesta_namelist_class
!
!*******************************************************************************
!  Control Flags
!>  Use helical perturbation.
      LOGICAL :: ladd_pert = .TRUE.
!>  Use resistive perturbaton.
      LOGICAL :: lresistive = .TRUE.
!>  Restart SIESTA from pervious run.
      LOGICAL :: lrestart = .FALSE.
!>  Produce output file for fieldline tracing.
      LOGICAL :: l_tracing = .FALSE.
!>  Produce silo output.
      LOGICAL :: l_silo_output = .FALSE.
!>  Produce silo 3D output.
      LOGICAL :: l_silo3D = .FALSE.
!>  Write output files on all iterations.
      LOGICAL :: l_output_alliter = .FALSE.
!>  FIXME: Unknown
      LOGICAL :: l_VMEC_Uniform
!>  If extended grid is to be used using an available vessel file
      LOGICAL :: l_vessel = .FALSE.

!  Algrothim Control Variables}
!>  Maximum number of iterations after diagonal prec.
      INTEGER  :: niter = 10
!>  Force tolarance.
      REAL(dp) :: ftol = 1.E-20_dp
!>  Resistivity value.
      REAL(dp) :: eta_factor = 1.E-2_dp
!>  Skip diagonal preconditioner if greater than zero.
      INTEGER  :: nprecon = 0

!  Island parameters
!>  Sizes of the helical perturbation.
      INTEGER, DIMENSION(20)  :: mres = 0
!>  Sizes of the helical perturbation.
      REAL(dp), DIMENSION(20) :: HelPert = 0.0
!>  Sizes of the helical perturbation.
      REAL(dp), DIMENSION(20) :: HelPertA = 0.0

!  Grid Sizes
!>  Radial size of the plasma grid.
      INTEGER :: nsin = 101
!>  Radial size of the extended grid.
      INTEGER :: nsin_ext = 0
!>  Number of poloidal modes.
      INTEGER :: mpolin = 12
!>  Number of toroidal modes.
      INTEGER :: ntorin = 3
!>  Number of field periods to use. -1 means set this to the value in the wout
!>  file
      INTEGER :: nfpin = 0

!  Output Grid Sizes
!>  Number of cylindrical phi planes.
      INTEGER :: nphis = 2
!>  Number of radial grid points.
      INTEGER :: nrs=200
!>  Number of vertical grid points.
      INTEGER :: nzs=200
!>  Number of flux space toroidal points.
      INTEGER :: nvs=150
!>  Number of flux space poloidal points.
      INTEGER :: nus=150
!>  Number of flux space radial points.
      INTEGER :: nss=100

!  File Names
!>  Filename of the VMEC woutfile.
      CHARACTER(LEN=siesta_namelist_name_length) :: wout_file = ''
!>  Name of the restart file extension.
      CHARACTER(LEN=siesta_namelist_name_length) :: restart_ext = ''
!>  Filename of the VMEC woutfile.
      CHARACTER(LEN=siesta_namelist_name_length) :: mgrid_file = ''
!>  Name of the restart file extension.
      CHARACTER(LEN=siesta_namelist_name_length) :: vessel_file = ''

! Deprecated Parameters
!>  @deprecated
!>  First M number of island resonance.
      INTEGER  :: mres1
!>  @deprecated
!>  First N number of island resonance.
      INTEGER  :: nres1
!>  @deprecated
!>  First radius of island resonance.
      REAL(dp) :: rad1
!>  @deprecated
!>  First size of the helical perturbation.
      REAL(dp) :: HelPert1
!>  @deprecated
!>  Second M number of island resonance.
      INTEGER  :: mres2
!>  @deprecated
!>  Second N number of island resonance.
      INTEGER  :: nres2
!>  @deprecated
!>  Second radius of island resonance.
      REAL(dp) :: rad2
!>  @deprecated
!>  Second size of the helical perturbation.
      REAL(dp) :: HelPert2

!  Declare Namelist
      NAMELIST/siesta_info/                                             &
!  Control flags
        ladd_pert, lresistive, lrestart, l_tracing, lcolscale,          &
        l_silo_output, l_silo_output, l_silo3D, l_output_alliter,       &
        l_VMEC_Uniform, lasym, lrecon, l_vessel,                        &
!  Algrothim Control Variables
        niter, ftol, mupar, levmarq_param, eta_factor, nprecon,         &
        ngmres_type, iortho,                                            &
!  Island parameters Island Parameters
        mres, HelPert, HelPertA,                                        &
!  Input grid sizes
        nsin, nsin_ext, mpolin, ntorin, nfpin,                          &
!  Output grid sizes
        nphis, nrs, nzs, nvs, nus, nss,                                 &
!  File names
        wout_file, restart_ext, mgrid_file, vessel_file,                &
!  Test controls
        hesspass_test, mupar_test,                                      &
!  Deprecated Parameters
        mres1, nres1, rad1, HelPert1, mres2, nres2, rad2, HelPert2

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Reads the namelist input file.
!>
!>  Reads the namelist input file.
!>
!>  @param[in] namelist_file The file name of the namelist input file.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_namelist_read(namelist_file)
      USE safe_open_mod
      USE v3_utilities

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: namelist_file

!  local variables
      INTEGER                       :: iou_mnli
      INTEGER                       :: status

!  Start of executable code

!  Initalize a default value of the I\O unit. SIESTA increments from there.
      iou_mnli = 0
      CALL safe_open(iou_mnli, status, TRIM(namelist_file),             &
                    'old', 'formatted')
      CALL assert_eq(0, status, 'siesta_namelist_read' //               &
        ': Safe_open of ' // TRIM(namelist_file) // ' failed')

!  Read the namelist input file.
      READ (iou_mnli, nml=siesta_info)
      CLOSE (iou_mnli, iostat=status)
      CALL assert_eq(0, status, 'siesta_namelist_read' //               &
        ': Error closing ' // TRIM(namelist_file) // ' failed')

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Writes the namelist input file.
!>
!>  Writes the namelist input file.
!>
!>  @param[in] namelist_file The file name of the namelist input file.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_namelist_write(namelist_file)
      USE safe_open_mod
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: namelist_file

!  local variables
      INTEGER                       :: iou_mnli
      INTEGER                       :: status

!  Start of executable code

!  Initalize a default value of the I\O unit. SIESTA increments from there.
      iou_mnli = 0
      CALL safe_open(iou_mnli, status, TRIM(namelist_file),                    &
     &               'replace', 'formatted', delim_in='quote')
      CALL assert_eq(0, status, 'siesta_namelist_write' //                     &
     &   ': Safe_open of ' // TRIM(namelist_file) // ' failed')

!  Write the namelist input file.
      WRITE (iou_mnli, nml=siesta_info)
      CLOSE (iou_mnli, iostat=status)
      CALL assert_eq(0, status, 'siesta_namelist_read' //                      &
     &   ': Error closing ' // TRIM(namelist_file) // ' failed')

      END SUBROUTINE

      END MODULE
