
!*******************************************************************************
!>  @file siesta.f90
!>  @brief Contains main routines for SIESTA.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  SIESTA is an ideal MHD equilibrium code that allows for islands and
!>  stocastic fields. SIESTA has been created with contributions from.
!>
!>  @author Steve Hirshman
!>  @author R. Sanchez
!>  @author Sudip Seal
!>  @author Mark Cianciosa
!>
!>  The theory behind SIESTA is outlined in
!>     S.P. Hirshman et. al. https://doi.org/10.1063/1.3597155
!>  The design of the parallelization method is outlined in
!>     S.K. Seal at. al. https://doi.org/10.1002/cpe.2919
!>
!>  Please report any bugs to ciancisamr@ornl.gov
!*******************************************************************************
!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief SIESTA Main program.
!>
!>  Highest level SIESTA routine.
!-------------------------------------------------------------------------------
      PROGRAM SIESTA
      USE siesta_run
      USE mpi_inc
      USE siesta_namelist, ONLY: lrestart

      IMPLICIT NONE

!  Local variables
      TYPE (siesta_run_class), POINTER :: context => NULL()

!  Local parameters
#if defined(MPI_OPT)
      INTEGER, PARAMETER               :: run_comm = MPI_COMM_WORLD
#else
      INTEGER, PARAMETER               :: run_comm = 0
#endif

!  Start of executable code.
!  Create a context and read the SIESTA namelist input file.
      context => siesta_run_construct(run_comm, .true.)

!  Check the restart flag in from the namelist input. If true, restart from an
!  existing restart file. Otherwise initalize from a VMEC file.
      IF (lrestart) THEN
         CALL siesta_run_set_restart(context)
      ELSE
         CALL siesta_run_set_vmec(context)
      END IF

      CALL siesta_run_converge(context)
      CALL siesta_run_destruct(context, .true., .true.)

      END PROGRAM SIESTA

     ! siesta_run_set_vmec 
     ! -> init_metric_elements (in metrics)
     ! -> grid_extender (in grid_extension)
     ! siesta_run_converge
     ! -> run_pchelms 
