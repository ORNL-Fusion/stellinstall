!>
!!  \mainpage SIESTA
!!  \brief Program for computing 3D MHD equilibria including magnetic islands
!!  \author S. P. Hirshman, R. Sanchez, and C. R. Cook (Physics of Plasmas 18, 062504, 2011)
!!  \author S. K. Seal, K. S. Perumalla and S. P. Hirshman (Concurrency and Computation: Practice and Experience, 25(15), p 2207-2223, 2013)
!!  \version 2.3 (June 2014)
!!  \bug Please report any bugs to hirshmansp@ornl.gov
      PROGRAM SIESTA
      USE siesta_run
      USE mpi_inc

      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      TYPE (siesta_run_class), POINTER :: context => NULL()

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
#if defined(MPI_OPT)
      INTEGER, PARAMETER :: run_comm = MPI_COMM_WORLD
#else
      INTEGER, PARAMETER :: run_comm = 0
#endif
!-----------------------------------------------

      context => siesta_run_construct(run_comm, .true.)
      CALL siesta_run_set_vmec(context)
      CALL siesta_run_converge(context)
      CALL siesta_run_destruct(context, .true., .true.)

      END PROGRAM SIESTA
