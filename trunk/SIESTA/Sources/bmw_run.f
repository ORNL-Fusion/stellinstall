!*******************************************************************************
!>  @file bmw.f
!>  @brief Contains the main routines for Biot-Savart Magnetic Vmec
!>  Vector-Potential (BMW) code.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  BMW is a code for extending fields belond the VMEC domain in a manner that
!>  ensures divergence free fields. BMW does a Biot-Savart volume integration of
!>  of the equilibrium current density to obtain a continous vector potential
!>  every where on the mgrid grid.
!>
!>  @author Mark Cianciosa
!>
!>  Below is a brief discription of the major top level objects of the code. For
!>  discriptions of lower level objects consult the referenced top level
!>  objects.
!>
!>  @ref m_grid        Object containing the vaccum field information.
!>  @ref primed_grid   Object containing the plasma currents and primed grid
!>                     positions.
!>  @ref unprimed_grid Object containing the plasma vector potential response.
!>  @ref bmw_context   Defines the main pointers to all allocated memory and
!>                     objects.
!*******************************************************************************
!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief BMW main program.
!>
!>  Highest level BMW routine. This computes the vector potential on the mgrid
!>  grid.
!-------------------------------------------------------------------------------
      MODULE bmw_run
      USE bmw_context
      USE bmw_commandline_parser
      USE, INTRINSIC :: iso_fortran_env, Only : output_unit
      USE bmw_state_flags

      IMPLICIT NONE

      CONTAINS 

      SUBROUTINE bmw_exec(mgrid_file_name, wout_file_name,                     &
     &                    r_grid, z_grid, dphi, surf_num,                      &
     &                    A_r, A_p, A_z)
      USE read_wout_mod, ONLY: nfp_vmec=>nfp
      USE island_params, ONLY: nfp_i

!  dummy variables
      CHARACTER (len=*), INTENT(in)                :: mgrid_file_name
      CHARACTER (len=*), INTENT(in)                :: wout_file_name
      REAL (rprec), DIMENSION(:,:,:), INTENT(in)   :: r_grid
      REAL (rprec), DIMENSION(:,:,:), INTENT(in)   :: z_grid
      REAL (rprec), INTENT(in)                     :: dphi
      INTEGER, INTENT(in)                          :: surf_num
      REAL (rprec), DIMENSION(:,:,:), INTENT(out)  :: A_r
      REAL (rprec), DIMENSION(:,:,:), INTENT(out)  :: A_p
      REAL (rprec), DIMENSION(:,:,:), INTENT(out)  :: A_z

!  local variables
      TYPE (bmw_parallel_context_class), POINTER   :: parallel => null()
      TYPE (bmw_context_class), POINTER            :: context => null()
      TYPE (bmw_commandline_parser_class), POINTER ::                          &
     &   cl_parser => null()
      INTEGER                                      :: flags
      INTEGER                                      :: num_p
      INTEGER                                      :: status
      INTEGER                                      :: i

!  Start of executable code
      CALL profiler_construct

      parallel => bmw_parallel_context_construct(MPI_COMM_WORLD)
      flags = bmw_state_flags_off

      CALL bmw_parallel_context_set_threads(parallel, -1)
      CALL bmw_parallel_context_report(parallel, output_unit)

!  BMW uses num_p to set the number of grid points in for the primed grid. In
!  doing so, it mutliplies by the number of field periods in the mgrid file. If
!  using a user defined number of field periods that is less than the amount in
!  in mgrid, we need to scale num_p.
      num_p = SIZE(A_r, 2)*2/(nfp_vmec/nfp_i)
      context => bmw_context_construct(mgrid_file_name,                        &
     &                                 wout_file_name, '', flags,              &
     &                                 num_p, parallel, output_unit)

      CALL bmw_context_set_up_grid(context,                                    &
     &                             r_grid(:,:,surf_num:surf_num),              &
     &                             z_grid(:,:,surf_num:surf_num),              &
     &                             dphi, parallel, output_unit)

      A_r = context%up_grid%a_r
      A_p = context%up_grid%a_p
      A_z = context%up_grid%a_z

      CALL bmw_context_destruct(context)
      CALL profiler_destruct

      CALL bmw_parallel_context_destruct(parallel                              &
#if defined (MPI_OPT)
     &                                   ,.false.                              &
#endif
     &                                  )

1000  FORMAT('BMW ',i4,' Series.')

      END SUBROUTINE bmw_exec

      END MODULE bmw_run
