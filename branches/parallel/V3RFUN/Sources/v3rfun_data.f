      MODULE v3rfun_data

      USE vsvd0, ONLY: nigroup
      USE stel_kinds
      USE bsc_T
      USE diagnostic_T

!  Declarations/Initializations of variables used in namelist module
      CHARACTER(len=120)           :: name_coils_dot=' '
      CHARACTER(len=120)           :: name_diagnostic_dot=' '
      CHARACTER(len=20)           :: idrfun=' '
      CHARACTER(len=20)           :: myid=' '
      LOGICAL                     :: lstell_sym=.false.
!  JDH 6.26.03 Changed default value of lstell_sym from .true to .false.
!  See Stellarator Symmetry Powerpoint of 6.10.03
      REAL(rprec)                 :: rmin = 0, rmax = 0 
      REAL(rprec)                 :: zmin = 0, zmax = 0
      REAL(rprec)                 :: len_integrate_ddc = 0.001
      INTEGER                     :: ir = 101, jz = 101, kp = 1
      INTEGER                     :: kp_store
      LOGICAL                     :: l_read_coils_dot=.true.
!  name_coils_dot        name of the 'coils.' file
!  name_diagnostic_dot   name of the 'diagnostic.' file
!  idrfun                identification for the run
!  myid                  supplemental id for user 
!  lstell_sym            Logical - True for stellarator symmetry
!  l_read_coils_dot      Logical, True to have v3rfun read a coils_dot file.
!                           If false, coil-diagnostic response not calculated.
!  rmin                  Minimum R for plasma grid
!  rmax                  Maximum R for plasma grid
!  zmin                  Minimum z for plasma grid
!  zmax                  Maximum z for plasma grid
!  ir                    number of grid points in R, plasma grid
!  jz                    number of grid points in z, plasma grid
!  kp                    number of phi planes per field period in plasma grid
!  kp_store              number of phi planes actually store in plasma grid
!                           (With lstell_sym = .true., don't need to store all planes)
!  len_integrate_ddc     length for integration over the diagnostic_dot coils
!                        Spacing for integration is determined by this length(m)
!                        At least one integration point per segment
!                        Also minimum 60 points for circular diagnostics 
!                        (set in bsc_flux_pos)
!  Numbers of coils and coil groups (Defined after input completed)
      INTEGER              :: n_diagn_c, n_field_cg

!  Number of field periods.
!     _nli - read in from the namelist input
!     _cd  - read in from the coils_dot file
!  Logic to determine which is used is in subroutine get_input
      INTEGER              :: n_field_periods
      INTEGER              :: n_field_periods_nli = 0
      INTEGER              :: n_field_periods_cd

!  Store the coil groups in an ALLOCATABLE array of TYPE bsc_coilcoll
!  - field_coils(i) holds all the coils in coil group i

!  JDH 2011-09-09. Use module biotsavart to parse coils file
!    Have field_coils point to module biotsavart variable
      TYPE (bsc_coilcoll), DIMENSION(:), POINTER :: field_coils

!  Store the diagnostic coils in a bsc_coilcoll
!  JDH 10.23.02. Compiler on logjam does not like this statement.
! '  2070-S: "rfun_module.f", line 105: diagnostic_coils shall have the SAVE attribute;'
! ' has type for which component initialization was specified in a MODULE.'
!  So, make it SAVE
!  2012-01-12. Eliminate diagnostic_coils. Information is in mddc_desc_array

!      TYPE (bsc_coilcoll), SAVE :: diagnostic_coils

!  Store the mddc_desc-s in an array.
      TYPE(mddc_desc), DIMENSION(10000), SAVE :: mddc_desc_array

!----------------------------------------------------------------------
!-- response functions due to external coils                         --
!--   rdiag_coilg         : diagnostic coils due to coil groups      --
!--   inductance          : inductance (normalized rdiag_coilg) 
!--                         due to coil groups                       --
!----------------------------------------------------------------------
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: rdiag_coilg
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: inductance

!  Plasma Response functions

!----------------------------------------------------------------------
!-- input/output unit numbers
!      Not parameters, in case get changed somewhere.                  
!--   iou_nli           : namelist input file                        --
!--   iou_out           : output file (text)                         --
!--   iou_list          : .LIST file for the diagnostic coils        --
!--   iou_mi            : .MI file for the mutual inductances        --
!--   iou_mdsig_ncdf    : mdsig files, netCDF                        --
!----------------------------------------------------------------------
      INTEGER :: iou_nli   = 11
      INTEGER :: iou_out   = 10
      INTEGER :: iou_list  = 35
      INTEGER :: iou_mi    = 37
      INTEGER :: iou_mdsig_ncdf  = 40

!----------------------------------------------------------------------
!-- Conducting Shell specification                                   --
!-- For now model this as a torus.                                   --
!--   use_con_shell : Computes the response function for a           --
!--                   conducting shell.                              --
!--   major_radius  : Major radius of the torus.                     --
!--   minor_radius  : Minor radius of the torus.                     --
!----------------------------------------------------------------------
      LOGICAL      :: use_con_shell = .false.
      REAL (rprec) :: major_radius = 0.0
      REAL (rprec) :: minor_radius = 0.0
      INTEGER      :: kp_shell = 1
      INTEGER      :: kp_shell_store

!----------------------------------------------------------------------
!-- version                                                          --
!--   rfvers            : version number                             --
!----------------------------------------------------------------------
!  NOTE: When changing this version changes, the select statement in
!  magnetic_construct_netcdf of magnetic.f in the V3FIT needs to be updated.
      CHARACTER(len=30) :: rfvers = 'MRC 2014-09-28'

!---------------------------------------------------------------------------
!  Variables for shifts and rotations of coil_groups of coils_dot file
!
!        (  : indicates dimension nextcur_dim)
!    cg_shift_1(:,3)        Vector to shift all the coils. (Before rotation)
!    cg_rot_theta(:)        Spherical polar angle to specify axis of rotation
!    cg_rot_phi(:)          Spherical azimuthal angle to specify axis of rotation
!    cg_rot_angle(:)        Angle to rotate about axis of rotation. 
!                             NB - LEFT HAND convention. Put left thumb along 
!                             axis of rotation, fingers indicate direction of 
!                             positive rotation.
!    cg_rot_xcent(:,3)      Position of center of rotation
!    l_rot_coil_center(:)   Logical. True - use current-averaged center of
!                             coil-group for center of rotation
!                             False - use position specified in cg_rot_xcent
!                             for center of rotation
!    cg_shift_2(:,3)        Vector to shift all the coils. (After rotation)    
!    
!---------------------------------------------------------------------------
      INTEGER, PARAMETER :: nextcur_dim = nigroup
      REAL(rprec), DIMENSION(nextcur_dim,3) ::  cg_shift_1,                    &
     &   cg_shift_2, cg_rot_xcent
      REAL(rprec), DIMENSION(nextcur_dim) :: cg_rot_theta,                     &
     &   cg_rot_phi, cg_rot_angle
      LOGICAL, DIMENSION(nextcur_dim) :: l_rot_coil_center
      
!----------------------------------------------------------------------
! Namelist v3r_coils.  
!----------------------------------------------------------------------
      NAMELIST/v3r_coils/ idrfun, name_coils_dot,                              &
     &   name_diagnostic_dot,                                                  &
     &   lstell_sym, rmin, rmax, zmin, zmax, ir, jz, kp,                       &
     &   len_integrate_ddc,                                                    &
     &   l_read_coils_dot, n_field_periods_nli,                                &
     &   cg_shift_1, cg_shift_2, cg_rot_xcent, cg_rot_theta,                   &
     &   cg_rot_phi, cg_rot_angle, l_rot_coil_center,                          &
     &   use_con_shell, major_radius, minor_radius, kp_shell

      END MODULE v3rfun_data
