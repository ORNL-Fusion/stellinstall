!*******************************************************************************
!>  @file ece_travis_parameters.f
!>  @brief Contains module @ref travis_parameters.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This modules contains parameters used by the ece_travis diagnostic
!*******************************************************************************
      MODULE travis_parameters

      USE stel_kinds
      IMPLICIT NONE
!>  Travis file name length 
      INTEGER, PARAMETER, PUBLIC :: trav_fnl = 256
!>  TRAVIS input string length
      INTEGER, PARAMETER, PUBLIC :: trav_isl = 32
!> Maximum number of TRAVIS antenna arrays
      INTEGER, PARAMETER :: v3fit_max_travis = 10
!********************************************************************************    
!     TRAVIS input parameters
!********************************************************************************      
!> Option to write travis output to temp folder    
      INTEGER :: travisOutTemp = 0
!> Option to write travis output to screen
      INTEGER :: travisScrOut = 0
!> Folder for TRAVIS OUPUT in travisOutTemp=1
      CHARACTER (len=trav_fnl) :: TempDirectory = "TEMP/"
!> Vessel Configuration file
      CHARACTER (len=trav_fnl) :: Vessel_Configuration = 'nofile'
!> Mirror Configuration file
      CHARACTER (len=trav_fnl) :: Mirror_Configuration = 'nofile'
!> Normalization of B
      CHARACTER (len=trav_isl) :: B0_type = 'do not scale'
      
      REAL(rprec) :: B0_ref = 1
      
      REAL(rprec) :: phi_ref = 0

      REAL(rprec) :: B_scale = 1

      REAL(rprec) :: hgrid = 0.021

      REAL(rprec) :: dphi = 2
!*****************************************************************************
!     TRAVIS antenna inputs
!*****************************************************************************
!>    Array of antenna positions
      REAL(rprec),DIMENSION(v3fit_max_travis,3)                         &
     &                     :: antenna_position=2
!>    Array of antenna positions
      REAL(rprec),DIMENSION(v3fit_max_travis,3)                         &
     &                     :: target_position =1
!>    Coordinate type of the antenna (assume all antenna use the same coord)
      CHARACTER (len=trav_isl) :: antenna_coordinate = 'cyl'
!>    Coordinate type of the targer (assume all targets use the same coord)
      CHARACTER (len=trav_isl) :: target_coordinate = 'cyl'
!>    Beam radiuses
      REAL(rprec),DIMENSION(v3fit_max_travis,2) :: rbeam = 0.02
!>    Beam-ellipse rotation with respect to the equator plan (deg)
      REAL(8),DIMENSION(v3fit_max_travis) :: phibx =0
!>    Beam focal lengths
      REAL(rprec),DIMENSION(v3fit_max_travis,2) :: rfocus = 10.0
!>    whether to emulate quasi-optical definition
      INTEGER(4),DIMENSION(v3fit_max_travis) :: QOemul =0
!>    whether to emulate quasi-optical definition
      INTEGER(4),DIMENSION(v3fit_max_travis) :: nra =0
!>    whether to emulate quasi-optical definition
      INTEGER(4),DIMENSION(v3fit_max_travis) :: nphi =8
 
      END MODULE
