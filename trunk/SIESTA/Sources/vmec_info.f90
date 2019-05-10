!*******************************************************************************
!>  @file vmec_info.f90
!>  @brief Contains module @ref vmec_info
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Holds memory references for splined and island VMEC quantities.
!*******************************************************************************
      MODULE vmec_info
      USE v3_utilities, ONLY: assert
      USE stel_kinds

      IMPLICIT NONE

!*******************************************************************************
!  vmec_info module variables
!*******************************************************************************
!>  Splined Fourier coefficents for VMEC R cosine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rmnc_spline
!>  Splined Fourier coefficents for VMEC Z sine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: zmns_spline
!>  Splined Fourier coefficents for VMEC lambda sine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: lmns_spline
!>  Splined Fourier coefficents for VMEC R sine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rmns_spline
!>  Splined Fourier coefficents for VMEC Z cosine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: zmnc_spline
!>  Splined Fourier coefficents for VMEC lambda cosine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: lmnc_spline
!>  Splined Fourier coefficents for VMEC J^u cosine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: currumnc_spline
!>  Splined Fourier coefficents for VMEC J^v cosine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: currvmnc_spline
!>  Splined Fourier coefficents for VMEC J^u sine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: currumns_spline
!>  Splined Fourier coefficents for VMEC J^v sine parity.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: currvmns_spline

!>  Fourier coefficents for VMEC R cosine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: rmnc_i
!>  Fourier coefficents for VMEC Z sine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zmns_i
!>  Fourier coefficents for VMEC lambda sine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: lmns_i
!>  Fourier coefficents for VMEC R sine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: rmns_i
!>  Fourier coefficents for VMEC Z cosine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zmnc_i
!>  Fourier coefficents for VMEC lambda cosine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: lmnc_i
!>  Fourier coefficents for VMEC J^u cosine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: jcurrumnc
!>  Fourier coefficents for VMEC J^v cosine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: jcurrvmnc
!>  Fourier coefficents for VMEC J^u sine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: jcurrumns
!>  Fourier coefficents for VMEC J^v sine parity.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: jcurrvmns

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Allocate spline buffers.
!>
!>  Allocated memory of VMEC quantities spline buffers. Spline arrays take the
!>  dimensions of the total number of Fourier modes and the size of the radial
!>  grid. Lambda is allocated with an extra radial grid point.
!>
!>  @param[in] mnmax   Total number of Fourier modes.
!>  @param[in] ns      Total number of radial surfaces.
!>  @param[in] is_asym Flag to mark if not using stellarator symmetry.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_construct_spline(mnmax, ns, is_asym)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: ns
      LOGICAL, INTENT(in) :: is_asym

!  Start of executable code

!  Stellarator symmetric quantities.
      ALLOCATE(rmnc_spline(mnmax,ns))
      ALLOCATE(zmns_spline(mnmax,ns))
      ALLOCATE(lmns_spline(mnmax,ns + 1))
      ALLOCATE(currumnc_spline(mnmax,ns))
      ALLOCATE(currvmnc_spline(mnmax,ns))

!  Asymmetric quantities.
      IF (is_asym) THEN
         ALLOCATE(rmns_spline(mnmax,ns))
         ALLOCATE(zmnc_spline(mnmax,ns))
         ALLOCATE(lmnc_spline(mnmax,ns + 1))
         ALLOCATE(currumns_spline(mnmax,ns))
         ALLOCATE(currvmns_spline(mnmax,ns))
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Allocate island buffers.
!>
!>  Allocated memory of VMEC quantities island buffers. Island arrays take the
!>  dimensions of 0:mpol, -ntor:ntor and the size of the radial grid. Lambda is
!>  allocated with an extra radial grid point.
!>
!>  @param[in] mpol    Number of poloidal modes.
!>  @param[in] ntor    Number of Toroidal modes.
!>  @param[in] ns      Number of radial surfaces.
!>  @param[in] is_asym Flag to mark if not using stellarator symmetry.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_construct_island(mpol, ntor, ns, is_asym)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in) :: mpol
      INTEGER, INTENT(in) :: ntor
      INTEGER, INTENT(in) :: ns
      LOGICAL, INTENT(in) :: is_asym

!  Start of executable code

!  Stellarator symmetric quantities.
      ALLOCATE(rmnc_i(0:mpol,-ntor:ntor,ns))
      ALLOCATE(zmns_i(0:mpol,-ntor:ntor,ns))
      ALLOCATE(lmns_i(0:mpol,-ntor:ntor,ns + 1))
      ALLOCATE(jcurrumnc(0:mpol,-ntor:ntor,ns))
      ALLOCATE(jcurrvmnc(0:mpol,-ntor:ntor,ns))

      rmnc_i = 0
      zmns_i = 0
      lmns_i = 0
      jcurrumnc = 0
      jcurrvmnc = 0

!  Asymmetric quantities.
      IF (is_asym) THEN
         ALLOCATE(rmns_i(0:mpol,-ntor:ntor,ns))
         ALLOCATE(zmnc_i(0:mpol,-ntor:ntor,ns))
         ALLOCATE(lmnc_i(0:mpol,-ntor:ntor,ns + 1))
         ALLOCATE(jcurrumns(0:mpol,-ntor:ntor,ns))
         ALLOCATE(jcurrvmns(0:mpol,-ntor:ntor,ns))

         rmns_i = 0
         zmnc_i = 0
         lmnc_i = 0
         jcurrumns = 0
         jcurrvmns = 0
      END IF

      END SUBROUTINE

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deallocate spline buffers.
!>
!>  Deallocated memory of VMEC quantities spline buffers.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_destruct_spline

      IMPLICIT NONE

!  Start of executable code

!  Stellarator symmetric quantities.
      DEALLOCATE(rmnc_spline)
      DEALLOCATE(zmns_spline)
      DEALLOCATE(lmns_spline)
      DEALLOCATE(currumnc_spline)
      DEALLOCATE(currvmnc_spline)

!  Asymmetric quantities.
      IF (ALLOCATED(rmns_spline)) THEN
         DEALLOCATE(rmns_spline)
         DEALLOCATE(zmnc_spline)
         DEALLOCATE(lmnc_spline)
         DEALLOCATE(currumns_spline)
         DEALLOCATE(currvmns_spline)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deallocate island buffers.
!>
!>  Deallocated memory of VMEC quantities island buffers.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_destruct_island

      IMPLICIT NONE

!  Start of executable code

!  Stellarator symmetric quantities.
      DEALLOCATE(rmnc_i)
      DEALLOCATE(zmns_i)
      DEALLOCATE(lmns_i)
      DEALLOCATE(jcurrumnc)
      DEALLOCATE(jcurrvmnc)

!  Asymmetric quantities.
      IF (ALLOCATED(rmns_i)) THEN
         DEALLOCATE(rmns_i)
         DEALLOCATE(zmnc_i)
         DEALLOCATE(lmnc_i)
         DEALLOCATE(jcurrumns)
         DEALLOCATE(jcurrvmns)
      END IF

      END SUBROUTINE

      END MODULE vmec_info
