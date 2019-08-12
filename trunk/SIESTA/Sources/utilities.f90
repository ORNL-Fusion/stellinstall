!*******************************************************************************
!>  @file utilities.f90
!>  @brief Contains module @ref utilities.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains utilities for converting the siesta grids.
!*******************************************************************************
      MODULE utilities
      USE stel_kinds
      USE island_params, ONLY: ohs=>ohs_i, ns=>ns_i, mpol=>mpol_i,             &
                               ntor=>ntor_i, nfp=>nfp_i
      USE v3_utilities, ONLY: assert, assert_eq
      USE fourier, ONLY: isym, iasym

      IMPLICIT NONE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface to convert between half and full mesh using either
!>  @ref to_full_mesh2 or @ref to_full_mesh3
!-------------------------------------------------------------------------------
      INTERFACE to_full_mesh
         MODULE PROCEDURE to_full_mesh2, to_full_mesh3
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to convert between full and half mesh using either
!>  @ref to_half_mesh or @ref to_half_mesh_p
!-------------------------------------------------------------------------------
      INTERFACE to_half_mesh
         MODULE PROCEDURE to_half_mesh, to_half_mesh_p
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to convert between full and half mesh using either
!>  @ref to_half_mesh or @ref to_half_mesh_p
!-------------------------------------------------------------------------------
      INTERFACE GradientHalf
         MODULE PROCEDURE GradientHalf, GradientHalf_p
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to take the curl on a half mesh quantity and return a full mesh
!>  using either @ref curl_htof or @ref curl_htof_targ
!-------------------------------------------------------------------------------
      INTERFACE curl_htof
         MODULE PROCEDURE curl_htof, curl_htof_targ
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculate the gradient to the half mesh from a full mesh quantity.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[inout] gradienth Gradient on the half grid.
!>  @param[in]    vecf      Full grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE GradientHalf(gradienth, vecf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: gradienth
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vecf

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(gradienth,3)
      nsmax = UBOUND(gradienth,3)

      gradienth(:,:,nsmin + 1:nsmax) = ohs*(vecf(:,:,nsmin + 1:nsmax)          &
                                     -      vecf(:,:,nsmin:nsmax - 1))
      gradienth(:,:,nsmin) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Calculate the gradient to the half mesh from a full mesh quantity.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>  Full vector is passed in as a pointer.
!>
!>  @param[inout] gradienth Gradient on the half grid.
!>  @param[in]    vecf      Full grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE GradientHalf_p(gradienth, vecf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: gradienth
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: vecf

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(gradienth,3)
      nsmax = UBOUND(gradienth,3)

      gradienth(:,:,nsmin + 1:nsmax) = ohs*(vecf(:,:,nsmin + 1:nsmax)          &
                                     -      vecf(:,:,nsmin:nsmax - 1))
      gradienth(:,:,nsmin) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Calculate the gradient to the full mesh from a half mesh quantity.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[inout] gradientf Gradient on the full grid.
!>  @param[in]    vech      Half grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE GradientFull(gradientf, vech)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: gradientf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vech

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(vech,3)
      nsmax = UBOUND(vech,3)

      gradientf(:,:,nsmin:nsmax - 1) = ohs*(vech(:,:,nsmin + 1:nsmax)          &
                                     -      vech(:,:,nsmin:nsmax - 1))

      IF (nsmin .eq. 1) THEN
         gradientf(:,:,nsmin) = 0
      END IF

      IF (nsmax .eq. ns .and. nsmax .ne. nsmin) THEN
         gradientf(:,:,nsmax) = gradientf(:,:,nsmax - 1)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Migrates a full mesh quantity to the half mesh.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[in]    vecf Full grid quantity.
!>  @param[inout] vech Half grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE to_half_mesh(vecf, vech)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vecf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: vech

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(vech,3)
      nsmax = UBOUND(vech,3)

      vech(:,:,nsmin + 1:nsmax) = 0.5*(vecf(:,:,nsmin:nsmax - 1)               &
                                +      vecf(:,:,nsmin + 1:nsmax))
      vech(:,:,nsmin) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Migrates a full mesh quantity to the half mesh.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>  Full vector grid is passed in as pointer.
!>
!>  @param[in]    vecf Full grid quantity.
!>  @param[inout] vech Half grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE to_half_mesh_p(vecf, vech)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: vecf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: vech

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(vech,3)
      nsmax = UBOUND(vech,3)

      vech(:,:,nsmin + 1:nsmax) = 0.5*(vecf(:,:,nsmin:nsmax - 1)               &
                                +      vecf(:,:,nsmin + 1:nsmax))
      vech(:,:,nsmin) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Migrate a three dimensional quantity form the half mesh to the full
!>         mesh.
!>
!>  First and last point are extrapolated linearly. The arrays must be passed in
!>  as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[in]    vech Half mesh qantity.
!>  @param[inout] vecf Full mesh quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh3(vech, vecf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vech
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: vecf

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(vecf, 3)
      CALL assert(nsmin .ge. LBOUND(vech,3), 'LBOUND WRONG IN to_full_mesh3')
      nsmax = UBOUND(vecf, 3)
      CALL assert(nsmax .le. UBOUND(vech,3), 'UBOUND WRONG IN to_full_mesh3')

      vecf(:,:,nsmin:nsmax - 1) = 0.5*(vech(:,:,nsmin + 1:nsmax)               &
                                +      vech(:,:,nsmin:nsmax - 1))
      vecf(:,:,nsmax) = 0                              
      IF (nsmin .eq. 1)  THEN
         vecf(:,:,1) = 1.5*vech(:,:,2) - 0.5*vech(:,:,3)
      END IF
      IF (nsmax .eq. ns) THEN
         vecf(:,:,ns)= 1.5*vech(:,:,ns) - 0.5*vech(:,:,ns - 1)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Migrate a two dimensional quantity form the half mesh to the full
!>         mesh.
!>
!>  First and last point are extrapolated linearly. The arrays must be passed in
!>  as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[in]    vech Half mesh qantity.
!>  @param[inout] vecf Full mesh quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh2(vech, vecf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:), ALLOCATABLE, INTENT(in)    :: vech
      REAL (dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: vecf

!  local variables
      INTEGER                                               :: nsmin
      INTEGER                                               :: nsmax

!  Start of executable code
      nsmin = LBOUND(vecf,2)
      CALL ASSERT(nsmin .ge. LBOUND(vech,2), 'LBOUND WRONG IN to_full_mesh2')
      nsmax = UBOUND(vecf,2)
      CALL ASSERT(nsmax .le. UBOUND(vech,2), 'UBOUND WRONG IN to_full_mesh2')

      vecf(:,nsmin:nsmax - 1) = 0.5*(vech(:,nsmin:nsmax - 1)                   &
                              +      vech(:,nsmin + 1:nsmax))

      vecf(:,nsmax) = 0                              
      IF (nsmin .EQ. 1)  THEN
         vecf(:,1) = vech(:,2)
      END IF
      IF (nsmax .EQ. ns) THEN
         vecf(:,ns)= vech(:,ns)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute curl of the full grid quantity a on the half grid.
!>
!>  The contravariant components of jac*Curl A are equal to
!>
!>    sqrt(g)*B^k = dA_j/du^i - dA_i/du^j                                    (j)
!>
!>  for cyclic permutations of i,j,k. This routine does not remove the jacobian
!>  factor. The arrays must be passed in as ALLOCATABLE or POINTER to preserve
!>  the array bounds.
!>
!>  @param[in]    asubsmnf  Covariant vector in the s direction full grid.
!>  @param[in]    asubsmnf  Covariant vector in the u direction full grid.
!>  @param[in]    asubsmnf  Covariant vector in the v direction full grid.
!>  @param[inout] jbsubsmnh Contravariant vector in the s direction with
!>                          jacobian factor half grid.
!>  @param[inout] jbsubsmnh Contravariant vector in the u direction with
!>                          jacobian factor half grid.
!>  @param[inout] jbsubsmnh Contravariant vector in the v direction with
!>                          jacobian factor half grid.
!>  @param[in]    iparity   Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE curl_ftoh(asubsmnf, asubumnf, asubvmnf,                       &
                           jbsupsmnh, jbsupumnh, jbsupvmnh, iparity,           &
                           nsmin, nsmax)
      USE nscalingtools, ONLY: startglobrow

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: asubsmnf
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: asubumnf
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: asubvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupvmnh
      INTEGER, INTENT(IN)                                     :: iparity
      INTEGER, INTENT(IN)                                     :: nsmin
      INTEGER, INTENT(IN)                                     :: nsmax

!  Local variables
      INTEGER                                                 :: m
      INTEGER                                                 :: n
      INTEGER                                                 :: mp
      INTEGER                                                 :: np
      INTEGER                                                 :: sparity
      INTEGER                                                 :: istat
      INTEGER                                                 :: nmin
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubumnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubvmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: da_uds
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: da_vds

!  Start of executable code
      CALL assert(nsmax .le. UBOUND(jbsupsmnh,3), 'UBOUND wrong in curl_ftoh')

      IF (iparity .EQ. isym) THEN
         sparity = 1
      ELSE
         sparity = -1
      END IF

      ALLOCATE(asubsmnh(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
               asubumnh(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
               asubvmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in curl_ftoh')

      CALL to_half_mesh(asubsmnf, asubsmnh)                       
      CALL to_half_mesh(asubumnf, asubumnh)
      CALL to_half_mesh(asubvmnf, asubvmnh)

!  Compute sqrt(g)B^X contravariant B-field components in Fourier space
      ALLOCATE(da_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                     &
               da_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation2 failed in curl_ftoh')

!  Radial gradients of asubmnu,v on half grid excluding origin and edge. Origin
!  and edge will be treated separately.
      CALL GradientHalf(da_uds, asubumnf)
      CALL GradientHalf(da_vds, asubvmnf)

      nmin = MAX(2, startglobrow)
      CALL assert(nmin .ge. LBOUND(jbsupsmnh,3), 'LBOUND wrong in curl_ftoh')

!  (sqrt(g)*B^X)mn
      DO m = 0, mpol
         mp = sparity*m
         DO n = -ntor, ntor
            np = sparity*n*nfp
            jbsupsmnh(m,n,nmin:nsmax) = np*asubumnh(m,n,nmin:nsmax)            &
                                      - mp*asubvmnh(m,n,nmin:nsmax)
            jbsupumnh(m,n,nmin:nsmax) = np*asubsmnh(m,n,nmin:nsmax)            &
                                      - da_vds(m,n,nmin:nsmax)
            jbsupvmnh(m,n,nmin:nsmax) = da_uds(m,n,nmin:nsmax)                 &
                                      - mp*asubsmnh(m,n,nmin:nsmax)
         END DO
      END DO

      DEALLOCATE (asubsmnh, asubumnh, asubvmnh, da_vds, da_uds, stat=istat)   
      CALL assert_eq(0, istat, 'Deallocation failed in curl_ftoh')

      IF (startglobrow .eq. 1) THEN
         jbsupsmnh(:,:,1) = 0
         jbsupumnh(:,:,1) = 0
         jbsupvmnh(:,:,1) = 0
      END IF

      jbsupsmnh(0,-ntor:-1,:) = 0
      jbsupumnh(0,-ntor:-1,:) = 0
      jbsupvmnh(0,-ntor:-1,:) = 0

      END SUBROUTINE


!-------------------------------------------------------------------------------
!>  @brief Compute curl of the half grid quantity a on the full grid.
!>
!>  The contravariant components of jac*Curl A are equal to
!>
!>    sqrt(g)*B^k = dA_j/du^i - dA_i/du^j                                    (j)
!>
!>  for cyclic permutations of i,j,k. This routine does not remove the jacobian
!>  factor. The arrays must be passed in as ALLOCATABLE or POINTER to preserve
!>  the array bounds.
!>
!>  @param[in]    Asubsmnh  Covariant vector in the s direction half grid.
!>  @param[in]    Asubsmnh  Covariant vector in the u direction half grid.
!>  @param[in]    Asubsmnh  Covariant vector in the v direction half grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the s direction with
!>                          jacobian factor full grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the u direction with
!>                          jacobian factor full grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the v direction with
!>                          jacobian factor full grid.
!>  @param[in]    ipar      Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE curl_htof(asubsmnh, asubumnh, asubvmnh,                       &
                           jbsupsmnf, jbsupumnf, jbsupvmnf, ipar, nsmin, nsmax)
      USE nscalingtools, ONLY: endglobrow
      USE descriptor_mod, ONLY: iam
      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: asubsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: asubumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: asubvmnh
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jbsupsmnf
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jbsupumnf
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jbsupvmnf
      INTEGER, INTENT(IN)                                  :: ipar
      INTEGER, INTENT(IN)                                  :: nsmin
      INTEGER, INTENT(IN)                                  :: nsmax

!  Local variables
      INTEGER                                              :: m
      INTEGER                                              :: n
      INTEGER                                              :: mp
      INTEGER                                              :: np
      INTEGER                                              :: mo
      INTEGER                                              :: no
      INTEGER                                              :: moff
      INTEGER                                              :: noff
      INTEGER                                              :: sparity
      INTEGER                                              :: istat
      INTEGER                                              :: nmax
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: asubsmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: asubumnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: asubvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: da_uds
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: da_vds

!  Start of executable code
      CALL assert(nsmin .ge. LBOUND(jbsupsmnf,3), 'LBOUND wrong in curl_htof')

      IF (ipar .eq. isym) THEN
         sparity = 1
      ELSE
         sparity = -1
      END IF

      moff = LBOUND(jbsupsmnf,1)
      noff = LBOUND(jbsupsmnf,2) + ntor

      ALLOCATE(asubsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               asubumnf(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               asubvmnf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in curl_htof')

      CALL to_full_mesh(asubsmnh, asubsmnf)
      CALL to_full_mesh(asubumnh, asubumnf)
      CALL to_full_mesh(asubvmnh, asubvmnf)

!  Compute sqrt(g)B^X == A^X contravariant current components in Fourier space.
      ALLOCATE(da_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                          &
               da_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation2 failed in curl_htof')

!  Radial gradients of asubmnu,v on full grid excluding origin and edge. Origin
!  and edge will be treated separately.
      CALL GradientFull(da_uds, asubumnh)
      CALL GradientFull(da_vds, asubvmnh)

!  Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
      IF (nsmin .eq. 1) THEN
         asubsmnf(2:,:,1) = 0
         asubumnf(:,:,1) = 0
         asubumnf(1,:,1) = asubumnh(1,:,2)
         da_uds(:,:,1) = 0
         da_uds(0:1,:,1) = 2*ohs*asubumnh(0:1,:,2)
         da_vds(:,:,1) = 0
         da_vds(1,:,1) = 2*ohs*asubvmnh(1,:,2)
         asubvmnf(2:,:,1) = 0
      END IF

      nmax = endglobrow
      CALL assert(nmax .le. UBOUND(jbsupsmnf,3), 'UBOUND wrong in curl_htof')

!  (sqrt(g)*J^X)mn
      DO m = 0, mpol
         mp = sparity*m
         mo = m + moff
         DO n = -ntor, ntor
            np = sparity*n*nfp
            no = n + noff
            jbsupsmnf(mo,no,nsmin:nmax) = np*asubumnf(m,n,nsmin:nmax)          &
                                        - mp*asubvmnf(m,n,nsmin:nmax)
            jbsupumnf(mo,no,nsmin:nmax) = np*asubsmnf(m,n,nsmin:nmax)          &
                                        - da_vds(m,n,nsmin:nmax)
            jbsupvmnf(mo,no,nsmin:nmax) = da_uds(m,n,nsmin:nmax)               &
                                        - mp*asubsmnf(m,n,nsmin:nmax)
         END DO
      END DO

      DEALLOCATE(asubsmnf, asubumnf,  asubvmnf, da_uds, da_vds)

!  Avoid fp failure in debug mode
      jbsupsmnf(:,:,nmax+1:UBOUND(jbsupsmnf,3)) = 0
      jbsupumnf(:,:,nmax+1:UBOUND(jbsupsmnf,3)) = 0
      jbsupvmnf(:,:,nmax+1:UBOUND(jbsupsmnf,3)) = 0

      jbsupsmnf(:,:,LBOUND(jbsupsmnf,3):nsmin-1) = 0
      jbsupumnf(:,:,LBOUND(jbsupsmnf,3):nsmin-1) = 0
      jbsupvmnf(:,:,LBOUND(jbsupsmnf,3):nsmin-1) = 0

      jbsupsmnf(moff,-ntor+noff:-1+noff,:) = 0
      jbsupumnf(moff,-ntor+noff:-1+noff,:) = 0
      jbsupvmnf(moff,-ntor+noff:-1+noff,:) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Convert allocatable to pointer for use in @ref curl_htof
!>
!>  @param[in]    Asubsmnh  Covariant vector in the s direction half grid.
!>  @param[in]    Asubsmnh  Covariant vector in the u direction half grid.
!>  @param[in]    Asubsmnh  Covariant vector in the v direction half grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the s direction with
!>                          jacobian factor full grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the u direction with
!>                          jacobian factor full grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the v direction with
!>                          jacobian factor full grid.
!>  @param[in]    ipar      Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE curl_htof_targ(asubsmnh, asubumnh, asubvmnh,                  &
                                jbsupsmnf, jbsupumnf, jbsupvmnf,               &
                                ipar, nsmin, nsmax)

!  Declare Arguments
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in) :: asubsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in) :: asubumnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in) :: asubvmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET, INTENT(inout) ::       &
         jbsupsmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET, INTENT(inout) ::       &
         jbsupumnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET, INTENT(inout) ::       &
         jbsupvmnf
      INTEGER, INTENT(IN)                                  :: ipar
      INTEGER, INTENT(IN)                                  :: nsmin
      INTEGER, INTENT(IN)                                  :: nsmax

!  Local variables
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jbsp
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jbup
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jbvp

!  Start of executable code
      jbsp => jbsupsmnf
      jbup => jbsupumnf
      jbvp => jbsupvmnf
      CALL curl_htof(asubsmnh, asubumnh, asubvmnh,                      &
                     jbsp, jbup, jbvp, ipar, nsmin, nsmax)

      END SUBROUTINE

      END MODULE
