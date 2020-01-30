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
      USE fourier, ONLY: m0, m1, m2, n0, n1, f_sin, f_cos

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

!-------------------------------------------------------------------------------
!>  Interface to set the boundary conditions on a full mesh quantity converted
!>  from the half mesh for either a scalar @ref set_bndy_half_to_full or
!>  vector @ref set_bndy_half_to_full_vec quantities.
!-------------------------------------------------------------------------------
      INTERFACE set_bndy_half_to_full
         MODULE PROCEDURE set_bndy_half_to_full, set_bndy_half_to_full_vec
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to set the boundary conditions on a full mesh quantity converted
!>  from the half mesh for either a scalar @ref set_bndy_half_to_full or
!>  vector @ref set_bndy_half_to_full_vec quantities.
!-------------------------------------------------------------------------------
      INTERFACE set_bndy_half_to_full_ds
         MODULE PROCEDURE set_bndy_half_to_full_ds, set_bndy_half_to_full_uv_ds
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to set the fouier conditions on a quantity for either scalar
!>  @ref set_bndy_fouier_m0 or vector @ref set_bndy_fouier_m0_vec quantities.
!-------------------------------------------------------------------------------
      INTERFACE set_bndy_fouier_m0
         MODULE PROCEDURE set_bndy_fouier_m0, set_bndy_fouier_m0_vec
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

      IF (UBOUND(gradientf,3) .ge. nsmax) THEN
         gradientf(:,:,nsmax:) = 0
      END IF

      IF (nsmin .eq. 1) THEN
         gradientf(:,:,1) = 0
      END IF

      IF (nsmax .eq. ns) THEN
         gradientf(:,:,ns) = gradientf(:,:,ns - 1)
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

      vech(:,:,nsmin + 1:nsmax) = (vecf(:,:,nsmin:nsmax - 1)                   &
                                +  vecf(:,:,nsmin + 1:nsmax))*0.5
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
      vech(:,:,nsmin + 1:nsmax) = (vecf(:,:,nsmin:nsmax - 1)                   &
                                +  vecf(:,:,nsmin + 1:nsmax))*0.5
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
         vecf(:,:,1) = vech(:,:,2)
      END IF
      IF (nsmax .eq. ns) THEN
         vecf(:,:,ns) = vech(:,:,ns)
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
      INTEGER                                                 :: fouruv
      INTEGER                                                 :: istat
      INTEGER                                                 :: nmin
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubumnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubvmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: da_uds
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: da_vds

!  Start of executable code
      CALL assert(nsmax .le. UBOUND(jbsupsmnh,3), 'UBOUND wrong in curl_ftoh')

      IF (iparity .EQ. f_sin) THEN
         sparity = 1
         fouruv = f_cos
      ELSE
         sparity = -1
         fouruv = f_sin
      END IF

      ALLOCATE(asubsmnh(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
               asubumnh(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
               asubvmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in curl_ftoh')

      CALL to_half_mesh(asubsmnf, asubsmnh)                       
      CALL to_half_mesh(asubumnf, asubumnh)
      CALL to_half_mesh(asubvmnf, asubvmnh)

      CALL set_bndy_fouier_m0(asubsmnh, asubumnh, asubvmnh, iparity)

!  Compute sqrt(g)B^X contravariant B-field components in Fourier space
      ALLOCATE(da_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                     &
               da_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation2 failed in curl_ftoh')

!  Radial gradients of asubmnu,v on half grid excluding origin and edge. Origin
!  and edge will be treated separately.
      CALL GradientHalf(da_uds, asubumnf)
      CALL GradientHalf(da_vds, asubvmnf)

      CALL set_bndy_fouier_m0(da_uds, fouruv)
      CALL set_bndy_fouier_m0(da_vds, fouruv)

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

      CALL set_bndy_fouier_m0(jbsupsmnh, jbsupumnh, jbsupvmnh, iparity)

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
!>  @param[in]    iparity   Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!>  @param[in]    extend    Extend the radial index of the curl.
!-------------------------------------------------------------------------------
      SUBROUTINE curl_htof(asubsmnh, asubumnh, asubvmnh,                       &
                           jbsupsmnf, jbsupumnf, jbsupvmnf, iparity,           &
                           nsmin, nsmax, extend)
      USE nscalingtools, ONLY: endglobrow

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: asubsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: asubumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: asubvmnh
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jbsupsmnf
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jbsupumnf
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jbsupvmnf
      INTEGER, INTENT(in)                                  :: iparity
      INTEGER, INTENT(in)                                  :: nsmin
      INTEGER, INTENT(in)                                  :: nsmax
      LOGICAL, INTENT(in)                                  :: extend

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
      INTEGER                                              :: fours
      INTEGER                                              :: fouruv

!  Start of executable code
      CALL assert(nsmin .ge. LBOUND(jbsupsmnf,3), 'LBOUND wrong in curl_htof')

      IF (iparity .eq. f_sin) THEN
         sparity = 1
         fours = f_sin
         fouruv = f_cos
      ELSE
         sparity = -1
         fours = f_cos
         fouruv = f_sin
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

!  FIXME: CHECK Boundary Conditions. asubs,subu,subv don't match older version.
      CALL set_bndy_half_to_full(asubsmnh, asubvmnh, nsmin, iparity,           &
                                 asubsmnf, asubumnf, asubvmnf)

!  Compute sqrt(g)B^X == A^X contravariant current components in Fourier space.
      ALLOCATE(da_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                          &
               da_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation2 failed in curl_htof')

!  Radial gradients of asubmnu,v on full grid excluding origin and edge. Origin
!  and edge will be treated separately.
      CALL GradientFull(da_uds, asubumnh)
      CALL GradientFull(da_vds, asubvmnh)

      CALL set_bndy_half_to_full_ds(asubumnh, da_uds, nsmin, fouruv)
      CALL set_bndy_half_to_full_ds(asubvmnh, da_vds, nsmin, fouruv)

      IF (extend) THEN
         nmax = MIN(ns, endglobrow + 1)
      ELSE
         nmax = MIN(ns, endglobrow)
      END IF
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

      CALL set_bndy_fouier_m0(jbsupsmnf, jbsupumnf, jbsupvmnf, iparity)
      IF (nsmin .eq. 1) THEN
         CALL set_bndy_full_origin(jbsupsmnf, jbsupumnf, jbsupvmnf)
      END IF

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
!>  @param[in]    iparity   Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!>  @param[in]    extend    Extend the radial index of the curl.
!-------------------------------------------------------------------------------
      SUBROUTINE curl_htof_targ(asubsmnh, asubumnh, asubvmnh,                  &
                                jbsupsmnf, jbsupumnf, jbsupvmnf,               &
                                iparity, nsmin, nsmax, extend)

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
      INTEGER, INTENT(IN)                                  :: iparity
      INTEGER, INTENT(IN)                                  :: nsmin
      INTEGER, INTENT(IN)                                  :: nsmax
      LOGICAL, INTENT(in)                                  :: extend

!  Local variables
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jbsp
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jbup
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jbvp

!  Start of executable code
      jbsp => jbsupsmnf
      jbup => jbsupumnf
      jbvp => jbsupvmnf
      CALL curl_htof(asubsmnh, asubumnh, asubvmnh,                             &
                     jbsp, jbup, jbvp, iparity, nsmin, nsmax, extend)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set fouier boundary conditions for full grid quantity converted from
!>         half grid.
!>
!>  At the center of the grid, the quantity should not be a function of u so all
!>  m > 0 modes should be zero. For the m = 0 terms, use the value from the half
!>  grid.
!>
!>  @param[in]    amnh    Half grid fourier quantity.
!>  @param[inout] amnf    Full grid fourier quantity.
!>  @param[in]    nsmin   Minimum radial index.
!>  @param[in]    iparity Parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_half_to_full(amnh, amnf, nsmin, iparity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in)    :: amnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: amnf
      INTEGER, INTENT(IN)                                     :: nsmin
      INTEGER, INTENT(in)                                     :: iparity

!  Start of executable code
      IF (nsmin .eq. 1) THEN
         amnf(:,:,1) = 0.0
         amnf(m0,:,1) = amnh(m0,:,2)
      END IF

      IF (iparity .eq. f_cos) THEN
         amnf(m0,-ntor:-n1,:) = 0.0
      ELSE
         amnf(m0,-ntor:n0,:) = 0.0
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set fouier boundary conditions for full grid vector quantity
!>         converted from half grid.
!>
!>  At the center of the grid, the a_u basis vector is zero. The a_s basis
!>  vector only has m=1 components. The a_v vector only has m=0 components.
!>
!>  @param[in]    asubsmnh Half grid covariant s fourier quantity.
!>  @param[in]    asubvmnh Half grid covariant v fourier quantity.
!>  @param[in]    nsmin    Minimum radial index.
!>  @param[in]    parity   Parity flag.
!>  @param[inout] asubsmnf Full grid covariant s fourier quantity.
!>  @param[inout] asubumnf Full grid covariant u fourier quantity.
!>  @param[inout] asubvmnf Full grid covariant v fourier quantity.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_half_to_full_vec(asubsmnh, asubvmnh,            &
                                                nsmin, parity,                 &
                                                asubsmnf, asubumnf, asubvmnf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in)    :: asubsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in)    :: asubvmnh
      INTEGER, INTENT(IN)                                     :: nsmin
      INTEGER, INTENT(in)                                     :: parity
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: asubsmnf
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: asubumnf
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: asubvmnf

!  Local Variables
      INTEGER                                                 :: moffh
      INTEGER                                                 :: mofff

!  Start of executable code
      moffh = LBOUND(asubsmnh, 1)
      mofff = LBOUND(asubsmnf, 1)

      IF (nsmin .eq. 1) THEN
         asubsmnf(:,:,1) = 0
         asubsmnf(m1 + mofff,:,1) = asubsmnh(m1 + moffh,:,2)

         asubumnf(:,:,1) = 0

         asubvmnf(:,:,1) = 0
         asubvmnf(m0 + mofff,:,1) = asubvmnh(m0 + moffh,:,2)
      END IF

      CALL set_bndy_fouier_m0(asubsmnf, asubumnf, asubvmnf, parity)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set fouier boundary conditions for full grid radial derivative
!>         quantity converted from half grid.
!>
!>  At the center of the grid the the m=1 quanity flips size across the s=0
!>  point. The results in a doubling of the m=1 component for the first half
!>  grid.
!>
!>  @param[in]    amnh   Half grid fourier quantity.
!>  @param[inout] amnf   Full grid fourier quantity.
!>  @param[in]    nsmin  Minimum radial index.
!>  @param[in]    parity Parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_half_to_full_ds(amnh, amnf_ds, nsmin, parity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: amnh
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnf_ds
      INTEGER, INTENT(IN)                        :: nsmin
      INTEGER, INTENT(in)                        :: parity

!  Local Variables
      INTEGER                                    :: moff

!  Start of executable code
      IF (nsmin .eq. 1) THEN
         moff = LBOUND(amnh, 1)

         amnf_ds(:,:,1) = 0.0
         amnf_ds(m1 + moff,:,1) = 2.0*ohs*amnh(m1 + moff,:,2)
      END IF

      CALL set_bndy_fouier_m0(amnf_ds, parity)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set fouier boundary conditions for full grid radial derivative
!>         quantity converted from half grid for uv vectors.
!>
!>  At the center of the grid, the quantity should not be a function of u so all
!>  m > 0 modes should be zero. For the m = 0 terms, 2*x/ds to preserve
!>  divergence free values.
!>
!>  @param[in]    aumnh  Half grid u fourier quantity.
!>  @param[in]    avmnh  Half grid v fourier quantity.
!>  @param[inout] aumnf  Full grid u fourier quantity.
!>  @param[inout] avmnf  Full grid u fourier quantity.
!>  @param[in]    nsmin  Minimum radial index.
!>  @param[in]    parity Parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_half_to_full_uv_ds(aumnh, avmnh,                &
                                                  aumnf_ds, avmnf_ds,          &
                                                  nsmin, parity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: aumnh
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: avmnh
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: aumnf_ds
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: avmnf_ds
      INTEGER, INTENT(IN)                        :: nsmin
      INTEGER, INTENT(in)                        :: parity

!  Local Variables
      INTEGER                                    :: moff

!  Start of executable code
      CALL set_bndy_half_to_full_ds(aumnh, aumnf_ds, nsmin, parity)
      CALL set_bndy_half_to_full_ds(avmnh, avmnf_ds, nsmin, parity)

      IF (nsmin .eq. 1) THEN
         moff = LBOUND(aumnh, 1)
         aumnf_ds(m0 + moff,:,1) = 2.0*ohs*aumnh(m0 + moff,:,2)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Apply origin vector boundary conditions.
!>
!>  At the center of the grid, a vector quantity should have only m=1 modes in
!>  the s direction, no modes in the u direction and only m=0 modes in the v
!>  direction.
!>
!>  @param[inout] amnsf  Full grid fourier quantity for s component.
!>  @param[inout] amnuf  Full grid fourier quantity for u component.
!>  @param[inout] amnvf  Full grid fourier quantity for v component.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_full_origin(amnsf, amnuf, amnvf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnsf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnuf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnvf

!  Local Variables
      INTEGER                                    :: moff

!  Start of executable code
      moff = LBOUND(amnsf, 1)

      amnsf(m0 + moff,:,1) = 0
      amnsf(m2 + moff:,:,1) = 0

      amnuf(:,:,1) = 0

      amnvf(m1 + moff:,:,1) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set vector fourier conditions.
!>
!>  For fouier modes for m=0 and for the negative n modes. The n=0 modes are
!>  zero for quantities with sin parity.
!>
!>  @param[inout] amn    Fourier quantity.
!>  @param[in]    parity Fouier parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_fouier_m0(amn, parity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amn
      INTEGER, INTENT(in)                        :: parity

!  Local Variables
      INTEGER                                    :: moff
      INTEGER                                    :: noff

!  Start of executable code
      moff = LBOUND(amn, 1)
      noff = LBOUND(amn, 2) + ntor

      IF (parity .eq. f_cos) THEN
         amn(m0 + moff,-ntor+noff:-n1+noff,:) = 0
      ELSE
         amn(m0 + moff,-ntor+noff:n0+noff,:) = 0
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set vector fourier conditions.
!>
!>  For fouier modes for m=0 and for the negative n modes. The n=0 modes are
!>  zero for quantities with sin parity.
!>
!>  @param[inout] amns   Fourier quantity for s component.
!>  @param[inout] amnu   Fourier quantity for u component.
!>  @param[inout] amnv   Fourier quantity for v component.
!>  @param[in]    parity Fouier parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_fouier_m0_vec(amns, amnu, amnv, parity)
      USE fourier, ONLY: f_cos, f_sin

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amns
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnu
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnv
      INTEGER, INTENT(in)                        :: parity

!  Start of executable code
      IF (parity .eq. f_cos) THEN
         CALL set_bndy_fouier_m0(amns, f_cos)
         CALL set_bndy_fouier_m0(amnu, f_sin)
         CALL set_bndy_fouier_m0(amnv, f_sin)
      ELSE
         CALL set_bndy_fouier_m0(amns, f_sin)
         CALL set_bndy_fouier_m0(amnu, f_cos)
         CALL set_bndy_fouier_m0(amnv, f_cos)
      END IF

      END SUBROUTINE

      END MODULE
