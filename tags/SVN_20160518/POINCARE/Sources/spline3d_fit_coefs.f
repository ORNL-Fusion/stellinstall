      MODULE spline3d_fit_coefs                                             ! Includes arrays used in 3D spline fits + geometric par
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE
      INTEGER, PARAMETER:: ks = 4, ku = 4, kv = 4                           ! KS,KU,KV: spline order
      INTEGER:: is, iphi, itht                                              ! Spline dimensions: IS(+KS), IPHI(+KU), ITHT(+KV)
      REAL(rprec):: smin, smax, thmin, thmax, phimin, phimax                ! Limits of spline box
      REAL(rprec), ALLOCATABLE, DIMENSION(:):: s_knots, u_knots,
     &   v_knots                                                            ! Knots of the spline in (s, u, v) space
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: br_coef,                 ! Spline coefficients for: BR, BZ, BPHI, R, Z, PHI
     &  bz_coef, bphi_coef, r_coef, z_coef, phi_coef
      LOGICAL, PARAMETER :: is_uniformx = .true.                            ! assume uniform spacing in (s, u, v) in supplied data f
      LOGICAL, PARAMETER :: is_uniformy = .true.
      LOGICAL, PARAMETER :: is_uniformz = .true.
      LOGICAL:: l_test_splines = .FALSE.                                    ! =T, computes maximum relative error at center of cell;

      CONTAINS

        SUBROUTINE ALLOC_SPLINES
        IMPLICIT NONE
        INTEGER:: istat

        ALLOCATE(s_knots(is+ks), u_knots(itht+ku),
     &    v_knots(iphi+kv), STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        s_knots = zero; u_knots = zero; v_knots = zero

        ALLOCATE(r_coef(is,itht,iphi), z_coef(is,itht,iphi),
     &    phi_coef(is,itht,iphi), br_coef(is,itht,iphi),
     &    bz_coef(is,itht,iphi), bphi_coef(is,itht,iphi),
     &    STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        r_coef = zero; z_coef = zero; phi_coef = zero
        br_coef = zero; bz_coef = zero; bphi_coef = zero

        END SUBROUTINE ALLOC_SPLINES

        SUBROUTINE DEALLOC_SPLINES
        IMPLICIT NONE

        DEALLOCATE(s_knots, u_knots, v_knots)
        DEALLOCATE(r_coef, z_coef, phi_coef, br_coef,
     &    bz_coef, bphi_coef)

        END SUBROUTINE DEALLOC_SPLINES

      END MODULE spline3d_fit_coefs
