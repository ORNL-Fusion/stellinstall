      MODULE lsode_quantities
      USE stel_kinds, ONLY: rprec
      USE stel_constants
      USE spline3d_fit_coefs
      IMPLICIT NONE
      INTEGER:: num_eqns, norbs, lrw, liw, itol, itask, iopt,
     &  istate, jt
      INTEGER, ALLOCATABLE, DIMENSION(:):: iwork
      REAL(rprec), ALLOCATABLE, DIMENSION(:):: y, f, rwork
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: s_nb,
     &  phi_nb, theta_nb, bz_s_nb, br_s_nb, bphi_s_nb, z_s_nb,
     &  r_s_nb, drds_s_nb, dzds_s_nb, drdphi_s_nb, dzdtht_s_nb,
     &  dzdphi_s_nb, drdtht_s_nb

      CONTAINS

        SUBROUTINE ALLOC_LSODE
        IMPLICIT NONE
        INTEGER:: istat

        ALLOCATE(y(num_eqns), f(num_eqns),  STAT=istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        lrw = 52+9*num_eqns+(num_eqns**2); liw = 50+num_eqns
        ALLOCATE(iwork(liw), rwork(lrw), STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        y = zero; f = zero; rwork = zero; iwork = 0

        ALLOCATE(s_nb(norbs), theta_nb(norbs),
     &    phi_nb(norbs), bz_s_nb(norbs), br_s_nb(norbs),
     &    bphi_s_nb(norbs), r_s_nb(norbs), z_s_nb(norbs),
     &    drds_s_nb(norbs), dzds_s_nb(norbs),
     &    drdtht_s_nb(norbs), dzdtht_s_nb(norbs), drdphi_s_nb(norbs),
     &    dzdphi_s_nb(norbs), STAT = istat)
        IF (istat .NE. 0) STOP 'Allocation error'
        s_nb = zero; theta_nb = zero; phi_nb = zero
        br_s_nb = zero; bz_s_nb = zero; bphi_s_nb = zero
        r_s_nb = zero; drds_s_nb = zero; dzds_s_nb = zero
        drdtht_s_nb = zero; dzdtht_s_nb = zero; z_s_nb = zero
        drdphi_s_nb = zero; dzdphi_s_nb = zero

        END SUBROUTINE ALLOC_LSODE

        SUBROUTINE DEALLOC_LSODE
        IMPLICIT NONE

        DEALLOCATE(y, f, iwork, rwork)
        DEALLOCATE(s_nb, theta_nb, phi_nb, bz_s_nb, br_s_nb,
     &    bphi_s_nb, r_s_nb, drds_s_nb, dzds_s_nb, drdtht_s_nb,
     &    dzdtht_s_nb, drdphi_s_nb, dzdphi_s_nb, z_s_nb)

        END SUBROUTINE DEALLOC_LSODE

      END MODULE lsode_quantities
