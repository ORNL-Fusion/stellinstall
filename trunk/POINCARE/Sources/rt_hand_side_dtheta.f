

      SUBROUTINE rt_hand_side_dtheta(neq, theta, yin, fout)                     ! Defines field-line equations for LSODE integrator for type="poloidal"
      USE stel_kinds                                                            ! The right hand sides are B^s/B^u  and B^v/B^u.
      USE stel_constants                                                        ! To get them, one simply uses:
      USE bspline                                                               ! B^R =B*grad R = B^sRs + B^uRu + B^vRv
      USE spline3d_fit_coefs                                                    ! B^z =B*grad Z = B^sZs + B^uZu + B^vZu
      USE lsode_quantities                                                      ! B^Phi = B*grad PHi = B^sPhi_s + B^u Phi_u + B^vPhi
      IMPLICIT NONE                                                             ! ... and inverts to get B^s and B^u.
      INTEGER, INTENT(IN):: neq                                                 ! Recall that it is assumed that Phi=v, thus Phi_s = Phi_u = 0
      REAL(KIND=rprec), INTENT(IN):: theta
      REAL(KIND=rprec), DIMENSION(neq), INTENT(IN) :: yin
      REAL(KIND=rprec), DIMENSION(neq), INTENT(OUT):: fout
      INTEGER:: i, isize
      REAL(KIND=rprec) :: s_eval, phi_eval, theta_eval,
     &  bs, bu, bv, jac

      isize = neq/2
      DO i = 1, isize                                                           !loop over field lines
        s_eval = yin(2*i-1); theta_eval = theta
        phi_eval = yin(2*i)
        CALL reset_coords(s_eval, theta_eval, phi_eval)
        s_nb(i) = s_eval; phi_nb(i) = phi_eval
        theta_nb(i) = theta_eval
      ENDDO

      CALL dbs3vd(0, 0, 0, isize, s_nb, theta_nb, phi_nb, ks, ku,               ! Br
     &  kv, s_knots, u_knots, v_knots, is, itht, iphi, br_coef,
     &  br_s_nb, is_uniformx, is_uniformy, is_uniformz )

      CALL dbs3vd(0, 0, 0, isize, s_nb, theta_nb, phi_nb, ks, ku,               ! Bphi
     &  kv, s_knots, u_knots, v_knots, is, itht, iphi, bphi_coef,
     &  bphi_s_nb, is_uniformx, is_uniformy, is_uniformz )

      CALL dbs3vd(0, 0, 0, isize, s_nb, theta_nb, phi_nb, ks, ku,               ! Bz
     &   kv, s_knots, u_knots, v_knots, is, itht, iphi, bz_coef,
     &   bz_s_nb, is_uniformx, is_uniformy, is_uniformz )

      CALL dbs3vd(0, 0, 0, isize, s_nb, theta_nb, phi_nb, ks, ku,               ! R
     &  kv, s_knots, u_knots, v_knots, is, itht, iphi, r_coef,
     &  r_s_nb, is_uniformx, is_uniformy, is_uniformz )

      DO i = 1, isize
        s_eval = s_nb(i); theta_eval = theta_nb(i)
        phi_eval = phi_nb(i)
        drds_s_nb(i) = dbs3dr(1, 0, 0, s_eval, theta_eval, phi_eval,            ! dR/ds
     &    ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi, r_coef)
        drdtht_s_nb(i) = dbs3dr(0, 1, 0, s_eval, theta_eval, phi_eval,          ! dR/du
     &    ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi, r_coef)
        drdphi_s_nb(i) = dbs3dr(0, 0, 1, s_eval, theta_eval, phi_eval,          ! dR/dv
     &    ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi, r_coef)
        dzds_s_nb(i) = dbs3dr(1, 0, 0, s_eval, theta_eval, phi_eval,            ! dZ/ds
     &    ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi, z_coef)
        dzdtht_s_nb(i) = dbs3dr(0, 1, 0, s_eval, theta_eval, phi_eval,          ! dZ/du
     &    ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi, z_coef)
        dzdphi_s_nb(i) = dbs3dr(0, 0, 1, s_eval, theta_eval, phi_eval,          ! dZ/dv
     &    ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi, z_coef)
      ENDDO

      DO i = 1, isize                                                           ! Field line ODEs in cylindrical geometry
        jac = drds_s_nb(i)*dzdtht_s_nb(i) - drdtht_s_nb(i)*dzds_s_nb(i)
        bs = r_s_nb(i)*(dzdtht_s_nb(i)*br_s_nb(i) -                             !B^s = (R(Zu*B^R-RuB^Z))+B^Phi*(RuZv-RvZu))/Jacob
     &    drdtht_s_nb(i)*bz_s_nb(i)) + bphi_s_nb(i)*(
     &    drdtht_s_nb(i)*dzdphi_s_nb(i) - drdphi_s_nb(i)*dzdtht_s_nb(i))
        bu = r_s_nb(i)*(drds_s_nb(i)*bz_s_nb(i) -                               !B^u = (R(RsB^z-ZsB^R)}B^Phi*(ZsRv-RsZu))/Jacob
     &    dzds_s_nb(i)*br_s_nb(i)) + bphi_s_nb(i)*(
     &    drdphi_s_nb(i)*dzds_s_nb(i) - drds_s_nb(i)*dzdphi_s_nb(i))
        bv = bphi_s_nb(i)
        IF (bu == zero) STOP 'BU becomes zero for some orbit!'
        fout(2*i-1) = bs/bu                                                     ! r.h.s. = B^s/B^u
        fout(2*i) = (jac*bv)/bu                                                 ! r.h.s. = B^v/B^u
      END DO

      END SUBROUTINE rt_hand_side_dtheta
