      SUBROUTINE reset_coords(s, theta, phi)                                 ! Keeps orbit within volume where spline fits are valid
      USE stel_kinds
      USE stel_constants
      USE spline3d_fit_coefs
      IMPLICIT NONE
      REAL(rprec), INTENT(INOUT) :: s, theta, phi

      IF(s < smin) s = smin
      IF(s > smax) s = smax

      IF(theta < zero) THEN
        DO
          theta = theta + twopi
          IF(theta > zero .AND. theta <= twopi) EXIT
        ENDDO
      ELSEIF(theta > twopi) THEN
        DO
          theta = theta - twopi
          IF(theta > zero .AND. theta <= twopi) EXIT
        ENDDO
      ENDIF

      IF(phi < zero) THEN
        DO
          phi = phi + twopi
          IF(phi > zero .AND. phi <= twopi) EXIT
        ENDDO
      ELSEIF(phi > twopi) THEN
        DO
          phi = phi - twopi
          IF(phi > zero .AND. phi <= twopi) EXIT
        ENDDO
      ENDIF

      END SUBROUTINE reset_coords
