      SUBROUTINE enforce_periodicity(ns, nu, nv, array0, i_phi)
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE
      INTEGER, INTENT(IN):: ns, nu, nv, i_phi
      REAL(kind=rprec), DIMENSION(ns, nu, nv), INTENT(INOUT):: array0

      IF (i_phi == 0) THEN
        array0(:,nu,nv) = array0(:,1,1)                                    ! First make sure that four corners are identical
        array0(:,1,nv) = array0(:,1,1)
        array0(:,nu,1) = array0(:,1,1)
        array0(:,nu,2:nv-1) = array0(:,1,2:nv-1)                           ! then, do full lines
        array0(:,2:nu-1,nv) = array0(:,2:nu-1,1)
      ELSE                                                                 ! PHI is not periodic on v!!
        array0(:,1:nu,nv) = twopi
        array0(:,nu,1:nv-1) = array0(:,1,1:nv-1)
      ENDIF

      END SUBROUTINE enforce_periodicity
