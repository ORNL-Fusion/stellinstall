      MODULE bivariate
      USE stel_kinds

      TYPE bivariate_type
         PRIVATE
         INTEGER :: nrgrid, nzgrid, ndim
         INTEGER, DIMENSION(:,:), ALLOCATABLE :: irz11_bi, irz12_bi
         INTEGER, DIMENSION(:,:), ALLOCATABLE :: irz21_bi, irz22_bi
         REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: w11_bi, w12_bi
         REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: w21_bi, w22_bi
      END TYPE

C-----------------------------------------------
!
!     Written by S. P. Hirshman (01/30/03)
!
!  Modified 12-20-2004 JDH
!    Changed logic so that if module arrays are already allocated, they are
!    first deallocated, and then allocated to the correct size. (Too lazy to 
!    worry about checking to see if they are already allocated to the
!    appropriate size.)
!
!     Local variables:
!
!     irz11_bi(ns*nu,kp) : lower left corner indices on a
!                          rectangular rgrid,zgrid enclosing 
!                          the point (s,u,v)
!     irz22_bi(ns*nu,kp) : upper right corner indices
!     wij_bi(ns*nu,kp)  : weight factors for 4 pt interpolation
!
C-----------------------------------------------
     
      CONTAINS

      SUBROUTINE setbivariate (this, rsu, zsu,                                 &
     &                         rgrid, zgrid, kin, kp, ib)
      IMPLICIT NONE
      TYPE (bivariate_type), INTENT(inout) :: this
      INTEGER, INTENT(in) :: kin, kp, ib
      REAL(rprec), INTENT(in) :: rsu(:,:), zsu(:,:), rgrid(:), zgrid(:)
      INTEGER :: ku, js, istat, ns, nu, index1d, ir, jz, ir1, jz1
      REAL(rprec) :: rad0, zee0, ri, zj, pr, qz, temp, delr, delz
C-----------------------------------------------
!
!     Compute the indices (ir_bi, jz_bi) and weight factors (wij_bi) 
!     for performing bivariate (4 pt) interpolation FROM a rectangular (R X Z)
!     grid TO a general (non-orthogonal) grid (s, u)
!
!     rgrid(1:nr)     :  r coordinate on a fixed, equally spaced grid 
!     zgrid(1:nz)     :  z coordinate on a fixed, equally spaced grid
!                        (may be different mesh size than rgrid)
!     rsu(1:ns,1:nu)  :  r coordinate at s,u
!     zsu(1:ns,1:nu)  :  z coordinate at s,u
!
!     kin             :  toroidal plane index (1 <= kin <= kp)
!     kp              :  number of toroidal planes
!     ib              :  starting radial s index
!
!     index1d         :  index for making indices/weight arrays 1D (in each toroidal plane)
!                        (improves speed by as much as 2X)
!     
C-----------------------------------------------
      this%nrgrid = SIZE(rgrid); this%nzgrid = SIZE(zgrid)
      ns = SIZE(rsu,1);  nu = SIZE(rsu,2)
      this%ndim = (ns-ib+1)*nu
      istat = 0

      delr = rgrid(2) - rgrid(1)
      delz = zgrid(2) - zgrid(1)

!  JDH 12-20-2004 Change Allocation - Deallocation logic.
      IF (ALLOCATED(this%irz11_bi)) THEN
         DEALLOCATE (this%irz11_bi, this%irz12_bi, this%irz21_bi,              &
     &               this%irz22_bi, this%w11_bi, this%w12_bi,                  &
     &               this%w21_bi, this%w22_bi, stat = istat)
         IF (istat .ne. 0) STOP 'Deallocation error in setbivariate'
      ENDIF

      ALLOCATE (this%irz11_bi(this%ndim,kp),                                   &
     &          this%irz12_bi(this%ndim,kp),                                   &
     &          this%irz21_bi(this%ndim,kp),                                   &
     &          this%irz22_bi(this%ndim,kp),                                   &
     &          this%w11_bi(this%ndim,kp), this%w12_bi(this%ndim,kp),          &
     &          this%w21_bi(this%ndim,kp), this%w22_bi(this%ndim,kp),          &
     &          stat = istat)
      IF (istat .ne. 0) STOP 'Allocation error in setbivariate'

!      IF (.not.ALLOCATED(irz11_bi)) ALLOCATE (irz11_bi(ndim,kp), 
!     &    irz12_bi(ndim,kp), irz21_bi(ndim,kp), irz22_bi(ndim,kp),
!     &    w11_bi(ndim,kp), w12_bi(ndim,kp), 
!     &    w21_bi(ndim,kp), w22_bi(ndim,kp),stat = istat)
!      IF (istat .ne. 0) STOP 'Allocation error in setbivariate'

      index1d = 0

      DO ku = 1, nu
      DO js = ib, ns

         index1d = index1d + 1
!
!       CHECK THAT BOUNDARY POINTS ARE INSIDE GRID.  IF NOT, STOP!
!
         rad0 = rsu(js, ku)
         zee0 = zsu(js, ku)

         IF (rad0.lt.rgrid(1) .or. rad0.gt.rgrid(this%nrgrid) .or.
     &       zee0.lt.zgrid(1) .or. zee0.gt.zgrid(this%nzgrid)) 
     &       STOP 'Plasma point outside response function grid!'
!
!       DETERMINE INTEGER INDICES (IR,JZ) FOR LOWER LEFT R, Z CORNER GRID POINT
!
         ir = INT((rad0 - rgrid(1))/delr) + 1
         jz = INT((zee0 - zgrid(1))/delz) + 1
         ir1 = MIN(this%nrgrid,ir+1)
         jz1 = MIN(this%nzgrid,jz+1)

!
!        STORE INDICES IN 1D ARRAYS
!
         this%irz11_bi(index1d,kin) = ir  + this%nrgrid*(jz-1)
         this%irz22_bi(index1d,kin) = ir1 + this%nrgrid*(jz1-1)
         this%irz12_bi(index1d,kin) = ir  + this%nrgrid*(jz1-1)
         this%irz21_bi(index1d,kin) = ir1 + this%nrgrid*(jz-1)
!
!       COMPUTE RI, ZJ AND PR , QZ AT GRID POINT (IR , JZ)
!       ALSO, COMPUTE WEIGHTS WIJ FOR 4 CORNER GRID POINTS
!
         ri = rgrid(ir)
         zj = zgrid(jz)
         pr = (rad0 - ri)/delr
         qz = (zee0 - zj)/delz
         temp = pr*qz
         this%w22_bi(index1d,kin) = temp                   !p*q               
         this%w21_bi(index1d,kin) = pr - temp              !p*(1-q)
         this%w12_bi(index1d,kin) = qz - temp              !q*(1-p)
         this%w11_bi(index1d,kin) = 1 + temp - (pr + qz)   !(1-p)*(1-q)

      END DO
      END DO

      END SUBROUTINE setbivariate

      SUBROUTINE bivariate4pt (this, resp_rz, resp_su, kin)
      IMPLICIT NONE
      TYPE (bivariate_type), INTENT(inout) :: this
      INTEGER, INTENT(in)      :: kin
      REAL(rprec), INTENT(in)  :: resp_rz(this%nrgrid*this%nzgrid)
      REAL(rprec), INTENT(out) :: resp_su(this%ndim)
      INTEGER :: jsu, irz11, irz12, irz21, irz22
!-----------------------------------------------
!
!     DETERMINE resp_uv ON THE S-U GRID BY 4 PT BIVARIATE INTERPOLATION
!     OF resp_rz WHICH IS ON AN R-Z MESH, ALL AT A FIXED TOROIDAL PLANE
!
!     USES 2-D INTERPOLATION BASED ON THE FOUR POINT FORMULA
!     IN ABRAMOWITZ AND STEGUN, EQ. 25.2.66 
!
!     resp_rz(1:nrgrid, 1:nzgrid) :   response on the rectangular R, Z grid
!
!     resp_su(1:ns,1:nu)          :   interpolated response on s,u grid (as 1d array)
!
!-----------------------------------------------
      DO jsu = 1, this%ndim
         
         irz11 = this%irz11_bi(jsu,kin)
         irz22 = this%irz22_bi(jsu,kin)
         irz12 = this%irz12_bi(jsu,kin)
         irz21 = this%irz21_bi(jsu,kin)
!
!       COMPUTE RESPONSE AT R, PHI (fixed kin), Z BY 4-PT INTERPOLATION
!
         resp_su(jsu) = this%w11_bi(jsu,kin)*resp_rz(irz11) 
     &                + this%w22_bi(jsu,kin)*resp_rz(irz22) 
     &                + this%w21_bi(jsu,kin)*resp_rz(irz21) 
     &                + this%w12_bi(jsu,kin)*resp_rz(irz12)

      END DO

      END SUBROUTINE bivariate4pt

      SUBROUTINE clear_bivar(this)

      IMPLICIT NONE
      TYPE (bivariate_type), INTENT(inout) :: this

      IF (ALLOCATED(this%irz11_bi)) DEALLOCATE (this%irz11_bi)
      IF (ALLOCATED(this%irz12_bi)) DEALLOCATE (this%irz12_bi)
      IF (ALLOCATED(this%irz21_bi)) DEALLOCATE (this%irz21_bi)
      IF (ALLOCATED(this%irz22_bi)) DEALLOCATE (this%irz22_bi)
      IF (ALLOCATED(this%w11_bi)) DEALLOCATE (this%w11_bi)
      IF (ALLOCATED(this%w12_bi)) DEALLOCATE (this%w12_bi)
      IF (ALLOCATED(this%w21_bi)) DEALLOCATE (this%w21_bi)
      IF (ALLOCATED(this%w22_bi)) DEALLOCATE (this%w22_bi)
  
      END SUBROUTINE clear_bivar

      END MODULE bivariate
