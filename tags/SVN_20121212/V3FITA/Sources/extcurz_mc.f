!*******************************************************************************
!  File extcurz_mc.f
!  Contains module extcurz_mc
!           SUBROUTINE extcurz_mc_model_compute
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine extcurz_mc_model_compute
!  Both an extcurz_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE extcurz_mc
!    (extcurz - EXTernal CURrent along Z Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
! SECTION IV.    _OLD (SUPERCEDED) CODE
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

      MODULE extcurz_mc

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      USE stel_kinds
      USE v3f_global
      USE v3_utilities
      USE extcurz_T
      USE model_T

      IMPLICIT NONE

      CONTAINS
          
!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************

      SUBROUTINE extcurz_mc_model_compute(a_extcurz,a_model,                   &
     &   mod_signal,mod_sigma)

      USE stel_constants, ONLY: twopi, mu0
      USE line_segment, ONLY: line_seg ! linear interpolation

      TYPE (extcurz_desc), INTENT (inout) :: a_extcurz
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:) :: mod_signal, mod_sigma

      CHARACTER(*), PARAMETER :: sub_name = 'extcurz_mc_model_compute'
      INTEGER :: ns, jj, istat
      REAL(rprec) :: ds, bsubv00c, bsubv00s
      REAL(rprec), ALLOCATABLE :: svec(:), bvec(:)

!-------------------------------------------------------------------------------
!  START OF EXECUTABLE CODE
!-------------------------------------------------------------------------------         

      !WRITE(*,*) ' *** IN SUBROUTINE ', sub_name
      !WRITE(*,*) ' extcurz diag:  in =', a_extcurz % s0, a_extcurz % u0

      ! allocate data and sigma
      IF (ASSOCIATED(mod_signal)) THEN
         DEALLOCATE(mod_signal,mod_sigma, stat=istat)
         CALL assert_eq(0,istat,sub_name //                                   &
     &                     'mod_signal, sigma dealloc')
      ENDIF
      ALLOCATE(mod_signal(1),mod_sigma(1),stat=istat)
      CALL assert_eq(0,istat,sub_name // 'mod_signal, sigma alloc')

      ! get ns and compute ds
      ns = a_model % eqstate % aux1 % ns
      ds = 1. / (ns-1)

      ! allocate half grid vectors
      ALLOCATE(svec(2:ns),bvec(2:ns),stat=istat)
      CALL assert_eq(0,istat,sub_name // 'svec, bvec alloc')

      ! store s values
      svec(2) = ds/2.
      DO jj=3, ns
         svec(jj) = svec(jj-1) + ds
      END DO

      ! compute cosine component bsubv00c
      bvec(2:ns) = a_model % eqstate % aux1 % bsubvmnc(1,2:ns)
      CALL line_seg(  x=a_extcurz % s0,  y=bsubv00c,                           &
     &               xx=svec(2:ns),     yy=bvec(2:ns), n=ns-1 )

      IF ( a_model % eqstate % fixp % lasym ) THEN
         ! compute sine component bsubv00s
         bvec(2:ns) = a_model % eqstate % aux1 % bsubvmns(1,2:ns)
         CALL line_seg(  x=a_extcurz % s0,  y=bsubv00s,                        &
     &                  xx=svec(2:ns),     yy=bvec(2:ns), n=ns-1 )
      ELSE
         bsubv00s = 0.
      END IF

      ! free memory
      DEALLOCATE(svec,bvec,stat=istat)
      CALL assert_eq(0,istat,sub_name // 'svec, bvec dealloc')

      ! compute model: current = circuitation / mu0
      mod_signal(1) = twopi * ( bsubv00c + bsubv00s ) / mu0
      mod_sigma (1) = 0.

      !WRITE(*,*) '               out =', mod_signal(1), mod_sigma(1)

      END SUBROUTINE extcurz_mc_model_compute

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  GLT 10-sep-2012
!     first draft by Gregorio L. Trevisan <gtrevisan@igi.cnr.it>
!
!  GLT 19-sep-2012
!     final version
!

      END MODULE extcurz_mc

