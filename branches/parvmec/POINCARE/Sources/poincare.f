      PROGRAM poincare
!
!=====================================================================================================================
!
!     Program = Poincare
!     Version = 1.10         03/30/10
!               1.20         05/10/13  (add mpi)
!
!     Authors: Don A. Spong (original) and Raul Sanchez (updated), Fusion Energy Division, Oak Ridge National Laboratory
!
!     For problems, please contact: sanchezferlr@ornl.gov or spongda@ornl.gov
!
!     Description:
!
!       This program computes Poincare plots for a 3D magnetic field provided using flux-coordinates (s, theta, phi).
!       It is assumed that s in [0,1], u(=pol) in [0,2pi) and v(=tor) in [0,2pi). Values at u=2pi or v=2pi should
!       not be provided, since they are automatically set to the values at u=0 and v=0 to ensure periodicity
!       inside the program. s and u may be arbitrary flux coordinates, but it is assumed that v= geom. PHI, as in VMEC.
!
!       Poincare plots can be done on poloidal or toroidal cross-section, depending on input TYPE on command line,
!       which can be set to "poloidal" (theta is the independent variable) or "toroidal" (phi is the independent variable).
!       Choice should be done depending on whether Bu or Bv go through a zero. I.e., if BV is nonzero, choose "toroidal".
!       Angle for cross section can also be chosen on command line via the ANGLE input variable.
!
! ** INPUT:The field information must be provided in an ASCII file whose first line contains the size of the 3D grid:
!              NV = No. of toroidal points;   NS = No. of radial points;   NU = No. of poloidal points
!       Then, data must follow in six colums, giving the cylindrical coordinates info:
!              R(s,u,v)    Z(s,u,v)    Phi(s,u,v)   Br(s,u,v)   Bz(s,u,v)   BPhi(s,u,v)
!       The order assumed is that u is the fastest running, and v is the slowest running coordinate in the file.
!
! ** EXECUTION: The program is run from the command line invoking:
!              xpoincare   FILENAME.dat    N.Toroidal transits    N.surfaces     N.orbits/surface    N.steps/transit
!                          Innermost Surface     Outermost Surface     TYPE    ANGLE
!
!        Here, FILENAME.dat = name of file containing 3D field info (like the one produced by SIESTA code),
!              N. Toroidal transits = No. of 2*pi toroidal periods that must be followed per orbit
!              N. surfaces = No. of different s-values where orbits will be started
!              N. orbits/surface = No. of orbits to initialize, uniformly distributed, on each surface
!              N. steps/transit = the numerical integrator will integrate each transits in steps of size =  2*pi/(N.steps/transit)
!              Surfaces will be distributed between [Innermost Surface, Outermost Surface]
!              TYPE = "poloidal" (theta is ind.var.; choose if Bv nonzero) or "toroidal" (phi is ind. var.; choose if Bu nonzero)
!              ANGLE = angular cross-section where puncture plot computed. It is poloidal or toroidal depending on TYPE.
!
! ** OUTPUT: The output is a 4-column ASCII file containing a list of the following values:
!             Theta    s      R     Z
!             of all the crossings of the orbits with the plane PHI = ANGLE if TYPE = "toroidal".
!             Phi    s      R    Phi 
!             of all the crossings of the orbits with the plane THETA = ANGLE if TYPE = "poloidal".
!             The name of the file is simply FILENAME_puncture_TYPE.dat
!
! ** HISTORY:
!
!   10/24/09: First version released, based on a prior version by Don Spong.                            (Version: 1.00)
!   03/30/10: Several minor bugs fixed.                                                                 (Version: 1.10)
!
!=====================================================================================================================
!
      USE stel_kinds
      USE stel_constants
      USE bspline                                                           ! Module in LIBSPLINE
      USE spline3d_fit_coefs
      USE lsode_quantities
      USE mpi_params
#if defined(MPI_OPT)
      USE mpi_inc
#endif
      IMPLICIT NONE

      INTEGER:: i, j, k, iargc, numargs, isize, it,
     &  steps, num_orbs, icount, nsurfs, mn, scndry_steps
      REAL(rprec)::  s_test, theta_test, phi_test, delta,
     &  ang_in, ang_out, rdum, tstart, tend, s1, s2
      REAL(rprec):: angle, phi_puncture, theta_puncture
      REAL(rprec), PARAMETER :: relerr=1.d-8, abserr=1.d-8,
     &  eps=1.d-10
      REAL(rprec),ALLOCATABLE, DIMENSION(:,:):: start_points
      REAL(rprec),ALLOCATABLE, DIMENSION(:):: theta_nb_tot,
     &  phi_nb_tot, s_nb_tot, r_s_nb_tot, z_s_nb_tot
      LOGICAL:: l_ex
      INTEGER:: i_truncate, ilen, iargs
      INTEGER :: local_num_orbs, norbs_tot, istat
      CHARACTER(LEN=8):: type_in                                               ! ="poloidal", THETA is ind. variable; = "toroidal", PHI is ind. variable
      CHARACTER(LEN=100):: outfile
      CHARACTER(LEN=60):: arg1, arg2, arg3, arg4, arg5,
     &  arg6, arg7, arg8, arg9, arg10
      CHARACTER(LEN=13), PARAMETER:: banner = 'Version 1.20',
     &  update = 'May 10, 2013'
      EXTERNAL rt_hand_side_dphi, rt_hand_side_dtheta, jacobian             ! Used by the integrator LSODE

      CALL cpu_time(tstart)

C
C Get NPES and myid.  Requires initialization of MPI.
C
#if defined(MPI_OPT)
      CALL MPI_INIT(ierr_mpi)
      IF (ierr_mpi .NE. 0) THEN
        WRITE(*,'("MPI_INIT returned IER =")') ierr_mpi
        STOP
        ENDIF
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr_mpi)
      IF (ierr_mpi .NE. 0) THEN
        WRITE(*,'("MPI_COMM_SIZE returned IER =")') ierr_mpi
        STOP
        ENDIF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr_mpi)
      IF (ierr_mpi .NE. 0) THEN
        WRITE(*,'("MPI_COMM_RANK returned IER =")') ierr_mpi
        STOP
        ENDIF
#else
      numprocs = 1
      myid = 0
#endif
C

      IF(myid .eq. 0) THEN
      WRITE(*,*) '=================================================='
     & , '=================================================='
      WRITE(*,*) ' This the POINCARE puncture code '
      WRITE(*,*) ' ', TRIM(banner), '          Last update: ',
     &  TRIM(update)
      WRITE(*,'("  Runnning on ",i5," processors")') numprocs
      WRITE(*,*)
      WRITE(*,*) ' Please report problems to: ',
     &   ' Don Spong (spongda@ornl.gov)'
     & , ' or Raul Sanchez (sanchezferlr@ornl.gov).'
      WRITE(*,*) '--------------------------------------------------'
     &  , '--------------------------------------------------'
      ENDIF      !IF(myid .eq. 0) THEN

!     Default values for command line arguments:

      arg1 = "bfield.dat"
      steps = 10; num_orbs = 10; scndry_steps = 50
      nsurfs = 10; s1 = eps; s2 = 1.d0 -eps
      angle = 0.d0
!
      numargs = iargc()         !for non-Cray X1 platforms; for Cray X1 use: numargs = ipxfargc()

      IF( numargs .lt. 1 .OR. TRIM(arg1) == '-h')THEN
      IF(myid .eq. 0) THEN
        WRITE(*,*)
        PRINT *,'   Run with arguments using:'
        WRITE(*,*)
        WRITE(*,*) '    xpoincare  file  transits  surfaces  orbts/surf'
     &   , '  steps/period  inner-surf(0->1)  outer-surf(0->1) '
     &   , '  type    cross-section (0->360)'
        WRITE(*,*)
        PRINT *,'   - FILE is the BFIELD file produced by SIESTA'
        PRINT *,'   - transits are the number orbit transits (tor/pol)'
        PRINT *,'   - steps/period are the steps used by integrator'
        PRINT *,'   - inner-surf, outer-surf must be between 0 and 1.'
        PRINT *,'   - type = poloidal or toroidal, to choose ind. var.'
        PRINT *,'   - choose type = toroidal if BPHI nonnegative.'
        PRINT *,'   - choose type = poloidal if BPHI =0 somewhere,'
        PRINT *,'      and BTHETA always nonnegative.'
        PRINT *,'   - cross-section = angle (deg) where puncture done.'
        WRITE(*,*)
        WRITE(*,*) ' Or run with arguments:'
        WRITE(*,*)
        WRITE(*,*) '    xpoincare  file  test_splines'
        WRITE(*,*)
        WRITE(*,*) ' to test splines quality with data from FILE.'
        WRITE(*,*)
        WRITE(*,*) '=================================================='
     &    , '=================================================='
      ENDIF      !IF(myid .eq. 0) THEN
        STOP
      ELSE IF (numargs == 2 .AND. TRIM(arg2) == "test_splines") THEN
        l_test_splines = .TRUE.
	ELSE
         DO iargs = 1, numargs
            IF (iargs == 1) CALL getarg(1, arg1)
            IF (iargs == 2) CALL getarg(2, arg2)
            IF (iargs == 3) CALL getarg(3, arg3)
            IF (iargs == 4) CALL getarg(4, arg4)
            IF (iargs == 5) CALL getarg(5, arg5)
            IF (iargs == 6) CALL getarg(6, arg6)
            IF (iargs == 7) CALL getarg(7, arg7)
            IF (iargs == 8) CALL getarg(8, arg8)
            IF (iargs == 9) CALL getarg(9, arg9)
         END DO
         DO iargs = numargs+1,9
            IF (iargs == 2) arg2 = "100"                    !toroidal periods (50-500)
            IF (iargs == 3) arg3 = "50"                     !nsurfs
            IF (iargs == 4) arg4 = "12"                     !num_orbs/surf
            IF (iargs == 5) arg5 = "40"                     !scndry_steps
            IF (iargs == 6) arg6 = "0.05"   !"zoom: 0.65"
            IF (iargs == 7) arg7 = "0.95"   !"zoom: 0.85"
            IF (iargs == 8) arg8 = "toroidal"
            IF (iargs == 9) arg9 = "0"
         END DO
      ENDIF
            !for non-Cray X1 platforms. Read spline fit data, allocate arrays, open output files

      read(arg2,'(i6)') steps
      read(arg3,'(i6)') nsurfs
      read(arg4,'(i6)') num_orbs
      read(arg5,'(i6)') scndry_steps
      read(arg6,'(f8.4)') s1
      read(arg7,'(f8.4)') s2
      IF (s1 > s2) THEN
        s_test = s1
        s1 = s2
        s2 = s1
      ENDIF
      IF (s1 == s2) nsurfs = 1
      type_in = TRIM(arg8)                                                  ! Type = poloidal or toroidal
	CALL tolower(type_in)                                                 ! LIBSTELL function
      read(arg9, *) angle
      angle = angle*(twopi/360.d0)                                          ! Move to radians
      angle = MOD(angle, twopi)                                             ! Make sure angle between 0. and 2*pi
      IF (angle < 0.d0) angle = twopi + angle                               ! Make sure is positive
      IF (type_in == 'poloidal') THEN
        theta_puncture = angle
       IF(myid .eq. 0) THEN
        WRITE(*,'(a,1p,e10.4,a)') 
     &  ' Poloidal cross section done at angle: ', angle,' rads.'
       ENDIF      !IF(myid .eq. 0) THEN
      ELSE IF (type_in == 'toroidal') THEN
        phi_puncture = angle
       IF(myid .eq. 0) THEN
        WRITE(*,'(2x,a,1p,e10.4,a)') 
     &  'Toroidal cross section done at angle: ', angle,' rads.'
       ENDIF      !IF(myid .eq. 0) THEN
      ELSE
        STOP 'Wrong type!'
      ENDIF
!
      CALL spline_fields(arg1)                                              ! Compute splines. Coefficients and knots in SPLINE modu
      IF (l_test_splines) THEN                                              ! When testing splines, do not run integration. Exit.
       IF(myid .eq. 0) THEN
        WRITE(*,*)
        WRITE(*,*) '  Splines tested.'
        WRITE(*,*)
        CALL CPU_TIME(tend)
        WRITE(*,'(a15, f12.5, a5)') '  Total time = ', tend - tstart,
     &    ' sec.'

        WRITE(*,*) '=================================================='
     &    , '=================================================='
      ENDIF      !IF(myid .eq. 0) THEN
        STOP

      ENDIF

      smin = s_knots(1) + eps                                               ! Make sure that request is within SPLINE bounds
      smax = (1.d0 - eps)*s_knots(is + ks)
      thmin = u_knots(1)                                                    ! Internally, THETA runs between 0. and 2*pi
      thmax = (1.d0 - eps)*u_knots(itht + ku)
      phimin = v_knots(1)                                                   ! Internally, PHI runs between 0. and 2*pi
      phimax = (1.d0 - eps) *v_knots(iphi + kv)
      IF (s1 < smin .OR. s1 > smax) s1 = smin
      IF (s2 < smin .OR. s2 > smax) s2 = smax
      IF (type_in == 'toroidal') THEN                    
        IF (angle < phimin .OR. angle > phimax)  THEN
         IF(myid .eq. 0) THEN
          WRITE(*, *) 'angle requested:', angle
          WRITE(*, *) 'PHIMIN, PHIMAX::', phimin, phimax
        ENDIF      !IF(myid .eq. 0) THEN
          STOP 'Not enough points in toroidal spline mesh!'
        ENDIF
      ELSEIF (type_in == 'poloidal') THEN
        IF (angle < thmin .OR. angle > thmax) THEN
         IF(myid .eq. 0) THEN
          WRITE(*, *) 'angle requested:', angle
          WRITE(*, *) 'THMIN, THMAX::', thmin, thmax
         ENDIF      !IF(myid .eq. 0) THEN
          STOP 'Not enough points in poloidal spline mesh!'
        ENDIF
      ENDIF

      IF(myid .eq. 0) THEN
      WRITE(*,*)
      WRITE(*,*) '--------------------------------------------------'
     &  , '--------------------------------------------------'
      WRITE(*,*)
      WRITE(*,*) ' RUN INFO:'
      WRITE(*,*) ' ========='
      WRITE(*,*)
      WRITE(*,*) ' Type of integration: ', type_in
      WRITE(*,'("  Number of toroidal periods = ",i6,2x,
     &   "using ",i6, " steps.")') steps, scndry_steps
      WRITE(*,'("  Number of surfaces = ",i5)') nsurfs
      WRITE(*,'("  Innermost radius of surfaces = ", f8.4)') s1
      WRITE(*,'("  Outermost radius of surfaces = ", f8.4)') s2
      WRITE(*,'("  Number of orbits per surface = ",i5)') num_orbs
      ENDIF      !IF(myid .eq. 0) THEN
      local_num_orbs = nsurfs*num_orbs/numprocs
      norbs_tot = nsurfs*num_orbs
      num_eqns = 2*local_num_orbs; norbs = num_eqns/2                      ! NUM_EQNS = Number of equations that will be advanced in time
      IF(numprocs .gt. norbs_tot) THEN
      WRITE(*,*)
       WRITE(*,'("*** More processors (",i5,") than field lines ("
     &      ,i5,")***")') numprocs,norbs_tot
      WRITE(*,*)
       STOP
      ENDIF

      IF(numprocs*int(norbs_tot/numprocs) .ne. norbs_tot) THEN
      WRITE(*,*)
       WRITE(*,'("*** Number of processors (",i5,
     &      ") does not divide evenly into number of field lines ("
     &      ,i5,")***")') numprocs,norbs_tot
      WRITE(*,*)
       STOP
      ENDIF

      ALLOCATE(start_points(norbs_tot,2), STAT=istat)
        IF (istat .NE. 0) STOP 'Allocation error'
      ALLOCATE(theta_nb_tot(norbs_tot), STAT=istat)
        IF (istat .NE. 0) STOP 'Allocation error'
      ALLOCATE(phi_nb_tot(norbs_tot), STAT=istat)
        IF (istat .NE. 0) STOP 'Allocation error'
      ALLOCATE(s_nb_tot(norbs_tot), STAT=istat)
        IF (istat .NE. 0) STOP 'Allocation error'
      ALLOCATE(r_s_nb_tot(norbs_tot), STAT=istat)
        IF (istat .NE. 0) STOP 'Allocation error'
      ALLOCATE(z_s_nb_tot(norbs_tot), STAT=istat)
        IF (istat .NE. 0) STOP 'Allocation error'
      IF(myid .eq. 0) THEN
      WRITE(*,'("  Total number of orbits = ", i6)') norbs
      WRITE(*,'("  Number equations solved by LSODE = ", i6)') num_eqns
      WRITE(*,*)
      WRITE(*,*) '--------------------------------------------------'
     &  , '--------------------------------------------------'
      WRITE(*,*)
      ENDIF      !IF(myid .eq. 0) THEN
!
      i_truncate = 0                                                        ! Remove .txt or .dat of input file to create output file
      ilen = LEN(TRIM(arg1))
      IF (arg1(ilen-3:ilen) == '.dat' .OR.
     &  arg1(ilen-3:ilen) == '.txt') i_truncate = 1
      IF (i_truncate == 0) THEN
        outfile = arg1(1:ilen)//'_puncture_'//TRIM(type_in)//'.dat'
      ELSE
        outfile = arg1(1:ilen-4)//'_puncture_'//TRIM(type_in)//'.dat'
      ENDIF

      INQUIRE(file=TRIM(outfile), exist= l_ex)
      IF (l_ex) THEN
        OPEN(unit = 21, file = TRIM(outfile), status = "old"
     &     , access = "append")                                             ! Open output file as APPEND, since then several runs can be done to complete Poincare plot
      ELSE
        OPEN(unit = 21, file = TRIM(outfile), status = "new")               ! Write out headers if new file
      IF(myid .eq. 0) THEN
        IF (type_in == 'poloidal') THEN
          WRITE(21,*) " &   phi    s   r   phi" 
        ELSEIF (type_in == 'toroidal') THEN
          WRITE(21,*) " &   tht   s   r   z" 
        ENDIF
      ENDIF      !IF(myid .eq. 0) THEN
      ENDIF

      CALL alloc_lsode                                                      ! Prepare LSODE

      icount = 0
      DO j= 1, nsurfs                                                       ! Distribute orbit starting points for all processors
        DO i= 1, num_orbs
          icount = icount + 1
          IF (nsurfs > 1) THEN
            start_points(icount,1) = s1 + (s2-s1)*REAL(j-1, KIND=rprec)/
     &        REAL(nsurfs-1, KIND=rprec)     
          ELSE
           start_points(icount,1) = s1
          ENDIF
           start_points(icount,2) = twopi*REAL(i-1, KIND=rprec)/                        ! If type='toroidal', even y's are THETA values; ='poloidal', they are PHI values.
     &      REAL(num_orbs, KIND=rprec)
        END DO
      END DO
!
      icount = 0
!      write(*,*) myid, myid*local_num_orbs+1,(myid+1)*local_num_orbs
      DO i = myid*local_num_orbs+1,(myid+1)*local_num_orbs                               ! Distribute initial orbit locations
        icount = icount + 1
        y(2*icount-1) = start_points(i,1)
	y(2*icount) = start_points(i,2)
      ENDDO
!      
      IF(myid .eq. 0) THEN
      WRITE(*,*) ' Starting integration....'
      WRITE(*,*)
      ENDIF      !IF(myid .eq. 0) THEN
      ang_out = 0.0d0                                                       ! Start toroidal stepping loop
      delta = twopi/scndry_steps                                            ! Integration takes place with steps in toroidal angle o
#if defined(MPI_OPT)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif
      DO it = 1, steps*scndry_steps                                         ! STEPS =No toroidal transits; SCNDRY_STEPS = steps per
        IF(MOD(it,scndry_steps) == 0) THEN
         IF(myid .eq. 0) THEN
          IF (type_in == 'poloidal')  
     &      WRITE(*,*) '  -> Poloidal transit ', 
     &      it/scndry_steps, ' completed'
          IF (type_in == 'toroidal')  
     &      WRITE(*,*) '  -> Toroidal transit ', 
     &      it/scndry_steps, ' completed'
         ENDIF      !IF(myid .eq. 0) THEN
        ENDIF
        ang_in = ang_out                                                    ! ANG_IN = starting angle value for this loop
        ang_out = ang_in + delta                                            ! ANG_OUT = final angle value for this loop
        rwork(5) = delta; iwork(6) = 5000                                   ! RWORK(5) =  1st step size to try; IWORK(6) =  max. int
        istate = 1; itol = 1; itask = 1; iopt = 0; jt = 10                  ! ITOL = 1, same tolerance for all components; JT = 10,

        IF (type_in == 'poloidal') THEN
          CALL lsode(rt_hand_side_dtheta, num_eqns, y, ang_in, ang_out,
     &      itol, relerr, abserr, itask, istate, iopt, rwork, lrw,
     &      iwork, liw, jacobian, jt)
        ELSEIF (type_in == 'toroidal') THEN
          CALL lsode(rt_hand_side_dphi, num_eqns, y, ang_in, ang_out,
     &      itol, relerr, abserr, itask, istate, iopt, rwork, lrw,
     &      iwork, liw, jacobian, jt)
        ELSE
          STOP 'Bad type!!'
        ENDIF

        IF (istate < 0) THEN
          IF (istate == -1) THEN
            WRITE(*,*) ' Excess work. Wrong JT?'
            WRITE(*,*) ' JT = 10. Non-stiff'
            WRITE(*,*) ' JT = 21-25. Stiff. May need JAC.'
          ELSEIF (istate == -2) THEN
            WRITE(*,*) ' Excess accuracy requested'
          ELSEIF (istate == -3) THEN
            WRITE(*,*) ' Illegal input'
          ELSEIF (istate == -4) THEN
            WRITE(*,*) ' Repeated error test failures'
          ELSEIF (istate == -5) THEN
            WRITE(*,*) ' Repeated convergence failures'
            WRITE(*,*) ' Maybe bad Jacobian?'
          ELSEIF (istate == -6) THEN
            WRITE(*,*) ' Error weight became zero!'
          ENDIF
          WRITE(*,*) 'ANGLE done so far:', ang_in
        ENDIF
!
        IF (type_in == 'poloidal') THEN
          DO i = 1, norbs                                                       ! Keep field lines inside spline-fit domain
            s_test = y(2*i-1); phi_test = y(2*i)
            CALL reset_coords(s_test, rdum, phi_test)
            y(2*i-1) = s_test; y(2*i)= phi_test
          ENDDO
        ELSEIF (type_in == 'toroidal') THEN
          DO i = 1, norbs                                                       ! Keep field lines inside spline-fit domain
            s_test = y(2*i-1); theta_test = y(2*i)
            CALL reset_coords(s_test, theta_test, rdum)
            y(2*i-1) = s_test; y(2*i)= theta_test
          ENDDO
        ENDIF

        IF (MOD(it, scndry_steps) == 0) THEN                                    ! Write out tab-delimited puncture data

          IF (type_in == 'poloidal') THEN
            DO i= 1, norbs
              s_nb(i) =  y(2*i-1)
              theta_nb(i) =  theta_puncture
              phi_nb(i) = y(2*i)
            ENDDO

            CALL dbs3vd(0, 0, 0, norbs, s_nb, theta_nb, phi_nb,
     &        ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi,
     &        r_coef, r_s_nb, is_uniformx, is_uniformy, is_uniformz )

#if defined(MPI_OPT)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
       call MPI_GATHER(phi_nb,norbs,MPI_REAL8,
     >  phi_nb_tot,norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
       call MPI_GATHER(s_nb,norbs,MPI_REAL8,
     >  s_nb_tot,norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
       call MPI_GATHER(r_s_nb,norbs,MPI_REAL8,
     >  r_s_nb_tot,norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
#else
      phi_nb_tot(:) = phi_nb(:)
      s_nb_tot(:) = s_nb(:)
      r_s_nb_tot(:) = r_s_nb(:)
#endif
#if defined(MPI_OPT)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif
           IF(myid .eq. 0) THEN
            DO i= 1, norbs_tot
!              WRITE(21,'(e15.7, 3(3x, e15.7))') y(2*i),                         ! Pairs (phi, s) and (R, phi) at THETA=theta_puncture are written out
!     &          y(2*i-1),  r_s_nb(i), y(2*i)
              WRITE(21,'(e15.7, 3(3x, e15.7))') phi_nb_tot(i),                         ! Pairs (phi, s) and (R, phi) at THETA=theta_puncture are written out
     &          s_nb_tot(i),  r_s_nb_tot(i), phi_nb_tot(i)
            END DO
          ENDIF      !IF(myid .eq. 0) THEN

          ELSEIF (type_in == 'toroidal') THEN

            DO i= 1, norbs
              s_nb(i) =  y(2*i-1)
              theta_nb(i) =  y(2*i)
              phi_nb(i) = phi_puncture                                         
            ENDDO
            CALL dbs3vd(0, 0, 0, norbs, s_nb, theta_nb, phi_nb,            
     &        ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi,
     &        r_coef, r_s_nb, is_uniformx, is_uniformy, is_uniformz )

            CALL dbs3vd(0, 0, 0, norbs, s_nb, theta_nb, phi_nb,
     &        ks, ku, kv, s_knots, u_knots, v_knots, is, itht, iphi,
     &        z_coef, z_s_nb, is_uniformx, is_uniformy, is_uniformz )

#if defined(MPI_OPT)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
       call MPI_GATHER(theta_nb,norbs,MPI_REAL8,
     >  theta_nb_tot,norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
       call MPI_GATHER(s_nb,norbs,MPI_REAL8,
     >  s_nb_tot,norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
       call MPI_GATHER(r_s_nb,norbs,MPI_REAL8,
     >  r_s_nb_tot,norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
       call MPI_GATHER(z_s_nb,norbs,MPI_REAL8,
     >  z_s_nb_tot,norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
#else
      theta_nb_tot(:) = theta_nb(:)
      s_nb_tot(:) = s_nb(:)
      r_s_nb_tot(:) = r_s_nb(:)
      z_s_nb_tot(:) = z_s_nb(:)
#endif
#if defined(MPI_OPT)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif
           IF(myid .eq. 0) THEN
            DO i= 1, norbs_tot        
              WRITE(21,'(e15.7, 3(3x, e15.7))') theta_nb_tot(i),                        ! Pairs (theta,s) and (R, Z) at PHI = phi_puncture are written out
     &          s_nb_tot(i), r_s_nb_tot(i), z_s_nb_tot(i)
            END DO
           ENDIF      !IF(myid .eq. 0) THEN

          ENDIF

        ENDIF

      END DO

      CALL dealloc_splines                                                     ! Liberate memory
      CALL dealloc_lsode
#if defined(MPI_OPT)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif
      DEALLOCATE(s_nb_tot,phi_nb_tot,theta_nb_tot,r_s_nb_tot,
     &   z_s_nb_tot,start_points)
      CLOSE(21)
!
      CALL cpu_time(tend)
      IF(myid .eq. 0) THEN
      WRITE(*,*)
      WRITE(*,*) ' .... completed.'
      WRITE(*,*)
      WRITE(*,'(a15, f12.5, a5)') '  Total time = ', tend - tstart,
     &   ' sec.'

      WRITE(*,*) '=================================================='
     &    , '=================================================='
      ENDIF      !IF(myid .eq. 0) THEN

#if defined(MPI_OPT)
      call MPI_FINALIZE(ierr_mpi)
#endif

      END PROGRAM poincare
