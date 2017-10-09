      MODULE hessian
      USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i,          &
          nfp=>nfp_i
      USE timer_mod
      USE stel_constants, ONLY: pi
      USE descriptor_mod
      USE shared_data, ONLY: col_scale, neqs, lcolscale, l_ApplyPrecon, &
                             l_linearize, fsq_total1, l_getfsq,         &
                             ndims, mblk_size, lverbose
#if defined(MPI_OPT)
      USE prof_mod
      USE ptrd_mod
#endif
#if defined(SKS)
      USE nscalingtools, ONLY: PARSOLVER, MPI_ERR, PARFUNCTISL,         &
        rcounts, disp, startglobrow, endglobrow, nranks, rank,          &
        SKSDBG, TOFU, UPPER, LOWER, DIAG, SAVEDIAG, nrecd,              &
        SYMMETRYCHECK, totRecvs, PACKSIZE, WriteTime, send, receive,    &
        SetBlockTriDataStruct, GetFullSolution
      USE blocktridiagonalsolver_s, ONLY: GetColSum, ParallelScaling,   &
        ForwardSolve, SetMatrixRHS, BackwardSolve, GetSolutionVector,   &
        CheckSymmetry, Dmin_TRI, MAXEIGEN_TRI, FindMinMax_Tri
      USE mpi_inc
#else
      USE nscalingtools, ONLY: PARSOLVER, PARFUNCTISL, SKSDBG, TOFU,    &
                               startglobrow, endglobrow
#endif
      IMPLICIT NONE
!-----------------------------------------------
!SPH 090116: EVENTUALLY MOVE TO nscalingtools? ASK SKS
      PUBLIC  :: gather_array
      INTERFACE gather_array
         MODULE PROCEDURE gather_array1, gather_array2, gather_array1A, &
                          gather_array2A, gather_array3
      END INTERFACE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      PRIVATE
      INTEGER, PARAMETER :: jstart(3) = (/1,2,3/)
      INTEGER, ALLOCATABLE :: ipiv_blk(:,:)

      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: ublk, dblk, lblk
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: ublk_s, dblk_s, lblk_s
      LOGICAL :: l_backslv = .FALSE. 

      REAL(dp), ALLOCATABLE, DIMENSION(:) :: DataItem

      INTEGER :: myStartBlock, myEndBlock, mynumBlockRows 
      INTEGER :: mystart(3), myend(3)
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: SavedDiag

      REAL(dp) :: eps, ton, toff, lmdamp
      REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: colsum
      REAL(dp), PARAMETER, PUBLIC :: eps_factor = 1._dp 
      REAL(dp), ALLOCATABLE, PUBLIC :: mupar_norm(:)
      REAL(dp), PUBLIC  :: levmarq_param, levmarq_param0, asym_index,   &
                           mupar, mupar0, maxeigen, mindiag
      LOGICAL, PUBLIC   :: l_Compute_Hessian, l_Diagonal_Only
      INTEGER, PUBLIC :: HESSPASS=0
      PUBLIC  :: apply_precond, dealloc_hessian, InitHess,              &
                 Compute_Hessian_Blocks, PadSides

      CONTAINS

      SUBROUTINE gather_array1A(gc)
      REAL(dp), DIMENSION(mblk_size*ns), INTENT(INOUT) :: gc
!-----------------------------------------------
      CALL gather_array1(gc)

      END SUBROUTINE gather_array1A


      SUBROUTINE gather_array2A(section, gc)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size*ns), INTENT(OUT) :: gc
      REAL(dp), DIMENSION(mblk_size*mynumBlockRows), INTENT(IN) :: section
!-----------------------------------------------
      CALL gather_array2(section, gc)

      END SUBROUTINE gather_array2A


      SUBROUTINE gather_array1(gc)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(INOUT) :: gc
      INTEGER :: istat
      REAL(dp), DIMENSION (:,:), ALLOCATABLE :: section
!-----------------------------------------------
      IF (PARSOLVER) THEN
#if defined(SKS)
        ALLOCATE (section(mblk_size,myStartBlock:myEndBlock), stat=istat)
        IF (istat .NE. 0) CALL ErrorAbort('Allocation error in gather_array')

        section(:,myStartBlock:myEndBlock) = gc(:,myStartBlock:myEndBlock)
        CALL gather_array2(section, gc)
        
        DEALLOCATE (section)
#endif
      END IF
      
      END SUBROUTINE gather_array1


      SUBROUTINE gather_array2(section, gc)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(OUT) :: gc
      REAL(dp), DIMENSION(mblk_size,myStartBlock:myEndBlock),           &
                INTENT(IN) :: section
!-----------------------------------------------
      IF (PARSOLVER) THEN
#if defined(SKS)
         CALL MPI_ALLGATHERV(section, mynumBlockRows*mblk_size,         &
              MPI_REAL8, gc,rcounts, disp, MPI_REAL8, SIESTA_COMM,      &
              MPI_ERR)
#endif
      ELSE
         gc = section
      ENDIF
     
      END SUBROUTINE gather_array2

      SUBROUTINE gather_array3(section, gather3d)
      USE nscalingtools, ONLY: mnspcounts, mnspdisps 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:mpol, -ntor:ntor, ns), INTENT(OUT) :: gather3d
      REAL(dp), DIMENSION(0:mpol, -ntor:ntor, myStartBlock:myEndBlock), &
      INTENT(IN)  :: section
      INTEGER :: M2
!-----------------------------------------------
       IF (PARSOLVER) THEN
#if defined(SKS)
          M2 = SIZE(gather3d,1)*SIZE(gather3d,2)
          CALL MPI_ALLGATHERV(section, mynumBlockRows*M2, MPI_REAL8,    &
                           gather3d, mnspcounts, mnspdisps, MPI_REAL8,  &
                           SIESTA_COMM, MPI_ERR)
#endif
       ELSE
          gather3d = section
       END IF
      END SUBROUTINE gather_array3

      SUBROUTINE InitHess
!-----------------------------------------------
      IF (.NOT.PARSOLVER) THEN
         startglobrow = 1; endglobrow = ns
      END IF
      myStartBlock=startglobrow
      myEndBlock=endglobrow
      mynumBlockRows=myEndBlock-myStartBlock+1

      END SUBROUTINE InitHess

      SUBROUTINE Compute_Hessian_Blocks (func, l_Hess_sym, ldiagonal)
      USE shared_data, ONLY: xc, gc, l_linearize, lasym
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(IN)     :: l_Hess_sym, ldiagonal
      EXTERNAL                :: func
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      INTEGER                 :: iunit
      CHARACTER(LEN=3)        :: label
      REAL(dp)                :: bsize
!-----------------------------------------------

      l_Compute_Hessian = .TRUE.
      l_Diagonal_Only   = ldiagonal

      l_linearize = .TRUE.
      IF (lColScale) col_scale = 1

      IF (iam .EQ. 0) THEN
         DO iunit = 6,33,27
            IF (.not.lverbose .and. iunit .eq. 6) CYCLE
            IF (l_Diagonal_Only) THEN
               WRITE (iunit, 90) levmarq_param, muPar
            ELSE
               WRITE (iunit, 100) levmarq_param, muPar, asym_index
            ENDIF
         END DO
         CALL FLUSH(33)
      ENDIF

 90   FORMAT (/,' Computing diagonal preconditioner - ',                &
                ' LM parameter:',1pe9.2,' mu||:',1pe9.2)
 100  FORMAT (/,' Computing block preconditioner - ',                   &
                ' LM parameter:',1pe9.2,' mu||:',1pe9.2,                &
                ' Asym index:',1pe9.2)

      IF (PARSOLVER) THEN
        IF (PARFUNCTISL) THEN
          IF(SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_With_No_Col_Redist'; IF(SKSDBG) CALL FLUSH(TOFU)
          IF (iam.EQ.0 .AND. HESSPASS.EQ.0 .and. lverbose)              &
             PRINT *,'HESSIAN AND INVERSE CALCULATION USING BCYCLIC WITH NO COLUMN REDISTRIBUTION'
          CALL Compute_Hessian_Blocks_With_No_Col_Redist (xc, gc, func, l_Hess_sym)
        ELSE
          IF (SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_With_Col_Redist'; IF(SKSDBG) CALL FLUSH(TOFU)
          IF (iam.EQ.0 .AND. HESSPASS.EQ.0 .and. lverbose)              &
             PRINT *,'HESSIAN AND INVERSE CALCULATION USING BCYCLIC WITH COLUMN REDISTRIBUTION'
          CALL Compute_Hessian_Blocks_With_Col_Redist (xc, gc, func, l_Hess_sym)
        END IF 
      ELSE 
        IF (SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_Thomas'; CALL FLUSH(TOFU)
        IF (iam.EQ.0 .AND. HESSPASS.EQ.0 .and. lverbose)                &
           PRINT *,'HESSIAN AND INVERSE CALCULATION USING SERIAL THOMAS ALGORITHM'
        CALL Compute_Hessian_Blocks_Thomas (xc, gc, func, l_Hess_sym)
      END IF

      IF (HESSPASS.EQ.0) THEN
         bsize = 3*ns*KIND(dblk)                                         !3 blocks per row
         bsize = bsize*REAL(mblk_size*mblk_size,dp)
         IF (bsize .LT. 1.E6_dp) THEN
            bsize = bsize/1.E3
            label = " Kb"
         ELSE IF (bsize .lt. 1.E9_dp) THEN
            bsize = bsize/1.E6_dp
            label = " Mb"
         ELSE
            bsize = bsize/1.E9_dp
            label = " Gb"
         END IF

         IF (iam .EQ. 0) THEN
            DO iunit = 6,33,27
               IF (.not.lverbose .and. iunit .eq. 6) CYCLE
               WRITE (iunit, '(1x,a,i4,a,f6.2,a)') 'Block dim: ', mblk_size, &
                  '^2  Preconditioner size: ', bsize, TRIM(label)
            END DO
            CALL FLUSH(33)
         ENDIF 
      END IF

      l_Compute_Hessian = .FALSE.
      l_linearize = .FALSE.
      HESSPASS = HESSPASS+1

      END SUBROUTINE Compute_Hessian_Blocks

      SUBROUTINE Compute_Hessian_Blocks_With_No_Col_Redist (xc, gc,     &
                 func, l_Hess_sym)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION((1+mpol)*(2*ntor+1)*ndims,ns) :: xc, gc
      LOGICAL, INTENT(IN) :: l_Hess_sym
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(dp), PARAMETER    :: zero=0, one=1, p5=0.5_dp
      INTEGER  :: n, m, js, mesh, ntype, istat, iunit, icol
      REAL(dp) :: eps, ton, toff, bsize, lmdamp
      CHARACTER(LEN=3) :: label
      EXTERNAL    :: func
      INTEGER  :: js1
#if defined(SKS)
      REAL(dp) :: starttime, endtime, usedtime
      REAL(dp) :: colstarttime, colendtime, colusedtime
      INTEGER :: nsmin, nsmax, dcnt
      REAL(dp) :: skston, skstoff, temp, DSaved(mblk_size)
      INTEGER :: icolmpi, i1, i2, i3, i4
!----------------------------------------------

      starttime=0; endtime=0; usedtime=0; 
      colstarttime=0; colendtime=0; colusedtime=0

      !------------------------------------------------------------
      !     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
      !------------------------------------------------------------
      CALL second0(ton)

      IF (mblk_size .NE. SIZE(gc,1)) CALL ErrorAbort('MBLK_SIZE IS WRONG IN HESSIAN')
      eps = SQRT(EPSILON(eps))*ABS(eps_factor)

      ALLOCATE (DataItem(mblk_size), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('DataItem allocation failed!')
      DataItem = 0
      ALLOCATE (SavedDiag(mblk_size), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('SavedDiag allocation failed!')
      SavedDiag = 0

      nsmin = MAX(1, startglobrow-1);  nsmax = MIN(ns, endglobrow+1)

      IF (HESSPASS .EQ. 0) THEN
        myend = nsmax
        DO mesh = 1, 3
           icol = MOD(jstart(mesh)-nsmin, 3)
           IF (icol .LT. 0) icol = icol+3
           mystart(mesh) = nsmin+icol
        END DO
      ENDIF

      !----------------------------
      ! Column generation begins
      !----------------------------
      CALL second0(colstarttime)

      icol=0
      xc = 0

      VAR_TYPE_LG: DO ntype = 1, ndims
        VAR_N_LG: DO n = -ntor, ntor
          VAR_M_LG: DO m = 0, mpol

            icol=icol+1

            MESH_3PT_LG: DO mesh = 1, 3

              DO js = mystart(mesh), myend(mesh), 3
                 xc(icol,js) = eps
              END DO

              INHESSIAN=.TRUE.
              CALL second0(skston)
              IF (m.EQ.0 .AND. n.LT.0) THEN
                 gc(:,nsmin:nsmax) = 0
              ELSE
                 CALL func
              END IF
              CALL second0(skstoff)
              hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
              INHESSIAN=.FALSE.

              DO js = mystart(mesh), myend(mesh), 3
                 xc(icol,js) = 0
              END DO

              SKIP3_MESH_LG: DO js = mystart(mesh), myend(mesh), 3

                !ublk(js-1)
                js1=js-1
                IF (startglobrow.LE.js1 .AND. js1.LE.endglobrow) THEN
                   DataItem = gc(:,js1)/eps
                   IF (l_diagonal_only) THEN
                      temp=DataItem(icol)
                      DataItem=0
                      DataItem(icol)=temp
                   END IF
                   CALL SetBlockTriDataStruct(UPPER,js1,icol,DataItem) 
                END IF

                !dblk(js)
                IF (startglobrow.LE.js .AND. js.LE.endglobrow) THEN
                   DataItem = gc(:,js)/eps

                   IF (.NOT. lColScale) THEN                               
                      IF (ALL(DataItem .EQ. zero)) DataItem(icol) = -one
                      lmdamp = levmarq_param*DataItem(icol)
                      IF (levmarq_param .NE. zero) THEN
                         IF (lmdamp .EQ. zero) THEN
                           lmdamp = levmarq_param
                         ELSE IF (ntype.NE.3 .AND. ntype.NE.6) THEN
                           lmdamp = lmdamp/10                               !Emulate v|| damping only
                         END IF
                      END IF
				   END IF

                   IF (l_diagonal_only) THEN
                      temp=DataItem(icol)
                      DataItem=0
                      DataItem(icol)=temp
                   END IF
                   CALL SetBlockTriDataStruct(SAVEDIAG,js,icol,DataItem)
                   IF (.NOT. lColScale) DataItem(icol) = DataItem(icol) &
	 				                 + SIGN(lmdamp,DataItem(icol))
                   CALL SetBlockTriDataStruct(DIAG,js,icol,DataItem)
                END IF

                !lblk(js+1) 
                js1 = js+1
                IF (startglobrow.LE.js1 .AND. js1.LE.endglobrow) THEN
                   DataItem = gc(:,js1)/eps
                   IF (l_diagonal_only) THEN
                      temp=DataItem(icol)
                      DataItem=0
                      DataItem(icol)=temp
                   END IF
                   CALL SetBlockTriDataStruct(LOWER,js1,icol,DataItem)
                END IF

              END DO SKIP3_MESH_LG
            END DO MESH_3PT_LG
          END DO VAR_M_LG
        END DO VAR_N_LG
      END DO VAR_TYPE_LG

      DEALLOCATE (DataItem, stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('DataItem deallocation error!')
      DEALLOCATE (SavedDiag, stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('SavedDiag deallocation error!')
      !----------------------------
      ! Column generation ends
      !----------------------------

      !IF(DIAGONALDONE) CALL WriteBlocks (HESSPASS)

      CALL second0(colendtime)
      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)

      IF (lColScale) THEN
         ALLOCATE (colsum(mblk_size,startglobrow:endglobrow))
         CALL GetColSum(colsum)
         CALL VectorCopyPar (colsum, col_scale)
         CALL ParallelScaling(levmarq_param,colsum)
         DEALLOCATE(colsum)
      ELSE 
         CALL FindMinMax_TRI (levmarq_param)
      END IF

      mindiag = Dmin_TRI
      maxeigen = MAXEIGEN_TRI

      !     CHECK SYMMETRY OF BLOCKS 

      IF (SYMMETRYCHECK) THEN
        CALL second0(skston)
        CALL CheckSymmetry(asym_index)
        CALL second0(skstoff)
        asymmetry_check_time=asymmetry_check_time+(skstoff-skston)
      ENDIF      

      !
      !     SYMMETRIZE BLOCKS
      !

      CALL second0(toff)
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

!      IF (l_Diagonal_Only) RETURN         !ForwardSolve will be called first time but safer this way!

      ton = toff
      skston = ton
!
!     FACTORIZE (Invert) HESSIAN
!
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE(ipiv_blk,stat=istat)
      CALL ForwardSolve

      CALL second0(toff)
      skstoff=toff

      IF (iam .EQ. 0 .and. lverbose) PRINT '(a,1p,e12.3)',' BLOCK FACTORIZATION TIME: ', toff-ton

      block_factorization_time=block_factorization_time+(skstoff-skston)
      time_factor_blocks = time_factor_blocks + (toff-ton)

#endif
      END SUBROUTINE Compute_Hessian_Blocks_With_No_Col_Redist

      SUBROUTINE PadSides(arrin, blksize, top, bot)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)     :: top, bot, blksize
      REAL(dp), INTENT(INOUT) :: arrin(blksize, ns)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
#if defined(SKS)
      INTEGER                 :: MPI_STAT(MPI_STATUS_SIZE)
      INTEGER                 :: M1, left, right, tag1=1, tag2=2,       &
                                 tr, t1l, tl, t1r
!-----------------------------------------------
      left=rank-1;  IF(rank.EQ.0) left=MPI_PROC_NULL
      right=rank+1; IF(rank.EQ.nranks-1) right=MPI_PROC_NULL

      tr = endglobrow;              tl = startglobrow
      t1r= MIN(endglobrow+1, ns);   t1l= MAX(1, startglobrow-1)

!Send top / Rec bottom elements
      M1 = blksize*top
      CALL MPI_Sendrecv(arrin(:,tr-(top-1):tr), M1, MPI_REAL8,          &
               right, tag1, arrin(:,t1l-(top-1):t1l), M1, MPI_REAL8,    &
               left,  tag1, SIESTA_COMM, MPI_STAT, MPI_ERR)

!Send bottom / Rec top boundary elements
      M1 = blksize*bot
      CALL MPI_Sendrecv(arrin(:,tl:tl+(bot-1)), M1, MPI_REAL8,          &
               left,  tag2, arrin(:,t1r:t1r+(bot-1)), M1, MPI_REAL8,    &
               right, tag2, SIESTA_COMM, MPI_STAT, MPI_ERR)
#endif
      END SUBROUTINE PadSides

#if defined(SKS)
      SUBROUTINE VectorCopyPar (colsum, colscale)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: colsum(mblk_size,startglobrow:endglobrow)
      REAL(dp), INTENT(OUT) :: colscale(mblk_size,ns)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: js
!-----------------------------------------------
      DO js = startglobrow, endglobrow 
        colscale(:,js) = colsum(:,js)
      END DO

      CALL PadSides(colscale, mblk_size, 1, 1)

      END SUBROUTINE VectorCopyPar
#endif

      SUBROUTINE Compute_Hessian_Blocks_With_Col_Redist (xc, gc, func,  &
                 l_Hess_sym)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION((1+mpol)*(2*ntor+1)*ndims,ns) :: xc, gc
      LOGICAL, INTENT(IN) :: l_Hess_sym
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(dp), PARAMETER    :: zero=0, one=1, p5=0.5_dp
      INTEGER  :: n, m, js, js1, mesh, ntype, istat, iunit, icol, icolmpi
      EXTERNAL    :: func
#if defined(SKS)
      INTEGER, ALLOCATABLE, DIMENSION (:) :: sendCount, recvCount
      REAL(dp), ALLOCATABLE, DIMENSION (:) :: recvBuf
      INTEGER :: procID
      REAL(dp) :: starttime, endtime, usedtime
      REAL(dp) :: colstarttime, colendtime, colusedtime

      LOGICAL :: PROBEFLAG
      REAL(dp) :: skston, skstoff, temp

      INTEGER :: it 
!----------------------------------------------
      starttime=zero
      endtime=zero
      usedtime=zero
!------------------------------------------------------------
!     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
!------------------------------------------------------------

      CALL second0(ton)

      IF (mblk_size .NE. SIZE(gc,1)) CALL ErrorAbort('MBLK_SIZE IS WRONG IN HESSIAN')

      eps = SQRT(EPSILON(eps))*ABS(eps_factor)

      ALLOCATE (DataItem(mblk_size), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('DataItem allocation failed!')
      DataItem = 0
      ALLOCATE (SavedDiag(mblk_size), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('SavedDiag allocation failed!')
      SavedDiag = zero

      icolmpi = 0
      icol=0
      xc = 0

      IF(.NOT.ALLOCATED(sendCount)) ALLOCATE(sendCount(nprocs))
      IF(.NOT.ALLOCATED(recvCount)) ALLOCATE(recvCount(nprocs))
      IF(.NOT.ALLOCATED(recvBuf)) ALLOCATE(recvBuf(PACKSIZE))
      sendCount=0
      recvCount=0
      nrecd=0
      totRecvs=0
      PROBEFLAG=.TRUE.

      CALL second0(colstarttime)
      VAR_TYPE: DO ntype = 1, ndims
         VAR_N: DO n = -ntor, ntor
            VAR_M: DO m = 0, mpol
                  
               icol=icol+1
               IF(MOD(icol-1,nprocs)==iam) THEN
                  icolmpi = icolmpi+1

                  MESH_3PT: DO mesh = 1, 3

                  DO js = jstart(mesh), ns, 3
                     xc(icol,js) = eps
                  END DO

                  INHESSIAN=.TRUE.
                  CALL second0(skston)
                  CALL func
                  CALL second0(skstoff)
                  hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
                  INHESSIAN=.FALSE.

                  DO js = mystart(mesh), myend(mesh), 3
                     xc(icol,js) = 0
                  END DO

                  SKIP3_MESH: DO js = jstart(mesh), ns, 3

                     !ublk(js-1)
                     js1 = js-1
                     IF (js1 .gt. 0) THEN
                        DataItem = gc(:,js1)/eps
                        !m=0 constraint for n<0
                        IF (m.eq.0 .and. n.lt.0) DataItem = 0
                        IF (l_diagonal_only) THEN
                           temp=DataItem(icol)
                           DataItem=0
                           DataItem(icol)=temp
                        END IF
                        istat = 1
                        CALL receive(PROBEFLAG)
                        CALL send(DataItem,js1,istat,icol,procID)
                        IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                        CALL receive(PROBEFLAG)
                     END IF !js>1

                     !dblk(js)
                     DataItem = gc(:,js)/eps
                     IF (m.eq.0 .and. n.lt.0) DataItem = 0

                     IF (.NOT. lColScale) THEN                               
                        IF (ALL(DataItem .EQ. zero)) DataItem(icol) = -one
                        lmdamp = levmarq_param*DataItem(icol)
                        IF (levmarq_param .NE. zero) THEN
                           IF (lmdamp .EQ. zero) THEN
                              lmdamp = levmarq_param
                           ELSE IF (ntype.NE.3 .AND. ntype.NE.6) THEN
                              lmdamp = lmdamp/10                               !Emulate v|| damping only
                           END IF
                        END IF
				     END IF
                     
                     IF (l_diagonal_only) THEN
                        temp=DataItem(icol)
                        DataItem=0
                        DataItem(icol)=temp
                     END IF
                     SavedDiag=DataItem; istat=4
                     CALL receive(PROBEFLAG)
                     CALL send(SavedDiag,js,istat,icol,procID)
                     IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                     CALL receive(PROBEFLAG)
!                    Boundary condition at js=1 and ns (CATCH ALL OF THEM HERE)
!                    and m=0,n<0: NEED THIS to avoid (near) singular Hessian
!                    ASSUMES DIAGONALS ARE ALL NEGATIVE
                     IF (.NOT. lColScale) DataItem(icol) = DataItem(icol) &
                                        + SIGN(lmdamp,DataItem(icol))
                     istat = 2; js1 = js
                     CALL receive(PROBEFLAG)
                     CALL send(DataItem,js1,istat,icol,procID)
                     IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                     CALL receive(PROBEFLAG)
                   
                     !lblk(js+1)
                     js1 =js+1 
                     IF (js .lt. ns) THEN
                       DataItem = gc(:,js1)/eps
                       !m=0 constraint for n<0
                       IF (m.eq.0 .and. n.lt.0) DataItem=0
                       IF (l_diagonal_only) THEN
                          temp=DataItem(icol)
                          DataItem=0
                          DataItem(icol)=temp
                       END IF
                       istat = 3
                       CALL receive(PROBEFLAG)
                       CALL send(DataItem,js1,istat,icol,procID)
                       IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                       CALL receive(PROBEFLAG)
                     END IF !JS < NS

                  END DO SKIP3_MESH

               END DO MESH_3PT
               ENDIF
            END DO VAR_M
         END DO VAR_N
      END DO VAR_TYPE
      CALL second0(colendtime)
      colusedtime=colendtime-colstarttime

      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)

      CALL second0(skston)
      IF (PARSOLVER) THEN

        PROBEFLAG=.NOT.PROBEFLAG

        CALL MPI_AllTOALL(sendCount,1,MPI_INTEGER,recvCount,1,          &
		                  MPI_INTEGER,SIESTA_COMM,MPI_ERR)
        totRecvs=0
        DO js=1,nprocs,1
          totRecvs=totRecvs+recvCount(js)
        END DO

        DO WHILE (nrecd.LT.totRecvs)
          CALL receive(PROBEFLAG)
        END DO
      END IF

      DEALLOCATE (DataItem, stat=istat)
      DEALLOCATE (SavedDiag, stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('DataItem deallocation error!')

      CALL second0(skstoff)
      construct_hessian_time=construct_hessian_time+(skstoff-skston)
      CALL MPI_Barrier(SIESTA_COMM, MPI_ERR)

!
!     CHECK SYMMETRY OF BLOCKS 
!
      IF (SYMMETRYCHECK) THEN
        IF (PARSOLVER) THEN 
          CALL second0(skston)
          CALL CheckSymmetry(asym_index)
          CALL second0(skstoff)
          asymmetry_check_time=asymmetry_check_time+(skstoff-skston)
        ENDIF
      ENDIF      

!
!     FACTORIZE (Invert) HESSIAN
!
      CALL second0(toff)
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

      ton=toff
      starttime=ton
      skston=ton
      CALL ForwardSolve
      CALL second0(toff)
      endtime=toff
      skstoff=toff

      IF (iam.EQ.0 .AND. lverbose) PRINT '(a,1p,e12.3)',' BLOCK FACTORIZATION TIME: ', toff-ton

      block_factorization_time=block_factorization_time+(skstoff-skston)
      time_factor_blocks = time_factor_blocks + (toff-ton)
      usedtime=endtime-starttime
#endif
      END SUBROUTINE Compute_Hessian_Blocks_With_Col_Redist

      SUBROUTINE Compute_Hessian_Blocks_Thomas (xc, gc, func,           &
                 l_Hess_sym)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION((1+mpol)*(2*ntor+1)*ndims,ns),                &
               INTENT(INOUT)  :: xc, gc
      LOGICAL, INTENT(IN) :: l_Hess_sym
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(dp), PARAMETER    :: zero=0, one=1, p5=0.5_dp
      INTEGER  :: n, m, js, js1, mesh, ntype, istat, iunit, icol, icolmpi
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)  :: gc1
      REAL(dp) :: skston, skstoff, temp
      EXTERNAL :: func
#if defined(MPI_OPT)
      INTEGER  :: neqtotal, numroc
      REAL(dp) :: starttime, endtime, usedtime
      REAL(dp) :: colstarttime, colendtime, colusedtime
      EXTERNAL :: numroc
!------------------------------------------------------------

      starttime=0; endtime=0; usedtime=0
      colstarttime=0; colendtime=0; colusedtime=0
#endif

!------------------------------------------------------------
!     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
!------------------------------------------------------------
      CALL second0(ton)

      IF (mblk_size .NE. SIZE(gc,1)) CALL ErrorAbort('MBLK_SIZE IS WRONG IN HESSIAN')
!      mblk_size = SIZE(gc,1)

      eps = SQRT(EPSILON(eps))*ABS(eps_factor)

      IF (ALLOCATED(ublk)) DEALLOCATE(ublk, dblk, lblk, stat=istat)
      IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
        CALL blacs_gridinfo(icontxt_1xp,nprow,npcol,myrow,mycol)
        mb = mblk_size
        nb = 1

        rsrc = 0
        csrc = 0
        Locq = numroc( mblk_size, nb, mycol, csrc, npcol )
        Locp = numroc( mblk_size, mb, myrow, rsrc, nprow )
        mblk_size2=max(1,Locq)
        lld = max(1,Locp)
        call descinit(descA_1xp,mblk_size,mblk_size,mb,nb,rsrc,csrc,    &
                      icontxt_1xp,lld,info)
        if (info.ne.0 .and. lverbose) then
          write(*,*) 'myrow,mycol,nprow,npcol,desc(LLD_),lld ',         &
             myrow,mycol,nprow,npcol,descA_1xp(LLD_),lld
          write(*,*) 'Locp,m,mb ', Locp,mblk_size,mb
        endif

        call assert(info.eq.0,'descinit descA_1xp ',info)
        ineed = max(1,lld*Locq)
        mm = mblk_size*ns
        mb=mm
        nb=nrhs1
        LocqR = numroc( nrhs1, nb, mycol,csrc,npcol)
        LocpR = numroc( mm, mb, myrow,rsrc,nprow)
        lld = max(1,LocpR)
        call descinit(descR_all,mm,nrhs1, mb,nb,rsrc,csrc,icontxt,lld,info)
        call assert(info.eq.0,'test_pdtrd:descinit return info=',info)

        call blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)

        mb = 50
        nb = mb

        rsrc = 0
        csrc = 0
        Locq = numroc( mblk_size, nb, mycol, csrc, npcol )
        Locp = numroc( mblk_size, mb, myrow, rsrc, nprow )

        lld = MAX(1,Locp)
        call descinit(descA,mblk_size,mblk_size,mb,nb,rsrc,csrc,icontxt,lld,info)
        call assert(info.eq.0,'descinit descA ',info)
        ineed = MAX(1,lld*Locq)
        IF (ALLOCATED(ublkp)) DEALLOCATE(ublkp, dblkp, lblkp, stat=istat)
        IF (.NOT. ALLOCATED(ublkp)) THEN
           ALLOCATE(ublkp(ineed,ns-1), dblkp(ineed,ns),                 &
                    lblkp(ineed,ns-1), stat=istat)
           IF (istat .NE. 0) THEN
             CALL ErrorAbort('Not enough memory to store a single block!')
           END IF
        END IF

        ALLOCATE(ublk(mblk_size,mblk_size2,ns),                         &
                 dblk(mblk_size,mblk_size2,ns),                         &
                 lblk(mblk_size,mblk_size2,ns),                         &
                 stat=istat)
        IF (istat .NE. 0) THEN
           CALL ErrorAbort('Not enough memory to store a single block!')
        END IF
#else
        CALL ErrorAbort('LSCALAPACK=T BUT NOT MPI!')
#endif
      ELSE   !.NOT.LSCALAPACK

        ALLOCATE(ublk(mblk_size,mblk_size,ns),                          &
                 dblk(mblk_size,mblk_size,ns),                          &
                 lblk(mblk_size,mblk_size,ns),                          &
                 stat=istat)
        IF (istat .NE. 0) CALL ErrorAbort('Unable to allocate blocks!')
      
      END IF       !END LSCALAPACK

      IF (ALLOCATED(ublk)) THEN
         ublk = 0; dblk = 0; lblk = 0
      END IF

      ALLOCATE (DataItem(mblk_size), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('DataItem allocation failed!')
      DataItem = 0

#if defined(MPI_OPT)
      icolmpi = 0
#endif
      icol=0
      xc = 0

#if defined(SKS)
      CALL second0(colstarttime) !Added by SKS for timing comparisons
#endif
      VAR_TYPE: DO ntype = 1, ndims
         VAR_N: DO n = -ntor, ntor
            VAR_M: DO m = 0, mpol
                  
               icol=icol+1
#if defined(MPI_OPT)
               IF(MOD(icol-1,nprocs)==iam .OR. .NOT.LSCALAPACK) THEN
                  IF (LSCALAPACK) THEN
                     icolmpi = icolmpi+1
                  ELSE 
                     icolmpi = icol
                  END IF
#else
                  icolmpi = icol
#endif
                  MESH_3PT: DO mesh = 1, 3

                  DO js = jstart(mesh), ns, 3
                     xc(icol,js) = eps
                  END DO

                  INHESSIAN=.TRUE.
                  CALL second0(skston)
                  CALL func
                  CALL second0(skstoff)
#if defined(SKS)
                  hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
#endif
                  INHESSIAN=.FALSE.
                  DO js = jstart(mesh), ns, 3
                     xc(icol,js) = 0
                  END DO


!              COMPUTE PRECONDITIONER (HESSIAN) ELEMENTS. LINEARIZED EQUATIONS
!              OF TRI-DIAGONAL (IN S) FORM, a=lblk, d=dblk, b=ublk
!
!              F(j-1) = a(j-1)x(j-2) + d(j-1)x(j-1) + b(j-1)x(j)
!              F(j)   =                a(j)x(j-1)   + d(j)  x(j)  + b(j)  x(j+1)
!              F(j+1) =                               a(j+1)x(j)  + d(j+1)x(j+1) + b(j+1)x(j+2)
!
!              HESSIAN IS H(k,j) == dF(k)/dx(j); aj == lblk; dj == dblk; bj = ublk

                  SKIP3_MESH: DO js = jstart(mesh), ns, 3

!FORCE RESPONSE AT mp,np,nptype TO m,n,js,ntype VELOCITY PERTURBATION

                     !ublk(js-1)
                     js1 = js-1
                     IF (js1 .GT. 0) THEN
                        DataItem = gc(:,js1)/eps

                     !m=0 constraint for n<0
                        IF (m.eq.0 .and. n.lt.0) DataItem = 0
                        IF (l_Diagonal_Only) THEN
                           temp=DataItem(icol)
                           DataItem=0
                           DataItem(icol)=temp
                        END IF
                        ublk(:,icolmpi,js1) = DataItem
                     END IF

                     !dblk(js)
                     DataItem = gc(:,js)/eps
                     IF (m.eq.0 .and. n.lt.0) DataItem = 0

!                    Boundary condition at js=1 and ns (CATCH ALL OF THEM HERE)
!                    and m=0,n<0: NEED THIS to avoid (near) singular Hessian
!                    ASSUMES DIAGONALS ARE ALL NEGATIVE

                     IF (.NOT. lColScale) THEN                               
                        IF (ALL(DataItem .EQ. zero)) DataItem(icol) = -one
                        lmdamp = levmarq_param*DataItem(icol)
                        IF (levmarq_param .NE. zero) THEN
                           IF (lmdamp .EQ. zero) THEN
                             lmdamp = levmarq_param
                           ELSE IF (ntype.NE.3 .AND. ntype.NE.6) THEN
                             lmdamp = lmdamp/10                               !Emulate v|| damping only
                           END IF
                        END IF
                        DataItem(icol) = DataItem(icol)                 &
						               + SIGN(lmdamp,DataItem(icol))
					 END IF

                     IF (l_Diagonal_Only) THEN
                        temp=DataItem(icol)
                        DataItem=0
                        DataItem(icol)=temp
                     END IF

                     dblk(:,icolmpi,js) = DataItem

                     !lblk(js+1) 
                     js1=js+1
                     IF (js .LT. ns) THEN
                        DataItem = gc(:,js1)/eps
                     !m=0 constraint for n<0
                        IF (m.eq.0 .and. n.lt.0) DataItem = 0
                        IF (l_Diagonal_Only) THEN
                           temp=DataItem(icol)
                           DataItem=0
                           DataItem(icol)=temp
                        ENDIF
                        lblk(:,icolmpi,js1) = DataItem
                     END IF

                  END DO SKIP3_MESH

               END DO MESH_3PT
#if defined(MPI_OPT)
               ENDIF
#endif
            END DO VAR_M
         END DO VAR_N
      END DO VAR_TYPE
      
      CALL second0(toff)
      
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

!SPH 071416 ADD COLUMN SCALING
      IF (lColScale) THEN
         IF (nprocs .GT. 1) CALL ErrorAbort('Thomas COL-Scaling not implemented for nprocs>1')
         CALL SerialScaling(mindiag, maxeigen)
      END IF
!
!SPH 053116 COMPUTE 2-NORMS OF EACH COLUMN FOR SCALING
!
      ton = toff
#if defined(SKS)
      CALL second0(colendtime)
      colusedtime=colendtime-colstarttime
      time_generate_blocks=time_generate_blocks+colusedtime
      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)
#endif

      DEALLOCATE (DataItem, stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('DataItem deallocation error!')

!
!     CHECK SYMMETRY OF BLOCKS 
!
      CALL second0(skston)
      IF (LSCALAPACK) THEN
         CALL check3d_symmetry(dblk,lblk,ublk,dblkp,lblkp,ublkp,        &
                               mblk_size,mblk_size2,ns,asym_index)
      ELSE
         CALL check3d_symmetry_serial(dblk,lblk,ublk,mpol,ntor,         & 
                               mblk_size,ns,asym_index)
      END IF
      CALL second0(skstoff)
#if defined(SKS)
      asymmetry_check_time=asymmetry_check_time+(skstoff-skston)
#endif

!
!     SYMMETRIZE BLOCKS
!
      IF (.NOT. l_Hess_sym) GOTO 1000
      icol = 0
      ALLOCATE(gc1(mblk_size,ns))
      
      DO m = 0, mpol
         DO n = -ntor, ntor
            DO ntype = 1, ndims
               icol = icol+1
               gc1(:,1:nsh) = (ublk(icol,:,1:nsh)                       &
                            +  lblk(:,icol,2:ns))/2
               ublk(icol,:,1:nsh) = gc1(:,1:nsh)
               lblk(:,icol,2:ns)  = gc1(:,1:nsh)
               gc1(:,:) = (dblk(icol,:,:)                               &
                       +  dblk(:,icol,:))/2
               dblk(icol,:,:) = gc1(:,:)
               dblk(:,icol,:) = gc1(:,:)
            END DO
         END DO
      END DO

      DEALLOCATE (gc1)
 1000 CONTINUE

!      IF (l_Diagonal_Only) RETURN


!      IF (nprecon .eq. 4) l_backslv = .TRUE.
      IF (l_backslv) THEN
         ALLOCATE(ublk_s(mblk_size,mblk_size,ns),                       &
                  dblk_s(mblk_size,mblk_size,ns),                       &
                  lblk_s(mblk_size,mblk_size,ns),                       &   
                  stat=istat)
         IF (istat .ne. 0) CALL ErrorAbort('Allocation error2 in Compute_Hessian_Blocks')
         ublk_s = ublk;  dblk_s = dblk; lblk_s = lblk
      END IF
!
!     FACTORIZE (Invert) HESSIAN
!
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE(ipiv_blk,stat=istat)
      CALL second0(skston)
      IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
         ineed = numroc(descA(M_),descA(MB_),0,0,nprow) + descA(MB_)
         ALLOCATE( ipiv_blk(ineed, ns), stat=istat )
         IF (istat .NE. 0) CALL ErrorAbort('ipiv_blk allocation error2 in Compute_Hessian_Blocks')
         CALL blk3d_parallel(dblk,lblk,ublk,dblkp,lblkp,ublkp,mblk_size,mblk_size2,ns)
         CALL blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)
         mm = mblk_size*ns
         LocqR = numroc( nrhs1, nb, mycol,csrc,npcol)
         LocpR = numroc( mm, mb, myrow,rsrc,nprow)
         lld = MAX(1,LocpR)
         CALL descinit(descX,mm,nrhs1, mb,nb,rsrc,csrc,icontxt,lld,info)
         CALL assert(info.eq.0,'test_pdtrd:descinit return info=',info)
         descR(:) = descX

         ineedR = MAX(1,lld*LocqR)
 
         CALL blacs_barrier(icontxt,'All')
         CALL profstart('factor call:123')
         CALL pdtrdf(mblk_size,ns,dblkp,lblkp,ublkp,ipiv_blk,descA)
         CALL blacs_barrier(icontxt,'All')
         CALL profend('factor call:123')
#else
         CALL ErrorAbort('LSCALAPACK=T BUT NOT IN MPI!')
#endif
      ELSE    !.NOT. LSCALAPACK
         ALLOCATE (ipiv_blk(mblk_size,ns), stat=istat)
         IF (istat .ne. 0) CALL ErrorAbort('ipiv_blk allocation error2 in Compute_Hessian_Blocks')
         CALL blk3d_factor(dblk, lblk, ublk, ipiv_blk, mblk_size, ns)
      END IF  ! LSCALAPACK

      CALL second0(skstoff)
	  toff = skstoff
      IF (iam.EQ.0 .AND. lverbose) PRINT '(a,1p,e12.3)',' BLOCK FACTORIZATION TIME: ', toff-ton

      block_factorization_time=block_factorization_time+(skstoff-skston)
      time_factor_blocks = time_factor_blocks + (toff-ton)

      END SUBROUTINE Compute_Hessian_Blocks_Thomas

      
      SUBROUTINE SerialScaling (mindiag, maxeigen)

      REAL(dp), INTENT(OUT)  :: mindiag, maxeigen
      REAL(dp), PARAMETER    :: minThreshold=1.E-6_dp
      REAL(dp), PARAMETER    :: zero=0, one=1, p5=0.5_dp
      INTEGER       :: js, icol
      REAL(dp)      :: eps, eps0, Dmin_TRI, DminL

      maxeigen = 0
      Dmin_TRI = HUGE(Dmin_TRI)

      ALLOCATE(colsum(mblk_size,ns))

      BLOCK_ROW: DO js = 1, ns    !sum over block rows
      COLS: DO icol = 1, mblk_size
         colsum(icol,js) = SUM(dblk(:,icol,js)**2)
         IF (js .GT. 1) THEN
            colsum(icol,js) = colsum(icol,js) + SUM(ublk(:,icol,js-1)**2)         
         END IF
         IF (js .LT. ns) THEN
            colsum(icol,js) = colsum(icol,js) + SUM(lblk(:,icol,js+1)**2)         
         END IF
      ENDDO COLS
         colsum(:,js) = SQRT(colsum(:,js))
      ENDDO BLOCK_ROW

      EPS0 = MAXVAL(colsum)
      DO js = 1, ns
         eps = 1.E-6_dp*MAXVAL(colsum(:,js))
         IF (eps .EQ. ZERO) eps = eps0
         WHERE (colsum(:,js) .LT. eps) colsum(:,js) = eps
      END DO

      COLSUM = one/SQRT(COLSUM)

      CALL VECTOR_COPY_SER (colsum, col_scale)

 !    SCALE BLOCKS
      BLOCK_ROW2: DO js = 1, ns
      COLS2: DO icol = 1, mblk_size
         dblk(:,icol,js) = colsum(:,js)*dblk(:,icol,js)*colsum(icol,js)
         IF (js .GT. 1) THEN
             lblk(:,icol,js) = colsum(:,js)*lblk(:,icol,js)*colsum(icol,js-1)
         END IF
         IF (js .LT. ns) THEN
             ublk(:,icol,js) = colsum(:,js)*ublk(:,icol,js)*colsum(icol,js+1)
         END IF
      END DO COLS2
      END DO BLOCK_ROW2

      BLOCK_ROW3: DO js = 1, ns
         mindiag = 0
         DO icol = 1, mblk_size
            mindiag = MAX(mindiag, ABS(dblk(icol,icol,js)))
         END DO

         mindiag = minThreshold*mindiag
         IF (minDiag .EQ. 0) minDiag = minThreshold

      COLS3: DO icol = 1, mblk_size
         eps = dblk(icol,icol,js) 
	 IF (ABS(eps) .LE. mindiag) THEN
	    dblk(icol,icol,js) = -mindiag
         ELSE IF (js .LT. ns) THEN
	    Dmin_TRI = MIN(MAX(1.E-12_dp,ABS(eps)), Dmin_TRI)
	 END IF
      END DO COLS3
      END DO BLOCK_ROW3

      mindiag=SIGN(Dmin_TRI, eps)
      DminL = levmarq_param*mindiag

!SPH061016: RESTORE LM PARAMETER
      DO js = 1, ns
	  DO icol = 1, mblk_size
          eps = dblk(icol,icol,js)
          dblk(icol,icol,js) = eps + SIGN(DminL,eps)
	  END DO
      END DO

      DEALLOCATE(colsum)

      END SUBROUTINE SerialScaling

      SUBROUTINE VECTOR_COPY_SER(colsum, colscale)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: COLSUM(mblk_size,ns)
      REAL(dp), INTENT(OUT) :: COLSCALE(mblk_size,ns)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: js
!-----------------------------------------------
 
      DO js = 1, ns
        COLSCALE(:,js) = COLSUM(:,js)
      END DO

      END SUBROUTINE

!#define SVD_ON
!#undef SVD_ON

      SUBROUTINE blk3d_factor(a, bm1, bp1, ipiv, mblk, nblocks)
      USE stel_constants, ONLY: one, zero
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER  :: bytes_per_rprec2 = 8
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(OUT) :: ipiv(mblk,nblocks)
      REAL(dp), TARGET, DIMENSION(mblk,mblk,nblocks),                   &
                     INTENT(INOUT) :: a, bm1, bp1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, k1, ier
      INTEGER, POINTER :: ipivot(:)
      REAL(dp), POINTER :: amat(:,:), bmat(:,:), cmat(:,:)
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: temp
#if defined(SVD_ON)
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: w, work
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: v, u
      REAL(dp), PARAMETER :: small = 1.E-300_dp
      INTEGER :: nw
      INTEGER :: lwork
      CHARACTER(LEN=1), PARAMETER :: jobu='A', jobvt='A'
#endif
!-----------------------------------------------
!  modified (June, 2003, ORNL):         S. P. Hirshman
!-----------------------------------------------------------------------
!
!  this subroutine solves for the Q factors of a block-tridiagonal system of equations.
!
!-----------------------------------------------------------------------
!  INPUT
!  mblk                : block size
!  nblocks             : number of blocks
!  a                   : diagonal blocks
!  bp1, bm1            : lower, upper blocks (see equation below)
!
!  OUTPUT
!  ipiv                : pivot elements for kth block
!  a                   : a-1 LU factor blocks
!  bm1                 : q = a-1 * bm1 matrix
!
!  LOCAL VARIABLES
!  iunit               : unit number for block-tridiagonal solution disk file.
!
!  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
!  equation is:
!
!           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
!
!     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
!
!     1. Start from row N and solve for x(N) in terms of x(N-1):
!
!        x(N) = -q(N)*x(N-1) + r(N)
!
!        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
!
!        where a(N)[-1] is the inverse of a(N)
!
!     2. Substitute for lth row to get recursion equation for q(l) and r(l):
!
!        x(l) = -q(l)*x(l-1) + r(l), in general, where:
!
!        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
!
!        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
!
!        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
!
!     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
!
!        x(1) = r(1)
!
!     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
!
!
!     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
!
!     IF USING SVD:
!
!     1. CALL dgesvd:   Perform SVD decomposition
!     2. CALL svbksb:   Solve Ax = bi, for i=1,mblk
!
!     OTHERWISE
!
!     1. CALL dgetrf:   Perform LU factorization of diagonal block (A) - faster than sgefa
!     2. CALL dgetrs:   With multiple (mblk) right-hand sides, to do block inversion
!                         operation, Ax = B  (stores result in B; here B is a matrix)
!
!      ndisk = mblk

!  create disk file for doing direct access i/o.

!      incnow = ndisk
!      irecl  = bytes_per_rprec2*incnow
!      incbu = 1 + (ndisk - 1)/incnow
!      ibuph = 0

!      iunit = 10
!      CALL safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted',
!     1     irecl, 'DIRECT')
!      IF (ier .ne. 0) CALL ErrorAbort('Error opening scratch file in blk3d')

!  main loop. load and process (backwards) block-rows nblocks to 1. 

#if defined(SVD_ON)
      lwork = 10*mblk
      ALLOCATE(v(mblk,mblk), w(mblk),   &
     &         u(mblk,mblk), work(lwork), stat=ier)
      IF (ier .NE. 0) CALL ErrorAbort('Allocation error in blk3d_factor!')
#endif
      ipiv = 0
      ALLOCATE (temp(mblk,mblk), stat=ier)
      IF (ier .NE. 0) CALL ErrorAbort('Allocation error in blk3d_factor!')

      BLOCKS: DO k = nblocks, 1, -1
!
!     Compute (and save) qblk(nblocks) = ablk(nblocks)[-1] * bml
!
         amat => a(:,:,k)

#if defined(SVD_ON)
!NOTE: v = vt coming out here
         CALL dgesvd (jobu, jobvt, mblk, mblk, amat, mblk, w, u, mblk,  &
        &             v, mblk, work, lwork, ier)
!Set SVD weights w to 0 for weights below the allowed threshold, so backsolver
!will compute pseudo-inverse
         WRITE (35, '(a,i4,a,1pe12.4)') 'Block: ',k, ' Condition #: ', SQRT(w(1)/w(mblk))
         DO nw = 2, mblk
            IF (w(nw) .le. small*w(1)) w(nw) = 0
         END DO
!
!        STORE svd pseudo-inverse IN AMAT since u,v,w will be deallocated at end
         CALL svdinv2 (amat, u, v, w, mblk)
#else
         ipivot => ipiv(:,k)
         CALL dgetrf (mblk, mblk, amat, mblk, ipivot, ier)
#endif
         IF (ier .ne. 0) GOTO 200
         IF (k .eq. 1) EXIT

         bmat => bm1(:,:,k)

#if defined(SVD_ON)
!
!        Use (pseudo) inverse stored in AMAT = V*1/w*Ut 
!
         temp = bmat
         bmat = MATMUL(amat, temp)
#else
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot,               &
     &                bmat, mblk, ier)
         IF (ier .ne. 0) GOTO 305
#endif
!         CALL wrdisk(iunit, ql, ndisk, incnow, ibuph, incbu, ier)
!         IF (ier .ne. 0) GOTO 302

!
!      Update effective diagonal "a" matrix. Use dgemm: faster AND doesn't overflow normal stack
!
         k1 = k-1 
         amat => bp1(:,:,k1)
         cmat => a(:,:,k1)
!         cmat = cmat - MATMUL(amat, bmat)
         CALL dgemm('N','N',mblk,mblk,mblk,-one,amat,mblk,              &
     &              bmat, mblk, one, cmat, mblk)

      END DO BLOCKS

#if defined(SVD_ON)
      DEALLOCATE(v, w, u, work)
#endif
!
!     COMPUTE TRANSPOSES HERE, SINCE REPEATEDLY CALLING MATMUL OPERATION
!     X*At IS FASTER THAN A*X DUE TO UNIT STRIDE
!

      DO k = 1, nblocks
         IF (k .ne. nblocks) THEN
            temp = TRANSPOSE(bp1(:,:,k))
            bp1(:,:,k) = temp
         END IF
         IF (k .ne. 1) THEN
            temp = TRANSPOSE(bm1(:,:,k))
            bm1(:,:,k) = temp
         END IF
      END DO

      DEALLOCATE (temp)

      GOTO 400

!  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6, '(2x,a,i4)') 'Error factoring matrix in blk3d: block = '&
                             , k
      IF (ier < 0)WRITE (6,'(i4, a)') ier, 'th argument has illegal value'
      IF (ier > 0)WRITE (6,'(i4, a)') ier, 'th diagonal factor exactly zero'
      CALL ErrorAbort('Error in blk3d_factor')
  301 CONTINUE
!      WRITE (6, '(a,i8)') ' BLK3D:   error in opening file:  ',
!     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine WRDISK'
  303 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine RDDISK'
      ier = -2
  305 CONTINUE
      WRITE (6, '(2/a,i10,2/)')                                         &
      ' BLK3D-FACTOR:   error detected:   ier =', ier
      CALL ErrorAbort('Error in blk3d_factor')

!  destroy disk file and return. ---------------------------------------
 
  400 CONTINUE

!      CLOSE (iunit)

      END SUBROUTINE blk3d_factor


      SUBROUTINE block_precond(gc)
      USE timer_mod
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(INOUT) :: gc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER  ::  zero = 0
      INTEGER :: m, n, js, ntype, istat, icol
      REAL(dp) :: t1, error, error_sum, b_sum
      REAL(dp), ALLOCATABLE :: gc_s(:,:)
#if defined(MPI_OPT)
      REAL(dp) :: alpha0,beta0
      INTEGER  :: i,j,k,nn,kk,ia,ja,ix,jx,ic,jc,ib,jb
#endif
#if defined(SKS)
      REAL(dp) :: starttime, endtime, usedtime
      INTEGER :: PRECONDPASS=0
      INTEGER :: globrow, cnt
      !Might need to make it allocatable - just a test right now.
      REAL(dp), DIMENSION (:,:), ALLOCATABLE :: solvec
!-----------------------------------------------
      starttime=0
      endtime=0
#endif
!-----------------------------------------------
!
!     Applies 3D block-preconditioner to forces vector (gc)
!
#if defined(SKS)
      IF (PARSOLVER) THEN

        DO globrow=myStartBlock,myEndBlock
          CALL SetMatrixRHS (globrow,gc(1:mblk_size,globrow))
        END DO
        
        CALL BackwardSolve
       
        ALLOCATE (solvec(mblk_size,mynumBlockRows), stat=istat)
        IF (istat .NE. 0) CALL ErrorAbort('Allocation error in block_precond before gather')

        cnt=1
        DO globrow=myStartBlock,myEndBlock
          CALL GetSolutionVector (globrow,solvec(1:mblk_size,cnt))
          cnt = cnt + 1
        END DO

        CALL gather_array(solvec, gc)
!        CALL MPI_Allgatherv(solvec,mynumBlockRows*mblk_size,MPI_REAL8,gc, &
!             rcounts, disp, MPI_REAL8, SIESTA_COMM, MPI_ERR)

!!        CALL GetFullSolution(temp,gc)
        
        DEALLOCATE (solvec)

      ELSE

        CALL second0(starttime)
#endif
         IF (l_backslv) THEN
            ALLOCATE (gc_s(mblk_size,ns), stat=istat)
            IF (istat .NE. 0) CALL ErrorAbort('Allocation error0 in block_precond')
            gc_s = gc
         END IF

!     Solve Hx = gc, using Hessian (H) LU factors stored in block_... matrices
!     Store solution "x" in gc on output
         IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
         ALLOCATE (tempp(ineedR), stat=istat)
         ia = 1
         ja = 1
         ib = 1
         jb = 1
         beta0 = 0.0d0
         alpha0 = 1.0d0
         mm0=descR(3)
         nrhs0=descR(4)
         call profstart('pdgeadd call:123')
         call pdgeadd( 'N', mm0,nrhs0, alpha0,gc,ia,ja,descR_all,       &
                       beta0,tempp,ib,jb,descR )
         call profend('pdgeadd call:123')
         ir = 1
         jr = 1
         call blacs_barrier(icontxt,'All')
         call profstart('solver call:123')
         call pdtrds(mblk_size,ns,dblkp,lblkp,ublkp,ipiv_blk,descA,     &
                     nrhs1,tempp,ir,jr,descR )
         call blacs_barrier(icontxt,'All')
         call profend('solver call:123')
         gc = 0
         ia = 1
         ja = 1
         ib = 1
         jb = 1
         beta0 = 0.0d0
         alpha0 = 1.0d0
         mm0=descR(3)
         nrhs0=descR(4)
!     descR_all(7)=-1
!     descR_all(8)=-1
         call profstart('pdgeadd2 call:123')
         call pdgeadd( 'N', mm0,nrhs0, alpha0, tempp,ia,ja,descR,       &
                       beta0,gc,ib,jb,descR_all )
         call dgsum2d( icontxt,'A', ' ', mm,1,gc,mm,       -1,-1)
         call profend('pdgeadd2 call:123')

         DEALLOCATE (tempp, stat=istat)
#else
         CALL ErrorAbort('MPI_OPT must be true for LSCALAPACK!')
#endif
         ELSE       !NOT LSCALAPACK

#if defined(SKS)
         CALL second0(endtime)
         usedtime=endtime-starttime
         IF (SKSDBG) WRITE(TOFU,*)                                      &
         'HessianTag-2 : Time to BackwardSolve:',usedtime,              &
         'PRECONDPASS',PRECONDPASS,'in native run'
         IF (SKSDBG) CALL FLUSH(TOFU)
#endif
!        Serial solver
         CALL blk3d_slv(dblk, lblk, ublk, gc, ipiv_blk, mblk_size, ns)
         END IF     !END LSCALAPACK

#if defined(SKS)
      END IF !END OF IF(PARSOLVER) CONDITIONAL
      PRECONDPASS = PRECONDPASS + 1
#endif

#if !defined(MPI_OPT)
      IF (l_backslv) THEN
         error_sum = 0;  b_sum = SUM(gc_s*gc_s)
         l_backslv = .false.
         icol = 0

         IF (lverbose) WRITE (6, *) ' Writing block Hessian check to unit 34'
         WRITE (34, *)
         WRITE (34, *) ' BLK3D FACTORIZATION CHECK: Ax = b ?'
        
         DO n = -ntor, ntor
            WRITE (34, *) ' N = ', n
            DO m = 0, mpol
               WRITE (34, *) ' M = ', m
               DO ntype = 1, ndims
                  icol = icol+1
                  WRITE (34, *) ' TYPE = ', ntype
                  WRITE (34, *)                                         &
                  '   js        x             Ax             b     '//  &
                  '     Ax - b        RelErr'
                  js = 1
                  t1 = SUM(dblk_s(icol,:,js)*gc(:,js)                   &
                     + ublk_s(icol,:,js)*gc(:,js+1))

                  error = t1 - gc_s(icol,js)
                  error_sum = error_sum + error*error
                  IF (t1 .eq. zero) t1 = EPSILON(t1)
                  IF (ABS(error) .gt. 1.E-10_dp) THEN
                  WRITE (34, 101) js, gc(icol,js), t1,                  &
                                  gc_s(icol,js), error, error/t1
                  ELSE
                  WRITE (34, 100) js, gc(icol,js), t1,                  &
                                  gc_s(icol,js), error, error/t1
                  END IF

               DO js = 2, ns-1
                  t1 = SUM(lblk_s(icol,:,js)*gc(:,js-1)                 &
                         + dblk_s(icol,:,js)*gc(:,js)                   &
                         + ublk_s(icol,:,js)*gc(:,js+1))
                  error = t1 - gc_s(icol,js)
                  error_sum = error_sum + error*error
                  IF (t1 .eq. zero) t1 = EPSILON(t1)
                  IF (ABS(error) .gt. 1.E-10_dp) THEN
                  WRITE (34, 101) js, gc(icol,js), t1,                  &
                         gc_s(icol,js), error, error/t1
                  ELSE
                  WRITE (34, 100) js, gc(icol,js), t1,                  &
                         gc_s(icol,js), error, error/t1
                  END IF
               END DO

                 js = ns        
                 t1 = SUM(lblk_s(icol,:,js)*gc(:,js-1)                  &
                        + dblk_s(icol,:,js)*gc(:,js))
                 error = t1 - gc_s(icol,js)
                 error_sum = error_sum + error*error
                 IF (t1 .eq. zero) t1 = EPSILON(t1)
                 IF (ABS(error) .gt. 1.E-10_dp) THEN
                 WRITE (34, 101) js, gc(icol,js), t1,                   &
                        gc_s(icol,js), error, error/t1
                 ELSE
                 WRITE (34, 100) js, gc(icol,js), t1,                   &
                        gc_s(icol,js), error, error/t1
                 END IF
              END DO
            END DO
         END DO

         DEALLOCATE(dblk_s, lblk_s, ublk_s, gc_s, stat=istat)

         IF (lverbose) PRINT *,' Hessian Error: ', SQRT(error_sum/b_sum)

      END IF
#endif

 100  FORMAT(i6,1p,5e14.4)
 101  FORMAT(i6,1p,5e14.4,'  *')

      END SUBROUTINE block_precond

 
      SUBROUTINE apply_precond(gc)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(inout) :: gc
!      LOGICAL, INTENT(in)                              :: ltype
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: ton, toff
!-----------------------------------------------
      CALL second0(ton)
!      IF (ltype) THEN
!         CALL diag_precond(gc)
!      ELSE
         CALL block_precond(gc)
!      END IF 
      CALL second0(toff)
      time_apply_precon = time_apply_precon+(toff-ton)

      END SUBROUTINE apply_precond


      SUBROUTINE blk3d_slv(ablk, qblk, bp1, source, ipiv, mblk, nblocks)
      USE stel_kinds
      USE stel_constants, ONLY: one
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: bytes_per_rprec2 = 8
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(dp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(IN) ::     &
                                     ablk, qblk, bp1
      REAL(dp), DIMENSION(mblk,nblocks), TARGET, INTENT(INOUT) :: source
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, POINTER  :: ipivot(:)
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, ier
      REAL(dp), POINTER :: amat(:,:), x1(:), source_ptr(:)
!-----------------------------------------------
!  modified (June, 2003, ORNL):         S. P. Hirshman
!-----------------------------------------------------------------------
!
!  this subroutine solves a block-tridiagonal system of equations, using 
!  the ABLK, QBLK factors from blk3d_factor,
!
!-----------------------------------------------------------------------
!  INPUT
!  mblk                : block size
!  nblocks             : number of blocks
!  bp1                 : upper blocks (see equation below)
!  ipiv                : pivot elements for kth block
!  ablk                : a-1 blocks
!  qblk                : q = a-1 * bm1
!  source              : input right side
!
!  OUTPUT
!  source              : Solution x of A x = source
! 
!  LOCAL VARIABLES
!  iunit               : unit number for block-tridiagonal solution disk file.
!
!  the tri-diagonal equation is:
!
!           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
!
!     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
!
!     1. Start from row N and solve for x(N) in terms of x(N-1):
!
!        x(N) = -q(N)*x(N-1) + r(N)
!
!        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
!
!        where a(N)[-1] is the inverse of a(N)
!
!     2. Substitute for lth row to get recursion equation fo q(l) and r(l):
!
!        x(l) = -q(l)*x(l-1) + r(l), in general, where:
!
!        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
!
!        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
!
!        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
!
!     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
!
!        x(1) = r(1)
!
!     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
!
!
!     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
!
!     1. CALL dgetrs:   With single right hand side (source) to solve A x = b (b a vector)
!                       Faster than dgesl
!      ndisk = mblk

!  create disk file for doing direct access i/o.

!      incnow = ndisk
!      irecl  = bytes_per_rprec2*incnow
!      incbu = 1 + (ndisk - 1)/incnow
!      ibuph = 0

!      iunit = 10
!      CALL safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted',
!     1     irecl, 'DIRECT')
!      IF (ier .ne. 0) CALL ErrorAbort('Error opening scratch file in blk3d')

!  main loop. load and process (backwards) block-rows nblocks to 1. 
!  note: about equal time is spent in calling dgetrs and in performing
!  the two loop sums: on ibm-pc, 2 s (trs) vs 3 s (sums); on linux (logjam),
!  2.4 s (trs) vs 3 s (sums).

!
!     Back-solve for modified sources first
!
      BLOCKS: DO k = nblocks, 1, -1

         source_ptr => source(:,k)
         amat => ablk(:,:,k)
#if defined(SVD_ON)
         !source_sp = MATMUL(amat, source_sp)
#else
         ipivot => ipiv(:,k);   
         CALL dgetrs('n', mblk, 1, amat, mblk, ipivot, source_ptr, mblk, ier)
         IF (ier .ne. 0) GOTO 305
#endif
         IF (k .eq. 1) EXIT
!
!        NOTE: IN BLK3D_FACTOR, BP1 AND BM1 WERE TRANSPOSED (AND STORED)
!        TO MAKE FIRST INDEX FASTEST VARYING IN THE FOLLOWING MATMUL OPS
!
         amat => bp1(:,:,k-1)
         x1 => source(:,k);  source_ptr => source(:,k-1)
!         source_ptr = source_ptr - MATMUL(x1,amat)  !USE THIS FORM IF TRANSPOSED bp1
!         source_ptr = source_ptr - MATMUL(amat,x1)  !UNTRANSPOSED FORM
         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,one,source_ptr,1)

      END DO BLOCKS
!
!  forward (back-substitution) solution sweep for block-rows k = 2 to nblocks
!  now, source contains the solution vector
!
      DO k = 2, nblocks
!         CALL rddisk (iunit, ql, ndisk, incnow, ibuph, ier)
!         IF (ier .ne. 0) GOTO 303
!         ibuph = ibuph - incbu

         amat => qblk(:,:,k)
         x1 => source(:,k-1);  source_ptr => source(:,k)
!         source_ptr = source_ptr - MATMUL(x1,amat)  !USE THIS FORM IF TRANSPOSED qblk
!         source_ptr = source_ptr - MATMUL(amat,x1)  !UNTRANSPOSED FORM
         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,one,source_ptr,1)

      END DO

      GOTO 400

!  error returns. ------------------------------------------------------

  301 CONTINUE
!      WRITE (6, '(a,i8)') ' BLK3D:   error in opening file:  ',
!     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine WRDISK'
  303 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine RDDISK'
      ier = -2
  305 CONTINUE
      WRITE (6, '(2/a,i10,2/)') ' BLK3D-SLV:   error detected:   ier =', ier
      CALL ErrorAbort('Error in blk3d_slv')

!  destroy disk file and return. ---------------------------------------

  400 CONTINUE

!      CLOSE (iunit)

      END SUBROUTINE blk3d_slv

      
      SUBROUTINE dealloc_hessian
      INTEGER :: istat

      IF (ALLOCATED(ublk)) DEALLOCATE(ublk, dblk, lblk, stat=istat)
      IF (ALLOCATED(ublkp))DEALLOCATE(dblkp,lblkp,ublkp,stat=istat)
      IF (ALLOCATED(mupar_norm)) DEALLOCATE(mupar_norm, stat=istat)
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE (ipiv_blk)

      END SUBROUTINE dealloc_hessian

#if defined(MPI_OPT)
      SUBROUTINE blk3d_parallel(a,b,c,ap,bp,cp,mblk,mblkc,nblocks)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: ipart,jpart,i,j,k
      INTEGER :: ia,ja,inc1,ib,jb,inc2,nsm
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: nblocks, mblk, mblkc
        REAL(dp) :: aij2,aij
      REAL(dp), TARGET, DIMENSION(mblk,mblkc,nblocks),                  &
               INTENT(INOUT) ::  a, b, c
      REAL(dp), TARGET, DIMENSION(Locq,Locp,nblocks),                   &
               INTENT(INOUT) ::  ap
      REAL(dp), TARGET, DIMENSION(Locq,Locp,nblocks-1),                 &
               INTENT(INOUT) ::  bp, cp
      REAL(dp), POINTER :: amat(:,:), bmat(:,:), cmat(:,:)
      REAL(dp), POINTER :: amatp(:,:), bmatp(:,:), cmatp(:,:)

      do k=1,nblocks
        call profstart('change context:123')
        call blacs_barrier(descA(CTXT_),'All')
         amat => a(:,:,k)
         amatp => ap(:,:,k)
          call pdgemr2d(mblk,mblk,amat,1,1,descA_1xp,                   &
                       amatp,1,1,descA,  icontxt_global)
        if(k<nblocks)then
         bmat => b(:,:,k+1)
         bmatp => bp(:,:,k)
          call pdgemr2d(mblk,mblk,bmat,1,1,descA_1xp,                   &
                        bmatp,1,1,descA,  icontxt_global)
         cmat => c(:,:,k)
         cmatp => cp(:,:,k)
          call pdgemr2d(mblk,mblk,cmat,1,1,descA_1xp,                   &
                        cmatp,1,1,descA,  icontxt_global)
       endif
        call blacs_barrier(descA(CTXT_),'All')
        call profend('change context:123')
!        do ja=1,mblk
!          call pdelget( 'A',' ',aij, amatp,ja,ja,descA )
!          call pdelget( 'A',' ',aij2, amat,ja,ja,descA_1xp )
!              write(*,*) 'ja,ja,aij ',ja,ja,aij,aij2
!        enddo
       end do
      END SUBROUTINE blk3d_parallel
#endif

      END MODULE hessian
