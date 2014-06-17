!*******************************************************************************
!>  @file profiler.f
!>  @brief Contains module @ref profiler.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines functions for measuring an tabulating performance of function and
!>  subroutine calls.
!*******************************************************************************

!  Change to 0 to turn off profiling.
#define PROFILE_ON 1

      MODULE profiler
      USE stel_kinds

      IMPLICIT NONE

!*******************************************************************************
!  profiler module parameters
!*******************************************************************************
!>  Max string length
      INTEGER, PARAMETER :: profiler_string_size = 68
!>  Max number of buckets
      INTEGER, PARAMETER :: profiler_bucket_size = 1000

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) bucket
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Full table of profiled functions.
!-------------------------------------------------------------------------------
      TYPE profiler_bucket
!>  Number of calls to a function.
         INTEGER                              :: number_of_calls = 0
!>  Total time spent in that function.
         REAL (rprec)                         :: total_time = 0.0
!>  Symbol name of that function.
         CHARACTER (len=profiler_string_size) :: symbol_name = ''
      END TYPE

!*******************************************************************************
!  profiler module variables
!*******************************************************************************
!>  OpenMP lock variable.
!$    INTEGER :: lock
!>  Array of buckets to hold the values.
      TYPE (profiler_bucket), DIMENSION(profiler_bucket_size) ::               &
     &      buckets

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref profiler_table object.
!>
!>  Allocates memory for a profiler table.
!>
!>  @returns A pointer to a constructed @ref profiler_table object.
!-------------------------------------------------------------------------------
      SUBROUTINE profiler_construct()

      IMPLICIT NONE

!  Start of executable code
!$    CALL OMP_INIT_LOCK(lock)

      END SUBROUTINE

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref profiler_table object.
!>
!>  Deallocates memory for a @ref profiler_table object.
!>
!>  @param[inout] this A @ref profiler_table instance.
!-------------------------------------------------------------------------------
      SUBROUTINE profiler_destruct()

      IMPLICIT NONE

!  Start of executable code
!$    CALL OMP_DESTROY_LOCK(lock)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Gets the end time of profiled function.
!>
!>  Sets the end time of profiled function and adds the result to the table. The
!>  symbol name is hashed to determine the index of the array it lives in. If
!>  the index is full, the next instance is used. When profiling is turned off,
!>  this does nothing. To avoid race conditons in multi threaded code, this must
!>  only be run by one thread or process at a time.
!>
!>  @param[in]    symbol_name Name of the symbol to profile.
!>  @param[in]    start_time  Starting time of the symbol call.
!-------------------------------------------------------------------------------
      SUBROUTINE profiler_set_stop_time(symbol_name, start_time)

!  Declare Arguments
      CHARACTER (len=*), INTENT(in)        :: symbol_name
      REAL (rprec), INTENT(in)             :: start_time

!  local variables
      INTEGER                              :: i, index
      REAL (rprec)                         :: time

!  Start of executable code.
#if PROFILE_ON
       CALL second0(time)
       time = time - start_time
       index = profiler_hash_function(symbol_name)

!$     CALL OMP_SET_LOCK(lock)
       IF (buckets(index)%symbol_name .eq. '') THEN
          buckets(index)%symbol_name = symbol_name
          buckets(index)%number_of_calls = 1
          buckets(index)%total_time = time
       ELSE
          DO i = index, profiler_bucket_size
             IF (buckets(i)%symbol_name .eq. symbol_name) THEN
                buckets(i)%number_of_calls = buckets(i)%number_of_calls        &
     &                                     + 1
                buckets(i)%total_time = buckets(i)%total_time + time
!$              CALL OMP_UNSET_LOCK(lock)
                RETURN
             END IF
          END DO

          DO i = 1, index - 1
             IF (buckets(i)%symbol_name .eq. symbol_name) THEN
                buckets(i)%number_of_calls = buckets(i)%number_of_calls        &
     &                                     + 1
                buckets(i)%total_time = buckets(i)%total_time + time
!$              CALL OMP_UNSET_LOCK(lock)
                RETURN
             END IF
          END DO

          WRITE (*,*) 'Profile table full: ' // TRIM(symbol_name) //           &
     &                ' could not be profiled.'
       END IF
!$     CALL OMP_UNSET_LOCK(lock)
#endif

       END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Gets the start time of profiled function.
!>
!>  Gets a starting time of a function. When profiling is turned off, this
!>  returns 0 instead.
!>
!>  @returns The start time of a function.
!-------------------------------------------------------------------------------
      FUNCTION profiler_get_start_time()

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: profiler_get_start_time

!  Start of executable code
#if PROFILE_ON
      CALL second0(profiler_get_start_time)
#else
      profiler_get_start_time = 0.0
#endif

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Computes a hash for the symbol name.
!>
!>  The hash is computed using the djb2 algorithm.
!>
!>      hash_0 = 5381
!>      hash_i+1 = hash_i * 33 ^ c_i
!>
!>  The magic numbers of 5381 and 33 help the hash function avoid collision and
!>  good distribution.
!>
!>  @param[in] symbol_name Name of the symbol to hash.
!>  @returns The hash value of the symbol.
!-------------------------------------------------------------------------------
      FUNCTION profiler_hash_function(symbol_name)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: profiler_hash_function
      CHARACTER (len=*), INTENT(in) :: symbol_name

!  local variables
      INTEGER(KIND=8)               :: hash
      INTEGER                       :: i

!  Start of executable code
      hash = 5381
      DO i = 1, LEN_TRIM(symbol_name)
         hash = (LSHIFT(hash, 5) + hash) + ICHAR(symbol_name(i:i))
      END DO

      profiler_hash_function = MOD(ABS(hash), profiler_bucket_size) + 1

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write out the profiled data to an output file.
!>
!>  This loops through the buckets and prints out the profiled data. Any empty
!>  buckets get skiped.
!>
!>  @param[in] iou  Input/output unit of the output file.
!-------------------------------------------------------------------------------
      SUBROUTINE profiler_write(iou)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)               :: iou

!  local variables
      INTEGER                           :: i

!  Start of executable code
#if PROFILE_ON
      WRITE (iou,*)
      WRITE (iou,*) ' *** Code Profile Information'

      WRITE (iou,1000)

      DO i = 1, profiler_bucket_size
         IF (buckets(i)%symbol_name .ne. '') THEN
            WRITE (iou, 1001) buckets(i)%symbol_name,                          &
     &                        buckets(i)%total_time,                           &
     &                        buckets(i)%number_of_calls
         END IF
      END DO

1000  FORMAT('Symbol Name',59x,'Time',10x,'Num Calls')
1001  FORMAT(a68,2x,es12.5,2x,i8)
#endif

      END SUBROUTINE

      END MODULE