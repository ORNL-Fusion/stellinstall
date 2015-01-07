!*******************************************************************************
!>  @file pprofile_T.f
!>  @brief Contains module @ref pprofile_T.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref pprofile_class. This module contains
!>  all the code necessary to define parameterized profiles.
!*******************************************************************************
      MODULE pprofile_T
      USE stel_kinds
      USE stel_constants
      USE v3_utilities
#if defined (MPI_OPT)
      USE mpi
#endif
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  model module parameters
!*******************************************************************************
!>  Maximum size for parameter profile name lengths.
      INTEGER, PARAMETER, PUBLIC :: p_type_len=20

!>  Lower array bound for function profiles.
      INTEGER, PARAMETER, PUBLIC :: ilb_b   = 0
!>  Upper array bound for function profiles.
      INTEGER, PARAMETER, PUBLIC :: iub_b   = 21
!>  Array size for spline profiles.
      INTEGER, PARAMETER, PUBLIC :: iub_asf = 101

!>  No profile type.
      INTEGER, PARAMETER, PRIVATE :: pprofile_none_type         = -1
!>  Two Power profile type.
      INTEGER, PARAMETER, PRIVATE :: pprofile_two_power_type    = 0
!>  Two Power with guassian profile type.
      INTEGER, PARAMETER, PRIVATE :: pprofile_two_power_gs_type = 1
!>  Power Series profile type.
      INTEGER, PARAMETER, PRIVATE :: pprofile_power_series_type = 2
!>  Cubic Spline profile type.
      INTEGER, PARAMETER, PRIVATE :: pprofile_cubic_spline_type = 3
!>  Akima Spline profile type.
      INTEGER, PARAMETER, PRIVATE :: pprofile_akima_spline_type = 4
!>  Line Segment profile type.
      INTEGER, PARAMETER, PRIVATE :: pprofile_line_segment_type = 5

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) model class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a parameterized profile.
!-------------------------------------------------------------------------------
      TYPE pprofile_class
!>  Parameterized profile type
         INTEGER :: p_type = pprofile_none_type
!>  Array of coefficients for the functional profiles.
         REAL(rprec), DIMENSION(ilb_b:iub_b) :: b = 0.0
!>  Array of radial position values for the segmented profiles.
         REAL(rprec), DIMENSION(iub_asf)     :: as = 0.0
!>  Array of profile values for the segmented profiles.
         REAL(rprec), DIMENSION(iub_asf)     :: af = 0.0
!>  Index of the last segemented value.
         INTEGER                             :: maxSplineIndex = 1
      END TYPE

!-------------------------------------------------------------------------------
!>  Pointer to a pprofile object. Used for creating arrays of pprofile pointers.
!>  This is needed because fortran does not allow arrays of pointers directly.
!-------------------------------------------------------------------------------
      TYPE pprofile_pointer
         TYPE (pprofile_class), POINTER :: p => null()
      END TYPE

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref model_class types using
!>  @ref model_construct_new or @ref signal_construct_mpi
!-------------------------------------------------------------------------------
      INTERFACE pprofile_construct
         MODULE PROCEDURE pprofile_construct_new,                              &
     &                    pprofile_construct_mpi
      END INTERFACE

      PRIVATE findMaxIndex

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref pprofile_class.
!>
!>  Allocates memory and initializes a @ref pprofile_class.
!>
!>  @param[in] p_type Profile type discription.
!>  @param[in] b_arg  Array of functional profile coefficients
!>  @param[in] as_arg Array of segmented radial positions.
!>  @param[in] af_arg Array of segmented profile values.
!>  @returns A pointer to a constructed @ref pprofile_class object.
!-------------------------------------------------------------------------------
      FUNCTION pprofile_construct_new(p_type_arg, b_arg, as_arg, af_arg)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile_class), POINTER :: pprofile_construct_new
      CHARACTER (len=*), INTENT(in)         :: p_type_arg
      REAL(rprec), DIMENSION(:), INTENT(in) :: b_arg
      REAL(rprec), DIMENSION(:), INTENT(in) :: as_arg
      REAL(rprec), DIMENSION(:), INTENT(in) :: af_arg

!  local variables
      CHARACTER (len = p_type_len)          :: p_type_lc
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(pprofile_construct_new)

! Store profile coeffients.
      pprofile_construct_new%b = b_arg

! Find and store the maximum s index if specifying splines.
      pprofile_construct_new%maxSplineIndex = findMaxIndex(as_arg)
      pprofile_construct_new%as = as_arg
      pprofile_construct_new%af = af_arg

!  Convert the type to lower case, and truncate
      p_type_lc = p_type_arg
      CALL tolower(p_type_lc)
      SELECT CASE(TRIM(p_type_lc))
      
         CASE ('two_power')
            pprofile_construct_new%p_type = pprofile_two_power_type

         CASE ('two_power_gs')
            pprofile_construct_new%p_type = pprofile_two_power_gs_type

         CASE ('cubic_spline')
            IF (pprofile_construct_new%maxSplineIndex .lt. 4) THEN
               WRITE(*,*) 'pprofile:cubic spline: too few as values'
               WRITE(*,*) 'maxSplineIndex, as = ',                             &
     &                    pprofile_construct_new%maxSplineIndex
               WRITE(*,*) pprofile_construct_new%as
               CALL EXIT(1)
            ENDIF
            pprofile_construct_new%p_type = pprofile_cubic_spline_type

         CASE ('akima_spline')
            IF (pprofile_construct_new%maxSplineIndex .lt. 4) THEN
               WRITE(*,*) 'pprofile:akima spline: too few as values'
               WRITE(*,*) 'maxSplineIndex, as = ',                             &
     &                    pprofile_construct_new%maxSplineIndex
               WRITE(*,*) pprofile_construct_new%as
               CALL EXIT(1)
            ENDIF
            pprofile_construct_new%p_type = pprofile_akima_spline_type

         CASE ('line_segment')
            pprofile_construct_new%p_type = pprofile_line_segment_type

         CASE ('none')
            pprofile_construct_new%p_type = pprofile_none_type

         CASE DEFAULT
            pprofile_construct_new%p_type = pprofile_none_type
            WRITE(*,*) 'Unrecognized p_type:', p_type_lc
            WRITE(*,*) ' *** CHECK YOUR INPUT ***'
         
      END SELECT

      CALL profiler_set_stop_time('pprofile_construct_new', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref pprofile_class object by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref pprofile_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref pprofile_class object.
!-------------------------------------------------------------------------------
      FUNCTION pprofile_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (pprofile_class), POINTER :: pprofile_construct_mpi
      INTEGER, INTENT(in)            :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                        :: error

!  Start of executable code
      ALLOCATE(pprofile_construct_mpi)

      CALL MPI_BCAST(pprofile_construct_mpi%p_type, 1, MPI_INTEGER, 0,         &
     &               intra_comm, error)
      CALL MPI_BCAST(pprofile_construct_mpi%b,                                 &
     &               SIZE(pprofile_construct_mpi%b), MPI_REAL8, 0,             &
     &               intra_comm, error)
      CALL MPI_BCAST(pprofile_construct_mpi%as, iub_asf, MPI_REAL8, 0,         &
     &               intra_comm, error)
      CALL MPI_BCAST(pprofile_construct_mpi%af, iub_asf, MPI_REAL8, 0,         &
     &               intra_comm, error)

      pprofile_construct_mpi%maxSplineIndex =                                  &
     &   findMaxIndex(pprofile_construct_mpi%as)
#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref pprofile_class object.
!>
!>  Deallocates memory and uninitializes a @ref pprofile_class object.
!>
!>  @param[inout] this A @ref pprofile_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_destruct(this)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile_class), POINTER :: this

!  Start of executable code
      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Gets the value of a profile at a radial position.
!>
!>  Evaluates the profile at the specifed radial position. If the radial
!>  position is out side the plasma, it returns the baseline value. If there is
!>  no profile type specified return the default value.
!>
!>  @param[in] this  A @ref pprofile_class instance.
!>  @param[in] s_arg Radial position to evaluate the profile at.
!>  @returns The value of the profile.
!>  @todo This function makes vmec specific assumptions of 1 being the maximum
!>  radial value.
!-------------------------------------------------------------------------------
      FUNCTION pprofile_get_value(this, s_arg)
      USE line_segment
      USE functions

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec)                        :: pprofile_get_value
      TYPE (pprofile_class), INTENT (in) :: this
      REAL(rprec), INTENT(in)            :: s_arg

!  local variables
      REAL(rprec)                        :: s_use
      REAL(rprec)                        :: s_01
      LOGICAL                            :: l_01
      INTEGER                            :: i, iflag
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Define argument and logical variable for s between zero and one
      s_01 = MAX(zero, MIN(one,s_arg))
      l_01 = (s_arg .ge. zero) .and. (s_arg .le. one)

      SELECT CASE(this%p_type)
      
         CASE (pprofile_none_type)
            pprofile_get_value = -1.E49_rprec

! v3fit uses an extra varible as an offset than vmec as a result the b array is
! one index shorter than what two_power takes.
         CASE (pprofile_two_power_type)
            pprofile_get_value = this%b(0)
            IF (l_01) THEN
               pprofile_get_value = pprofile_get_value                         &
     &                            + two_power(s_arg, this%b(1:iub_b))
            ENDIF

! v3fit uses an extra varible as an offset than vmec as a result the b array is
! one index shorter than what two_power_gs takes. As a result only 5 guassian
! peaks maybe specified instead of the 6 vmec can use.
         CASE (pprofile_two_power_gs_type)
            pprofile_get_value = this%b(0)
            IF (l_01) THEN
               pprofile_get_value = pprofile_get_value                         &
     &                            + two_power_gs(s_arg,                        &
     &                                           this%b(1:iub_b))
            ENDIF
      
         CASE (pprofile_power_series_type)
            pprofile_get_value = zero
            DO i = iub_b, ilb_b, -1      ! all of b array, backwards
               pprofile_get_value = pprofile_get_value*s_01                    &
     &                            + this%b(i)
            END DO

         CASE (pprofile_cubic_spline_type)
            s_use = MIN(this%as(this%maxSplineIndex - 1),                      &
     &                  MAX(s_arg, this%as(1)))
            CALL spline_cubic(s_use, pprofile_get_value,                       &
     &                        this%as, this%af,                                &
     &                        this%maxSplineIndex, iflag)
            IF (iflag .ge. 0) THEN

               CALL profiler_set_stop_time('pprofile_get_value',               &
     &                                     start_time)

               RETURN
            ELSE IF (iflag .eq. -1) THEN
               WRITE (*,*) 'ERROR: pprofile: outside value from ' //           &
     &                     'spline_cubic'
            ELSE IF(iflag .eq. -2) THEN
               WRITE (*,*) 'ERROR: pprofile:  decreasing s values ' //         &
     &                     'in spline_cubic'
            ELSE
               WRITE (*,*) 'ERROR: pprofile: unknown error from ' //           &
     &                     'spline_cubic'
            ENDIF
!  Splines shuld not reach this point without triggering an error.
            CALL EXIT(1)

         CASE (pprofile_akima_spline_type)
            s_use = MIN(this%as(this%maxSplineIndex - 1),                      &
     &                  MAX(s_arg,this%as(1)))
            CALL spline_akima(s_use, pprofile_get_value,                       &
     &                        this%as, this%af,                                &
     &                        this%maxSplineIndex, iflag)
            IF (iflag < 0) THEN
               WRITE (*,*) 'ERROR: pprofile: bad value from ' //               &
     &                     'spline_akima requested'
               CALL EXIT(1)
            ENDIF

         CASE (pprofile_line_segment_type)
            CALL line_seg(s_arg, pprofile_get_value, this%as, this%af,         &
     &                    this%maxSplineIndex)
         
      END SELECT

      CALL profiler_set_stop_time('pprofile_get_value', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the name of the profile type.
!>
!>  Convertes a profile type to a string discription.
!>
!>  @param[in] this  A @ref pprofile_class instance.
!>  @returns The name of the profile type.
!-------------------------------------------------------------------------------
      FUNCTION pprofile_get_p_type_name(this)

!  Declare Arguments
      CHARACTER(len=p_type_len)         :: pprofile_get_p_type_name
      TYPE (pprofile_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (this%p_type)

         CASE (pprofile_two_power_type)
            pprofile_get_p_type_name = 'two_power'

         CASE (pprofile_two_power_gs_type)
            pprofile_get_p_type_name = 'two_power_gs'

         CASE (pprofile_power_series_type)
            pprofile_get_p_type_name = 'power_series'

         CASE (pprofile_cubic_spline_type)
            pprofile_get_p_type_name = 'cubic_spline'

         CASE (pprofile_akima_spline_type)
            pprofile_get_p_type_name = 'akima_spline'

         CASE (pprofile_line_segment_type)
            pprofile_get_p_type_name = 'line_segment'

         CASE DEFAULT
            pprofile_get_p_type_name = 'none'

      END SELECT

      CALL profiler_set_stop_time('pprofile_get_p_type_name',                  &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out the profile to an output file.
!>
!>  This formats and writes out information about a parameterized profile.
!>
!>  @param[in] this A @ref pprofile_class instance.
!>  @param[in] iou  Input/output unit of the output file.
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_write(this, id, iou)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (pprofile_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)     :: id
      INTEGER, INTENT(in)               :: iou

!  local variables
      INTEGER                           :: i
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Actual write. Could do a select-case on p_type.
      WRITE(iou,1100) id
      WRITE(iou,1200) TRIM(pprofile_get_p_type_name(this))
      
      SELECT CASE(this%p_type)
         
         CASE (pprofile_two_power_type)
            WRITE(iou,1210)
            WRITE(iou,1300) this%b(0:3)

         CASE (pprofile_power_series_type,                                     &
     &         pprofile_two_power_gs_type)
            WRITE(iou,1220)
            WRITE(iou,1300) this%b(ilb_b:iub_b)
            
         CASE (pprofile_cubic_spline_type,                                     &
     &         pprofile_akima_spline_type,                                     &
     &         pprofile_line_segment_type)
            WRITE(iou,1230)
            WRITE(iou,1231) (i, this%as(i), this%af(i),                        &
     &                       i = 1, this%maxSplineIndex)
         
      END SELECT

      CALL profiler_set_stop_time('pprofile_write', start_time)

1100  FORMAT(/' Parameterized Profile Write: id = ',a)
1200  FORMAT('   pp_type = ',a)
1210  FORMAT(' b_0 + Th(s)Th(1-s)(b_1 (1 - s ** b_2) ** b_3).',                &
     &   '   b(0:3) = ')
1220  FORMAT(' Th(s)Th(1-s)[Sum_0_n b_i s** i].   b(0:n) = ')
1230  FORMAT(' i       as(i)           af(i)')
1231  FORMAT(1x,i3,2x,es15.8,2x,es15.8)
1300  FORMAT(4(2x,es15.8))

      END SUBROUTINE

!*******************************************************************************
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize data to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref pprofile_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (pprofile_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                  :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                              :: error
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL MPI_BCAST(this%p_type, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%b, SIZE(this%b), MPI_REAL8, 0,                &
     &               intra_comm, error)
      CALL MPI_BCAST(this%as, iub_asf, MPI_REAL8, 0, intra_comm, error)
      CALL MPI_BCAST(this%af, iub_asf, MPI_REAL8, 0, intra_comm, error)

      this%maxSplineIndex = findMaxIndex(this%as)

      CALL profiler_set_stop_time('pprofile_sync_children', start_time)
#endif

      END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Finds the index of the maximum radial position.
!>
!>  This assumes that the radials positions are specifed in increasing order.
!>
!>  @param[in] s_array Array of radial positions
!>  @returns Index of the last position.
!-------------------------------------------------------------------------------
      FUNCTION findMaxIndex(s_array)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                               :: findMaxIndex
      REAL(rprec), DIMENSION(:), INTENT(in) :: s_array

!  local variables
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      findMaxIndex = MINLOC(s_array(2:), dim=1)

      CALL profiler_set_stop_time('findMaxIndex', start_time)

      END FUNCTION

      END MODULE
