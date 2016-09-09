!*******************************************************************************
!>  @file ece_travis.f
!>  @brief Contains module @ref ece_travis.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref ece_travis_class.
!>  @par Super Class:
!>  @ref diagnostic
!*******************************************************************************

      MODULE ece_travis

      USE stel_kinds, only: rprec
      USE travis_parameters
      USE model

      IMPLICIT NONE
      

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) ece_travis base class
!  2) ece_travis context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing an ECE signal.
!>  @par Super Class:
!>  @ref diagnostic
!-------------------------------------------------------------------------------
      TYPE ece_travis_class
!>  ID of the Antenna for this ECE signal 
         INTEGER                     :: antenna_id
!>  ECE mode type (0 is O mode a 1 is X mode)
         INTEGER                     :: wmode
!>  ECE Frequency in GHz
         REAL (rprec)                :: freq
!> I need to figure out what this does
         REAL (rprec)                :: band
!> I need to figure out what this does
         INTEGER                     :: disc
      END TYPE ece_travis_class


!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!> @TODO these functions were defined for OLD_TRAVIS
 !     PRIVATE ::  travis_set_tracingoutput, travis_set_VesselMirrors,          &
 !    &            travis_set_Antenna, travis_initEquilibr

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref ece_travis_class object.
!>
!>  Allocates memory and initializes a @ref ece_travis_class object.
!>
!>  @param[in] data_file TRAVIS data file where ece signal is defined.
!>  @param[in] index Index of ECE signal in TRAVIS data file.
!>  @returns A pointer to a constructed @ref ece_travis_class object.
!-------------------------------------------------------------------------------
      FUNCTION ece_travis_construct(antenna_id, wmode,freq, band, disc)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (ece_travis_class), POINTER        :: ece_travis_construct
      INTEGER, INTENT(in)                     :: antenna_id
      INTEGER, INTENT(in)                     :: wmode
      REAL (rprec), INTENT(in)                :: freq
      REAL (rprec), INTENT(in)                :: band
      INTEGER, INTENT(in)                     :: disc
!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(ece_travis_construct)

      ece_travis_construct%antenna_id = antenna_id
      ece_travis_construct%wmode = wmode
      ece_travis_construct%freq = freq
      ece_travis_construct%band = band
      ece_travis_construct%disc = disc

      CALL profiler_set_stop_time('ece_travis_construct', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref ece_travis_class object.
!>
!>  Deallocates memory and uninitializes a @ref ece_class object.
!>
!>  @param[inout] this A @ref ece_travis_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE ece_travis_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (ece_travis_class), POINTER :: this

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  Calculates the ECE signal by searching along the chord for the resonance.
!>  The search criteria is provided by the @ref ece_function.
!>
!>  @param[in]    this    A @ref ece_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION ece_travis_get_modeled_signal(this, a_model, sigma)

#if defined (TRAV)
      USE head
      USE inout_data
      USE travis_inout
#endif
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: ece_travis_get_modeled_signal
      TYPE (ece_travis_class), INTENT(inout)         :: this
      TYPE (model_class), TARGET, INTENT(in)  :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

! local variables
      REAL (rprec)                            :: start_time
      CHARACTER(len=trav_fnl) :: wout_fnm
!> @todo read this variables for a_model
      INTEGER(4) :: n,ns
      REAL (8), ALLOCATABLE :: rho(:), ne(:), Te(:), Zeff(:)
      REAL(8) :: r,dr,df !TODO TEMP?
      CHARACTER(Len=256) :: labelType !TODO TEMP
      CHARACTER(len=5) :: hedi
      REAL(8) :: Tece
      INTEGER(4):: nresmax
      REAL(8) :: tau
      REAL(8) :: rhocy
      REAL(8) :: Rcy
      REAL(8),DIMENSION(3) :: rhoeceb
      REAL(8),DIMENSION(3) :: Receb
      REAL(8),DIMENSION(3) :: ueceb
!> @todo clenan local variables list
      INTEGER(4) :: ier
      INTEGER(8) :: mconf8=0
      INTEGER(8) :: mVessel=0

#if defined (TRAV)
      TYPE(Scenario) :: this_scen
      TYPE(type_EC_Antenna) :: this_antenna
#endif

!  Start of executable code
      start_time = profiler_get_start_time()
      ece_travis_get_modeled_signal = 0.0


#if defined (TRAV)

!> @todo I need to think about parallelization
!
      hedi='ece'
!> Set the travis hearder
      wout_fnm = a_model%equilibrium%vmec%wout_file_name
!> @TODO DEBUG 
      ier = mcLoad(mconf8,wout_fnm)
      write(*,*) "ier = ", ier
      IF (.NOT.(ier==1)) THEN
         write(*,*) "Error in Travis MCLOAD"
         write(*,*) "ier = ", ier
         ece_travis_get_modeled_signal = 0.0
         RETURN
       ENDIF

      
!> Set the TRAVIS OUTPUT (setting these to 0 disables the output)      
      CALL set_TracingOutput_f77(travisOutTemp,travisScrOut) 
!> Set Travis vessel and mirror      
      CALL set_VesselMirrors_f77(Vessel_configuration,
     &                            Mirror_Configuration)
!> Set The Travis Antenna
      CALL set_Antenna_f77(rbeam(this%antenna_id,:),                    &
     &                 rfocus(this%antenna_id,:),                       &
     &                 phibx(this%antenna_id),                          &
     &                 QOemul(this%antenna_id),                         &
     &                 nra(this%antenna_id),                            &
     &                 nphi(this%antenna_id),                           &
     &                 antenna_position(this%antenna_id,:),             &
     &                 target_position(this%antenna_id,:),              &
     &                 antenna_coordinate,                              &
     &                 target_coordinate)
!> Setup the equilibrium in TRAVIS  
      CALL initEquilibr_f77(hgrid,dphi,b0_ref,phi_ref,b_scale,B0_type)

      CALL set_magconfig_f77(mConf8,mVessel)

!> @todo ns should be changed to an input 
      ns =50
      ALLOCATE(rho(ns),ne(ns),Te(ns),Zeff(ns))
      dr = 1.0/(1.0*ns-1.0)
      DO n=1,ns
         r=(n-1)*dr !r/a
         rho(n) = r
         ne(n) = model_get_ne(a_model,rho(n))
         Te(n) = model_get_te(a_model,rho(n))/1.0e3_rprec
!> @todo write statements are for debugging
         write(*,*) "rho    ne     Te(keV)"
         write(*,*) rho(n), ne(n),  te(n)
         Zeff(n)=1.0
      ENDDO

!> @todo make sure v3fit label is toroidal flux
!      labelType = 'tor_rho'
      labelType = 'tor'
!
      CALL set_nTZ_profs_f77(ns,rho,ne,Te,Zeff,labelType)
      
!> @todo Fix to use all of TRAVIS output
! Is there a reason not to use all of the beam info?
! I can make a flag or set to true
!      IF (.false.) THEN
      IF (.true.) THEN
         CALL run_ECE_BEAM_f77c(1,this%freq,this%wmode,Tece,nresmax,           &
     &            tau,Rcy,rhocy,Receb,rhoeceb,ueceb)
         Tece = Tece*1.0e3_rprec ! covernt KeV to eV
!> @todo Debug output all information for debugging
         write(*,*) "After run ece beam"
         write(*,*) "freq, wmode, Tece, nresmax"
         write(*,*) this%freq, this%wmode, Tece, nresmax
         write(*,*) "tau, Rcy, rhocy"
         write(*,*) tau, Rcy, rhocy
      ELSE
         CALL run_ECE_BEAM_f77m(1,this%freq,this%wmode, Tece)
         Tece = Tece*1.0e3_rprec ! covernt KeV to eV
      ENDIF

       DEALLOCATE(rho,ne,Te,Zeff)
       ece_travis_get_modeled_signal = Tece
       
#else
!> Exit if V3FIT is not linked to Travis
         WRITE(*,*) "V3FIT must be linked to the TRAVIS library to use "       &
     &            //"the ece Travis diagnostic! This build of V3FIT is "       &
     &           // " not linked to TRAVIS."
         STOP 
#endif

!> @todo is there a good way to get sigma from Travis?
      sigma = 0.0
      
      CALL profiler_set_stop_time('ece_travis_get_modeled_signal',             &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the ece type.
!>
!>  Returns a description of the ece type for use when writting output files.
!>
!>  @param[in] this A @ref ece_class instance.
!>  @returns A string describing the ece type.
!-------------------------------------------------------------------------------
      FUNCTION ece_travis_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: ece_travis_get_signal_type
      TYPE (ece_travis_class), INTENT(in)     :: this

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ece_travis_get_signal_type = 'ece_travis'

      CALL profiler_set_stop_time('ece_get_signal_type', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  Returns a description of the array indices for use when writting output
!>  files.
!>
!>  @param[in] this A @ref ece_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION ece_travis_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &                                 ece_travis_get_header
      TYPE (ece_travis_class), INTENT(in)                   :: this

!  local variables
      REAL (rprec)                                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!> @todo These hearders are wrong!
      ece_travis_get_header(1) = 'x (m)'
      ece_travis_get_header(2) = 'y (m)'
      ece_travis_get_header(3) = 'z (m)'

      ece_travis_get_header(4) = 'model_sig(1)'
      ece_travis_get_header(5) = 'model_sig(2)'
      ece_travis_get_header(6) = 'model_sig(3)'
      ece_travis_get_header(7) = 'model_sig(4)'

      CALL profiler_set_stop_time('ece_travis_get_header', start_time)

      END FUNCTION

      END MODULE
