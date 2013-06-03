!*******************************************************************************
!  File vmec_equilibrium.f
!  Contains module vmec_equilibrium
!  Defines the base class of the type vmec:
!
!  Super Class: equilibrium
!  Sub Classes: None
!
!*******************************************************************************

      MODULE vmec_equilibrium
      USE stel_kinds
      USE data_parameters
      USE pprofile_T

      IMPLICIT NONE

!*******************************************************************************
!  vmec equilibrium module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: vmec_profile_type_length = 30

!  Define id's only for values that can change. All others shouldn't be needed
!  outside of the vmec interface. Parameter id's start at 3 since 0-2 are
!  reserved for non equilibrium model parameters. These numbers will need to be
!  update if any new model parameters are added.

!  Variable Parameters
      INTEGER, PARAMETER :: vmec_ac_id          =  3
      INTEGER, PARAMETER :: vmec_ac_aux_s_id    =  4
      INTEGER, PARAMETER :: vmec_ac_aux_f_id    =  5
      INTEGER, PARAMETER :: vmec_ai_id          =  6
      INTEGER, PARAMETER :: vmec_ai_aux_s_id    =  7
      INTEGER, PARAMETER :: vmec_ai_aux_f_id    =  8
      INTEGER, PARAMETER :: vmec_am_id          =  9
      INTEGER, PARAMETER :: vmec_am_aux_s_id    = 10
      INTEGER, PARAMETER :: vmec_am_aux_f_id    = 11
      INTEGER, PARAMETER :: vmec_bloat_id       = 12
      INTEGER, PARAMETER :: vmec_rbc_id         = 13
      INTEGER, PARAMETER :: vmec_zbs_id         = 14
      INTEGER, PARAMETER :: vmec_rbs_id         = 15
      INTEGER, PARAMETER :: vmec_zbc_id         = 16
      INTEGER, PARAMETER :: vmec_extcur_id      = 17
      INTEGER, PARAMETER :: vmec_curtor_id      = 18
      INTEGER, PARAMETER :: vmec_phiedge_id     = 19
      INTEGER, PARAMETER :: vmec_pres_scale_id  = 20

!  Auxilary Parameters
      INTEGER, PARAMETER :: vmec_rmnc_id        = 21
      INTEGER, PARAMETER :: vmec_zmns_id        = 22
      INTEGER, PARAMETER :: vmec_lmns_id        = 23
      INTEGER, PARAMETER :: vmec_gmnc_id        = 24
      INTEGER, PARAMETER :: vmec_bsubumnc_id    = 25
      INTEGER, PARAMETER :: vmec_bsubvmnc_id    = 26
      INTEGER, PARAMETER :: vmec_bsupsmns_id    = 27
      INTEGER, PARAMETER :: vmec_bsupumnc_id    = 28
      INTEGER, PARAMETER :: vmec_bsupvmnc_id    = 29
      INTEGER, PARAMETER :: vmec_rmns_id        = 30
      INTEGER, PARAMETER :: vmec_zmnc_id        = 31
      INTEGER, PARAMETER :: vmec_lmnc_id        = 32
      INTEGER, PARAMETER :: vmec_gmns_id        = 33
      INTEGER, PARAMETER :: vmec_bsubumns_id    = 34
      INTEGER, PARAMETER :: vmec_bsubvmns_id    = 35
      INTEGER, PARAMETER :: vmec_bsupsmnc_id    = 36
      INTEGER, PARAMETER :: vmec_bsupumns_id    = 37
      INTEGER, PARAMETER :: vmec_bsupvmns_id    = 38
      INTEGER, PARAMETER :: vmec_phi_id         = 39
      INTEGER, PARAMETER :: vmec_iotaf_id       = 40
      INTEGER, PARAMETER :: vmec_rzl_id         = 41

!  Auxilary Parameters 2
      INTEGER, PARAMETER :: vmec_bsubu_id       = 42
      INTEGER, PARAMETER :: vmec_bsubv_id       = 43
      INTEGER, PARAMETER :: vmec_rsuv_id        = 44
      INTEGER, PARAMETER :: vmec_zsuv_id        = 45
      INTEGER, PARAMETER :: vmec_gsuv_id        = 46
      INTEGER, PARAMETER :: vmec_currusuv_id    = 47
      INTEGER, PARAMETER :: vmec_currvsuv_id    = 48
      INTEGER, PARAMETER :: vmec_rusuv_id       = 49
      INTEGER, PARAMETER :: vmec_zusuv_id       = 50
      INTEGER, PARAMETER :: vmec_rvsuv_id       = 51
      INTEGER, PARAMETER :: vmec_zvsuv_id       = 52
      INTEGER, PARAMETER :: vmec_phitot_id      = 53

!  Auxilary Model Parameters
      INTEGER, PARAMETER :: vmec_pp_ne_b_id     = 54
      INTEGER, PARAMETER :: vmec_pp_ne_as_id    = 55
      INTEGER, PARAMETER :: vmec_pp_ne_af_id    = 56
      INTEGER, PARAMETER :: vmec_pp_sxrem_b_id  = 57
      INTEGER, PARAMETER :: vmec_pp_sxrem_as_id = 58
      INTEGER, PARAMETER :: vmec_pp_sxrem_af_id = 59
      INTEGER, PARAMETER :: vmec_pp_te_b_id     = 60
      INTEGER, PARAMETER :: vmec_pp_te_as_id    = 61
      INTEGER, PARAMETER :: vmec_pp_te_af_id    = 62

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) vmec base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Define the vmec base class
!
!  needs_reconvergence  External vmec only needs to be reconverged after a vmec
!                       reconstruction parameter has been changed.
!
!  ne     Parameterized profile describing a flux surface constant electron
!         denisty profile.
!
!  te     Parameterized profile describing a flux surface constant electron
!         temperature profile.
!
!  sxrem  Parameterized profile describing a flux surface conatsnt soft x-ray
!         emissivity profile.
!-------------------------------------------------------------------------------
      TYPE vmec_class
         CHARACTER (len=data_filename_length) :: wout_file_name
         INTEGER                  :: wout_file_name_end_index
         TYPE (pprofile), POINTER :: ne => null()
         TYPE (pprofile), POINTER :: te => null()
         TYPE (pprofile), POINTER :: sxrem => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION vmec_construct()

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER :: vmec_construct

!  Start of executable code
      ALLOCATE(vmec_construct)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE vmec_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%ne)) THEN
         CALL pprofile_destroy(this%ne)
         this%ne => null()
      END IF

      IF (ASSOCIATED(this%te)) THEN
         CALL pprofile_destroy(this%te)
         this%te => null()
      END IF

      IF (ASSOCIATED(this%sxrem)) THEN
         CALL pprofile_destroy(this%sxrem)
         this%sxrem => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
      FUNCTION vmec_set_param(this, id, i_index, j_index, value)
      USE vmec_input, only: rbc, zbs, rbs, zbc,                                &
     &                      ac, ac_aux_s, ac_aux_f,                            &
     &                      ai, ai_aux_s, ai_aux_f,                            &
     &                      am, am_aux_s, am_aux_f,                            &
     &                      extcur, curtor, phiedge, pres_scale,               &
     &                      bloat
      IMPLICIT NONE

!  Declare Sruments
      LOGICAL                          :: vmec_set_param
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(in)              :: id
      INTEGER, INTENT(in)              :: i_index
      INTEGER, INTENT(in)              :: j_index
      REAL (rprec), INTENT(in)         :: value

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_rbc_id)
            vmec_set_param = .true.
            rbc(i_index, j_index) = value

         CASE (vmec_zbs_id)
            vmec_set_param = .true.
            zbs(i_index, j_index) = value

         CASE (vmec_rbs_id)
            vmec_set_param = .true.
            rbs(i_index, j_index) = value

         CASE (vmec_zbc_id)
            vmec_set_param = .true.
            zbc(i_index, j_index) = value

         CASE (vmec_ac_id)
            vmec_set_param = .true.
            ac(i_index) = value

         CASE (vmec_ac_aux_s_id)
            vmec_set_param = .true.
            ac_aux_s(i_index) = value

         CASE (vmec_ac_aux_f_id)
            vmec_set_param = .true.
            ac_aux_f(i_index) = value

         CASE (vmec_ai_id)
            vmec_set_param = .true.
            ai(i_index) = value

         CASE (vmec_ai_aux_s_id)
            vmec_set_param = .true.
            ai_aux_s(i_index) = value

         CASE (vmec_ai_aux_f_id)
            vmec_set_param = .true.
            ai_aux_f(i_index) = value

         CASE (vmec_am_id)
            vmec_set_param = .true.
            am(i_index) = value

         CASE (vmec_am_aux_s_id)
            vmec_set_param = .true.
            am_aux_s(i_index) = value

         CASE (vmec_am_aux_f_id)
            vmec_set_param = .true.
            am_aux_f(i_index) = value

         CASE (vmec_pp_ne_b_id)
            vmec_set_param = .false.
            this%ne%b(i_index) = value

         CASE (vmec_pp_ne_as_id)
            vmec_set_param = .false.
            this%ne%as(i_index) = value

         CASE (vmec_pp_ne_af_id)
            vmec_set_param = .false.
            this%ne%af(i_index) = value

         CASE (vmec_pp_te_b_id)
            vmec_set_param = .false.
            this%te%b(i_index) = value

         CASE (vmec_pp_te_as_id)
            vmec_set_param = .false.
            this%te%as(i_index) = value

         CASE (vmec_pp_te_af_id)
            vmec_set_param = .false.
            this%te%af(i_index) = value

         CASE (vmec_pp_sxrem_b_id)
            vmec_set_param = .false.
            this%sxrem%b(i_index) = value

         CASE (vmec_pp_sxrem_as_id)
            vmec_set_param = .false.
            this%sxrem%as(i_index) = value

         CASE (vmec_pp_sxrem_af_id)
            vmec_set_param = .false.
            this%sxrem%af(i_index) = value

         CASE (vmec_extcur_id)
            vmec_set_param = .true.
            extcur(i_index) = value

         CASE (vmec_curtor_id)
            vmec_set_param = .true.
            curtor = value

         CASE (vmec_phiedge_id)
            vmec_set_param = .true.
            phiedge = value

         CASE (vmec_pres_scale_id)
            vmec_set_param = .true.
            pres_scale = value

         CASE (vmec_bloat_id)
            vmec_set_param = .true.
            bloat = value

         CASE DEFAULT
            vmec_set_param = .true.

      END SELECT

      END FUNCTION

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the id number for a vmec parameter
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_param_id(this, param_name)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: vmec_get_param_id
      TYPE (vmec_class), INTENT(in)                :: this
      CHARACTER (len=data_name_length), INTENT(in) :: param_name

!  Start of executable code
      SELECT CASE (TRIM(ADJUSTL(param_name)))

         CASE ('ac')
            vmec_get_param_id = vmec_ac_id

         CASE ('ac_aux_s')
            vmec_get_param_id = vmec_ac_aux_s_id

         CASE ('ac_aux_f')
            vmec_get_param_id = vmec_ac_aux_f_id

         CASE ('ai')
            vmec_get_param_id = vmec_ai_id

         CASE ('ai_aux_s')
            vmec_get_param_id = vmec_ai_aux_s_id

         CASE ('ai_aux_f')
            vmec_get_param_id = vmec_ai_aux_f_id

         CASE ('am')
            vmec_get_param_id = vmec_am_id

         CASE ('am_aux_s')
            vmec_get_param_id = vmec_am_aux_s_id

         CASE ('am_aux_f')
            vmec_get_param_id = vmec_am_aux_f_id

         CASE ('bloat')
            vmec_get_param_id = vmec_bloat_id

         CASE ('rbc')
            vmec_get_param_id = vmec_rbc_id

         CASE ('zbs')
            vmec_get_param_id = vmec_zbs_id

         CASE ('rbs')
            vmec_get_param_id = vmec_rbs_id

         CASE ('zbc')
            vmec_get_param_id = vmec_zbc_id

         CASE ('extcur')
            vmec_get_param_id = vmec_extcur_id

         CASE ('curtor')
            vmec_get_param_id = vmec_curtor_id

         CASE ('phiedge')
            vmec_get_param_id = vmec_phiedge_id

         CASE ('pres_scale')
            vmec_get_param_id = vmec_pres_scale_id

         CASE ('rmnc')
            vmec_get_param_id = vmec_rmnc_id

         CASE ('zmns')
            vmec_get_param_id = vmec_zmns_id

         CASE ('lmns')
            vmec_get_param_id = vmec_lmns_id

         CASE ('gmnc')
            vmec_get_param_id = vmec_gmnc_id

         CASE ('bsubumnc')
            vmec_get_param_id = vmec_bsubumnc_id

         CASE ('bsubvmnc')
            vmec_get_param_id = vmec_bsubvmnc_id

         CASE ('bsupsmns')
            vmec_get_param_id = vmec_bsupsmns_id

         CASE ('bsupumnc')
            vmec_get_param_id = vmec_bsupumnc_id

         CASE ('bsupvmnc')
            vmec_get_param_id = vmec_bsupvmnc_id

         CASE ('rmns')
            vmec_get_param_id = vmec_rmns_id

         CASE ('zmnc')
            vmec_get_param_id = vmec_zmnc_id

         CASE ('lmnc')
            vmec_get_param_id = vmec_lmnc_id

         CASE ('gmns')
            vmec_get_param_id = vmec_gmns_id

         CASE ('bsubumns')
            vmec_get_param_id = vmec_bsubumns_id

         CASE ('bsubvmns')
            vmec_get_param_id = vmec_bsubvmns_id

         CASE ('bsupsmnc')
            vmec_get_param_id = vmec_bsupsmnc_id

         CASE ('bsupumns')
            vmec_get_param_id = vmec_bsupumns_id

         CASE ('bsupvmns')
            vmec_get_param_id = vmec_bsupvmns_id

         CASE ('phi')
            vmec_get_param_id = vmec_phi_id

         CASE ('iotaf')
            vmec_get_param_id = vmec_iotaf_id

         CASE ('rzl')
            vmec_get_param_id = vmec_rzl_id

         CASE ('bsubu')
            vmec_get_param_id = vmec_bsubu_id

         CASE ('bsubv')
            vmec_get_param_id = vmec_bsubv_id

         CASE ('rsuv')
            vmec_get_param_id = vmec_rsuv_id

         CASE ('zsuv')
            vmec_get_param_id = vmec_zsuv_id

         CASE ('gsuv')
            vmec_get_param_id = vmec_gsuv_id

         CASE ('currusuv')
            vmec_get_param_id = vmec_currusuv_id

         CASE ('currvsuv')
            vmec_get_param_id = vmec_currvsuv_id

         CASE ('rusuv')
            vmec_get_param_id = vmec_rusuv_id

         CASE ('zusuv')
            vmec_get_param_id = vmec_zusuv_id

         CASE ('rvsuv')
            vmec_get_param_id = vmec_rvsuv_id

         CASE ('zvsuv')
            vmec_get_param_id = vmec_zvsuv_id

         CASE ('phitot')
            vmec_get_param_id = vmec_phitot_id

         CASE ('pp_ne_b')
            vmec_get_param_id = vmec_pp_ne_b_id

         CASE ('pp_ne_as')
            vmec_get_param_id = vmec_pp_ne_as_id

         CASE ('pp_ne_af')
            vmec_get_param_id = vmec_pp_ne_af_id

         CASE ('pp_sxrem_b')
            vmec_get_param_id = vmec_pp_sxrem_b_id

         CASE ('pp_sxrem_as')
            vmec_get_param_id = vmec_pp_sxrem_as_id

         CASE ('pp_sxrem_af')
            vmec_get_param_id = vmec_pp_sxrem_af_id

         CASE ('pp_te_b')
            vmec_get_param_id = vmec_pp_te_b_id

         CASE ('pp_te_as')
            vmec_get_param_id = vmec_pp_te_as_id

         CASE ('pp_te_af')
            vmec_get_param_id = vmec_pp_te_af_id

         CASE DEFAULT
            WRITE(*,*) TRIM(ADJUSTL(param_name)) //                            &
     &                 ' is not a valid vmec parameter'
            STOP

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!   Get Scaler Real value for an id number
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_scaler_value(this, id)
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: vmec_get_scaler_value
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      IF (.not. vmec_is_scaler_value(this, id)) THEN
         WRITE(*,*) vmec_get_param_name(this, id) //                           &
     &              ' Not scaler type.'
         STOP
      END IF

      vmec_get_scaler_value = 0.0 ! Set to zero for now.

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get vmec parameter name for an id number
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_param_name(this, id)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=data_name_length) :: vmec_get_param_name
      TYPE (vmec_class), INTENT(in)   :: this
      INTEGER, INTENT(in)             :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_ac_id)
            vmec_get_param_name = 'ac'

         CASE (vmec_ac_aux_s_id)
            vmec_get_param_name = 'ac_aux_s'

         CASE (vmec_ac_aux_f_id)
            vmec_get_param_name = 'ac_aux_f'

         CASE (vmec_ai_id)
            vmec_get_param_name = 'ai'

         CASE (vmec_ai_aux_s_id)
            vmec_get_param_name = 'ai_aux_s'

         CASE (vmec_ai_aux_f_id)
            vmec_get_param_name = 'ai_aux_f'

         CASE (vmec_am_id)
            vmec_get_param_name = 'am'

         CASE (vmec_am_aux_s_id)
            vmec_get_param_name = 'am_aux_s'

         CASE (vmec_am_aux_f_id)
            vmec_get_param_name = 'am_aux_f'

         CASE (vmec_bloat_id)
            vmec_get_param_name = 'bloat'

         CASE (vmec_rbc_id)
            vmec_get_param_name = 'rbc'

         CASE (vmec_zbs_id)
            vmec_get_param_name = 'zbs'

         CASE (vmec_rbs_id)
            vmec_get_param_name = 'rbs'

         CASE (vmec_zbc_id)
            vmec_get_param_name = 'zbc'

         CASE (vmec_extcur_id)
            vmec_get_param_name = 'extcur'

         CASE (vmec_curtor_id)
            vmec_get_param_name = 'curtor'

         CASE (vmec_phiedge_id)
            vmec_get_param_name = 'phiedge'

         CASE (vmec_pres_scale_id)
            vmec_get_param_name = 'pres_scale'

         CASE (vmec_rmnc_id)
            vmec_get_param_name = 'rmnc'

         CASE (vmec_zmns_id)
            vmec_get_param_name = 'zmns'

         CASE (vmec_lmns_id)
            vmec_get_param_name = 'lmns'

         CASE (vmec_gmnc_id)
            vmec_get_param_name = 'gmnc'

         CASE (vmec_bsubumnc_id)
            vmec_get_param_name = 'bsubumnc'

         CASE (vmec_bsubvmnc_id)
            vmec_get_param_name = 'bsubvmnc'

         CASE (vmec_bsupsmns_id)
            vmec_get_param_name = 'bsupsmns'

         CASE (vmec_bsupumnc_id)
            vmec_get_param_name = 'bsupumnc'
 
         CASE (vmec_bsupvmnc_id)
            vmec_get_param_name = 'bsupvmnc'

         CASE (vmec_rmns_id)
            vmec_get_param_name = 'rmns'

         CASE (vmec_zmnc_id)
            vmec_get_param_name = 'zmnc'

         CASE (vmec_lmnc_id)
            vmec_get_param_name = 'lmnc'

         CASE (vmec_gmns_id)
            vmec_get_param_name = 'gmns'
 
         CASE (vmec_bsubumns_id)
            vmec_get_param_name = 'bsubumns'

         CASE (vmec_bsubvmns_id)
            vmec_get_param_name = 'bsubvmns'

         CASE (vmec_bsupsmnc_id)
            vmec_get_param_name = 'bsupsmnc'

         CASE (vmec_bsupumns_id)
            vmec_get_param_name = 'bsupumns'

         CASE (vmec_bsupvmns_id)
            vmec_get_param_name = 'bsupvmns'

         CASE (vmec_phi_id)
            vmec_get_param_name = 'phi'

         CASE (vmec_iotaf_id)
            vmec_get_param_name = 'iotaf'

         CASE (vmec_rzl_id)
            vmec_get_param_name = 'rzl'

         CASE (vmec_bsubu_id)
            vmec_get_param_name = 'bsubu'

         CASE (vmec_bsubv_id)
            vmec_get_param_name = 'bsubv'

         CASE (vmec_rsuv_id)
            vmec_get_param_name = 'rsuv'

         CASE (vmec_zsuv_id)
            vmec_get_param_name = 'zsuv'

         CASE (vmec_gsuv_id)
            vmec_get_param_name = 'gsuv'

         CASE (vmec_currusuv_id)
            vmec_get_param_name = 'currusuv'

         CASE (vmec_currvsuv_id)
            vmec_get_param_name = 'currvsuv'

         CASE (vmec_rusuv_id)
            vmec_get_param_name = 'rusuv'
   
         CASE (vmec_zusuv_id)
            vmec_get_param_name = 'zusuv'

         CASE (vmec_rvsuv_id)
            vmec_get_param_name = 'rvsuv'

         CASE (vmec_zvsuv_id)
            vmec_get_param_name = 'zvsuv'

         CASE (vmec_phitot_id)
            vmec_get_param_name = 'phitot'

         CASE (vmec_pp_ne_b_id)
            vmec_get_param_name = 'pp_ne_b'

         CASE (vmec_pp_ne_as_id)
            vmec_get_param_name = 'pp_ne_as'

         CASE (vmec_pp_ne_af_id)
            vmec_get_param_name = 'pp_ne_af'

         CASE (vmec_pp_sxrem_b_id)
            vmec_get_param_name = 'pp_sxrem_b'

         CASE (vmec_pp_sxrem_as_id)
            vmec_get_param_name = 'pp_sxrem_as'

         CASE (vmec_pp_sxrem_af_id)
            vmec_get_param_name = 'pp_sxrem_af'

         CASE (vmec_pp_te_b_id)
            vmec_get_param_name = 'pp_te_b'

         CASE (vmec_pp_te_as_id)
            vmec_get_param_name = 'pp_te_as'

         CASE (vmec_pp_te_af_id)
            vmec_get_param_name = 'pp_te_af'

         CASE DEFAULT
            STOP 'Invalid VMEC parameter id.'

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get vmec electron density at position
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_ne(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_ne
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      IF (ASSOCIATED(this%ne)) THEN
         vmec_get_ne = pprofile_evaluate(this%ne,                              &
     &                                   vmec_get_s(this, x_cart))
      ELSE
         vmec_get_ne = 0.0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get vmec electrocn temperature at position
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_te(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_te
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      IF (ASSOCIATED(this%te)) THEN
         vmec_get_te = pprofile_evaluate(this%te,                              &
     &                                   vmec_get_s(this, x_cart))
      ELSE
         vmec_get_te = 0.0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get vmec soft x-ray emissivity at position
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_sxrem(this, x_cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_sxrem
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Start of executable code
      IF (ASSOCIATED(this%te)) THEN
         vmec_get_sxrem = pprofile_evaluate(this%sxrem,                        &
     &                                      vmec_get_s(this, x_cart))
      ELSE
         vmec_get_sxrem = 0.0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get vmec plasma pressure at position
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_p(this, x_cart)
      USE stel_constants, only: mu0

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_p
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  Forward declare the interface for pmass. This explicitly defined the type for
!  pmass to avoid a compiler error. If profile_functions.f was implimented as a
!  module this would not be needed.
      INTERFACE
         FUNCTION pmass(s)
         USE stel_kinds
         REAL (rprec)             :: pmass
         REAL (rprec), INTENT(in) :: s
         END FUNCTION
      END INTERFACE

!  Start of executable code
      vmec_get_p = pmass(vmec_get_s(this, x_cart))/mu0

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get vmec b field at position
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_B_vec(this, x_cart)
      USE coordinate_utilities, ONLY : cyl_to_cart_vec, cart_to_cyl
      USE vmec_utils, ONLY : GetBcyl_WOUT

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: vmec_get_B_vec
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  local variables
      REAL (rprec), DIMENSION(3)             :: r_cyl

!  Start of executable code
      r_cyl = cart_to_cyl(x_cart)

      CALL GetBcyl_WOUT(r_cyl(1), r_cyl(2), r_cyl(3),                          &
     &                  vmec_get_B_vec(1), vmec_get_B_vec(2),                  &
     &                  vmec_get_B_vec(3))

      vmec_get_B_vec = cyl_to_cart_vec(r_cyl, vmec_get_B_vec)

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get vmec line integrated b at position
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_Int_B_dphi(this, r, theta)
      USE line_segment, only: line_seg
      USE stel_constants, only: twopi
      USE read_wout_mod, only: bsubvmnc, ns

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: vmec_get_Int_B_dphi
      TYPE (vmec_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)      :: r
      REAL (rprec), INTENT(in)      :: theta

!  local variables
      REAL (rprec)                  :: bsubv00c, ds
      INTEGER                       :: i
      REAL (rprec), DIMENSION(ns)   :: s

!  Start of executable code
      ds = 1.0/(ns - 1)
      DO i = 0, ns - 1
         s(i + 1) = i*ds
      END DO

      CALL line_seg(r, bsubv00c, s, bsubvmnc(1,:), ns)

      vmec_get_Int_B_dphi = twopi*bsubv00c

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get vmec s value at position
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_s(this, x_cart)
      USE cyl_flux, only: cyl2flx
      USE coordinate_utilities, only: cart_to_cyl
      USE read_wout_mod, only: rzl_local, ns, ntmax, lthreed
      USE vmec_input, only: ntor, mpol, nfp, lasym

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: vmec_get_s
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart

!  local variables
      REAL (rprec), DIMENSION(3)             :: r_cyl
      REAL (rprec), DIMENSION(3)             :: s_flux = 0.0
      REAL (rprec)                           :: fmin = 0.0
      INTEGER                                :: info = 0
      INTEGER                                :: nfe = 0

!  Start of executable code
      r_cyl = cart_to_cyl(x_cart)
      r_cyl(2) = nfp*r_cyl(2)

      CALL cyl2flx(rzl_local, r_cyl, s_flux, ns, ntor, mpol, ntmax,            &
     &             lthreed, lasym, info, nfe, fmin)

      vmec_get_s = s_flux(1)

      IF (vmec_get_s .le. 1.0) THEN
         IF (fmin .gt.  1.0E-8_rprec) THEN
            WRITE(*,*) "Warn: VMEC could not converge to flux position."       &
            WRITE(*,*) info, fmin, vmec_get_s
         END IF
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get plasma edge at phi position
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_plasma_edge(this, phi, r, z)
      USE vmec_input, only: mpol
      USE read_wout_mod, only: xm, xn, rmnc, zmns, rmns, zmnc, lasym, ns
      USE stel_constants, only: twopi

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: vmec_get_plasma_edge
      TYPE (vmec_class), INTENT(in)                    :: this
      REAL (rprec), INTENT (in)                        :: phi
      REAL (rprec), DIMENSION(:), POINTER, INTENT(out) :: r
      REAL (rprec), DIMENSION(:), POINTER, INTENT(out) :: z

!  local variables
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: cosv, sinv
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: cosu, sinu
      INTEGER                                 :: i

!  Start of executable code
      vmec_get_plasma_edge = 2*mpol + 6

      IF (.not.ASSOCIATED(r)) THEN
         ALLOCATE(r(vmec_get_plasma_edge), z(vmec_get_plasma_edge))
      END IF
      ALLOCATE(cosv(SIZE(xn, 1)), sinv(SIZE(xn, 1)))
      ALLOCATE(cosu(SIZE(xm, 1)), sinu(SIZE(xm, 1)))

      cosv(:) = COS(xn*phi)
      sinv(:) = SIN(xn*phi)
      DO i = 1, vmec_get_plasma_edge
         cosu(:) = COS(xm*(i - 1)*twopi/vmec_get_plasma_edge)
         sinu(:) = SIN(xm*(i - 1)*twopi/vmec_get_plasma_edge)

         r(i) = SUM(rmnc(:, ns)*(cosu(:)*cosv(:) + sinu(:)*sinv(:)))
         z(i) = SUM(zmns(:, ns)*(sinu(:)*cosv(:) - cosu(:)*sinv(:)))

         IF (lasym) THEN
            r(i) = r(i) +                                                      &
     &             SUM(rmns(:, ns)*(sinu(:)*cosv(:) - cosu(:)*sinv(:)))
            z(i) = z(i) +                                                      &
     &             SUM(zmnc(:, ns)*(cosu(:)*cosv(:) + sinu(:)*sinv(:)))
         END IF
      END DO

      DEALLOCATE(cosu, cosv, sinu, sinv)

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Query variable id for scaler type
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_scaler_value(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_scaler_value
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_bloat_id, vmec_curtor_id, vmec_phiedge_id,                 &
     &         vmec_pres_scale_id, vmec_phitot_id)
            vmec_is_scaler_value = .true.

         CASE DEFAULT
            vmec_is_scaler_value = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 1D array
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_1d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_1d_array
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_ac_id, vmec_ac_aux_s_id, vmec_ac_aux_f_id,                 &
     &         vmec_ai_id, vmec_ai_aux_s_id, vmec_ai_aux_f_id,                 &
     &         vmec_am_id, vmec_am_aux_s_id, vmec_am_aux_f_id,                 &
     &         vmec_pp_ne_b_id, vmec_pp_ne_as_id, vmec_pp_ne_af_id,            &
     &         vmec_pp_te_b_id, vmec_pp_te_as_id, vmec_pp_te_af_id,            &
     &         vmec_pp_sxrem_b_id, vmec_pp_sxrem_as_id,                        &
     &         vmec_pp_sxrem_af_id, vmec_extcur_id, vmec_phi_id,               &
     &         vmec_iotaf_id)
            vmec_is_1d_array = .true.

         CASE DEFAULT
            vmec_is_1d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 2D array
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_2d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_2d_array
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_rbc_id, vmec_zbs_id, vmec_rbs_id, vmec_zbc_id,             &
     &         vmec_rmnc_id, vmec_zmns_id, vmec_lmns_id, vmec_gmnc_id,         &
     &         vmec_bsubumnc_id, vmec_bsubvmnc_id,                             &
     &         vmec_bsupsmns_id, vmec_bsupumnc_id, vmec_bsupvmnc_id,           &
     &         vmec_rmns_id, vmec_zmnc_id, vmec_lmnc_id, vmec_gmns_id,         &
     &         vmec_bsubumns_id, vmec_bsubvmns_id,                             &
     &         vmec_bsupsmnc_id, vmec_bsupumns_id, vmec_bsupvmns_id)
            vmec_is_2d_array = .true.

         CASE DEFAULT
            vmec_is_2d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 3D array
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_3d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_3d_array
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_bsubu_id, vmec_bsubv_id,                                   &
     &         vmec_rsuv_id, vmec_zsuv_id, vmec_gsuv_id,                       &
     &         vmec_currusuv_id, vmec_currvsuv_id, vmec_rusuv_id,              &
     &         vmec_zusuv_id, vmec_rvsuv_id, vmec_zvsuv_id)
            vmec_is_3d_array = .true.

         CASE DEFAULT
            vmec_is_3d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id for 4D array
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_4d_array(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_4d_array
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_rzl_id)
            vmec_is_4d_array = .true.

         CASE DEFAULT
            vmec_is_4d_array = .false.

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query variable id is a valid reconstruction parameter
!-------------------------------------------------------------------------------
      FUNCTION vmec_is_recon_param(this, id)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_is_recon_param
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(in)           :: id

!  Start of executable code
      SELECT CASE (id)

         CASE (vmec_rbc_id, vmec_zbs_id, vmec_rbs_id, vmec_zbc_id,             &
     &         vmec_ac_id, vmec_ac_aux_s_id, vmec_ac_aux_f_id,                 &
     &         vmec_ai_id, vmec_ai_aux_s_id, vmec_ai_aux_f_id,                 &
     &         vmec_am_id, vmec_am_aux_s_id, vmec_am_aux_f_id,                 &
     &         vmec_pp_ne_b_id, vmec_pp_ne_as_id, vmec_pp_ne_af_id,            &
     &         vmec_pp_te_b_id, vmec_pp_te_as_id, vmec_pp_te_af_id,            &
     &         vmec_pp_sxrem_b_id, vmec_pp_sxrem_as_id,                        &
     &         vmec_pp_sxrem_af_id, vmec_extcur_id, vmec_curtor_id,            &
     &         vmec_phiedge_id, vmec_pres_scale_id, vmec_bloat_id)
            vmec_is_recon_param = .true.

         CASE DEFAULT
            vmec_is_recon_param = .false.

      END SELECT

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
      FUNCTION vmec_converge(this, num_iter)
      USE vmec_input, only: delt, niter
      USE vmec_params, only: timestep_flag, output_flag,                       &
     &                       reset_jacdt_flag, jac75_flag,                     &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag
      USE xstuff, only: xcdot
      USE v3f_global, only: iou_runlog
      USE read_wout_mod, only: read_wout_file, LoadRZL
      USE vmec_main, only: iterc

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: vmec_converge
      TYPE (vmec_class), INTENT(in) :: this
      INTEGER, INTENT(inout)        :: num_iter

!  local variables
!  delt_local   Value of delt on entry. Store, for restoration after error
!               recovery.
!  ictrl_array  Array containing various vmec control parameters.
      REAL(rprec)                   :: delt_local
      INTEGER, DIMENSION(5)         :: ictrl_array = 0
      INTEGER                       :: recovery_index = 0
      INTEGER                       :: error, num_iter_init

!  Start of executable code
      delt_local = delt
      num_iter_init = iterc

!  Run vmec until successful or an unrecoverable error has occured.
      DO

         ictrl_array(1) = timestep_flag + output_flag +                        &
     &                    reset_jacdt_flag
         ictrl_array(3) = niter
         ictrl_array(4) = 1

!  Zero the xc arrays.
         xcdot = zero
         CALL runvmec(ictrl_array, this%wout_file_name, .false.)

!  Set error flag. (Flags declared in vmec_params.f)
!
!  Check error flags for error recovery modes.
!  If a recoverable error has occured, rerun vmec with altered eq_state. Error
!  recovery is attempted a maximum of 2 times.
!  If a non-recoverable error has occured, do not rerun vmec and report error.
!  If vmec terminated sucessfully continue.
         SELECT CASE (ictrl_array(2))
            CASE (jac75_flag)
               IF (recovery_index .le. 2) THEN
                  WRITE (*,*) 'VMEC CONVERGENCE ERROR: jac75_flag'
                  WRITE (*,*) 'delt chnage: old, new:', delt, delt/2.0
                  WRITE (iou_runlog,*)                                         &
     &               'VMEC CONVERGENCE ERROR: jac75_flag'
                  WRITE (iou_runlog,*)                                         &
     &               'delt chnage: old, new:', delt, delt/2.0
                  delt = delt/2.0
               ELSE
                  EXIT
               END IF

            CASE DEFAULT
               EXIT

         END SELECT

         WRITE (*,*) 'Attempting error recovery'
         WRITE (iou_runlog,*) 'Attempting error recovery'

         recovery_index = recovery_index + 1
      END DO

!  Check for errors, if vmec was successful, load vmec data for use with the
!  signals. Other wise report the error.
      vmec_converge = .false.
      SELECT CASE (ictrl_array(2))

         CASE (successful_term_flag, norm_term_flag, more_iter_flag)
            CALL read_wout_file(this%wout_file_name, error)
            IF (error .eq. 0) THEN
               CALL LoadRZL
               vmec_converge = .true.
            END IF

         CASE DEFAULT
!  FIXME: Put in error reporting here. For now leave this blank until I can
!  figure how v3fit decides if vmec has converged.

      END SELECT

!  Increment the number of iterations by the number of iterations used by
!  runvmec.
      num_iter = num_iter + iterc - num_iter_init
      delt = delt_local

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************

      END MODULE
