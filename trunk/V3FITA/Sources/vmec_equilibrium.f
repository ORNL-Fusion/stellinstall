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

!  Auxiliary Parameters
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

!  Auxiliary Model Parameters
      INTEGER, PARAMETER :: vmec_pp_ne_b_id     = 42
      INTEGER, PARAMETER :: vmec_pp_ne_as_id    = 43
      INTEGER, PARAMETER :: vmec_pp_ne_af_id    = 44
      INTEGER, PARAMETER :: vmec_pp_sxrem_b_id  = 45
      INTEGER, PARAMETER :: vmec_pp_sxrem_as_id = 46
      INTEGER, PARAMETER :: vmec_pp_sxrem_af_id = 47
      INTEGER, PARAMETER :: vmec_pp_te_b_id     = 48
      INTEGER, PARAMETER :: vmec_pp_te_as_id    = 49
      INTEGER, PARAMETER :: vmec_pp_te_af_id    = 50

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) vmec magnetic cache
!  2) vmec base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Define an object to cache data needed to compute magnetic signals.
!-------------------------------------------------------------------------------
      TYPE vmec_magnetic_cache
!  Number of toroidal planes.
         INTEGER                                 :: kv

!  Volume differential elements
         REAL (rprec)                            :: ds
         REAL (rprec)                            :: du
         REAL (rprec)                            :: dv

!  R and Z positions on the suv gird points.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: rsuv => null()
         REAL (rprec), DIMENSION(:,:,:), POINTER :: zsuv => null()

!  The following have a 1/Sqrt(g) factor embeded.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jrsuv => null()
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jphisuv => null()
         REAL (rprec), DIMENSION(:,:,:), POINTER :: jzsuv => null()
      END TYPE

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
         CHARACTER (len=data_filename_length) :: vmec_file_name
         TYPE (pprofile), POINTER             :: ne => null()
         TYPE (pprofile), POINTER             :: te => null()
         TYPE (pprofile), POINTER             :: sxrem => null()
         TYPE (vmec_magnetic_cache), POINTER  ::                               &
     &      magnetic_cache => null()
         REAL (rprec), DIMENSION(:), POINTER  :: xc => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE vmec_set_magnetic_cache
         MODULE PROCEDURE vmec_set_magnetic_cache_init,                        &
     &                    vmec_set_magnetic_cache_calc
      END INTERFACE

!  Forward declare the interface for runvmec. This explicitly defined the type
!  for runvmec to avoid a compiler bugs. If runvmec.f was implimented as a
!  module this would not be needed.
      INTERFACE
         SUBROUTINE runvmec(control_array, vmec_input_file,                    &
     &                      use_screen, reset_file_name)
         INTEGER, INTENT(inout), TARGET :: control_array(5)
         LOGICAL, INTENT(in)            :: use_screen
         CHARACTER(LEN=*), INTENT(in)   :: vmec_input_file
         CHARACTER(LEN=*), OPTIONAL     :: reset_file_name
         END SUBROUTINE
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION vmec_construct(file_name)
      USE vmec_params, only: restart_flag, readin_flag, timestep_flag,         &
     &                       successful_term_flag, norm_term_flag,             &
     &                       more_iter_flag
      USE vmec_input, only: l_v3fit

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), POINTER    :: vmec_construct
      CHARACTER (len=*), INTENT(in) :: file_name

!  local variables
      INTEGER, DIMENSION(5)         :: ictrl_array = 0
      INTEGER                       :: index_dat, index_end

!  Start of executable code
      WRITE (*,*) ' *** Initializing VMEC equilibrium from file ' //           &
     &            TRIM(file_name)

      ALLOCATE(vmec_construct)

!  Initialize VMEC by taking a single step.
      ictrl_array(1) = restart_flag + readin_flag + timestep_flag
      ictrl_array(3) = 1
      ictrl_array(4) = 1 ! Limit VMEC to single grid.

      l_v3fit = .true.

      vmec_construct%vmec_file_name = TRIM(file_name)
      CALL runvmec(ictrl_array, file_name, .false.)

!  Check for errors in runvmec.
      SELECT CASE (ictrl_array(2))

         CASE (successful_term_flag, norm_term_flag, more_iter_flag)
!  Set the wout file name.
            index_dat = INDEX(file_name,'input.')
            index_end = LEN_TRIM(file_name)
            IF (index_dat .gt. 0) THEN
               vmec_construct%wout_file_name = 'wout_' //                      &
     &            TRIM(file_name(index_dat + 6:index_end)) // '.nc'
            ELSE
               vmec_construct%wout_file_name = 'wout_' //                      &
     &            TRIM(file_name(1:index_end)) // '.nc'
            END IF

         CASE DEFAULT
            CALL err_fatal('vmec_construct runvmec command failed ',           &
     &                     int = ictrl_array(2))

      END SELECT

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

      IF (ASSOCIATED(this%magnetic_cache)) THEN
         DEALLOCATE(this%magnetic_cache%rsuv)
         DEALLOCATE(this%magnetic_cache%zsuv)
         DEALLOCATE(this%magnetic_cache%jrsuv)
         DEALLOCATE(this%magnetic_cache%jphisuv)
         DEALLOCATE(this%magnetic_cache%jzsuv)

         this%magnetic_cache%rsuv => null()
         this%magnetic_cache%zsuv => null()
         this%magnetic_cache%jrsuv => null()
         this%magnetic_cache%jphisuv => null()
         this%magnetic_cache%jzsuv => null()
         this%magnetic_cache%kv = 0

         DEALLOCATE(this%magnetic_cache)
         this%magnetic_cache => null()
      END IF

      IF (ASSOCIATED(this%xc)) THEN
         DEALLOCATE(this%xc)
         this%xc => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Set reconstruction parameter
!-------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------
!  Set magnetic cache initialy
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_magnetic_cache_init(this, total_tor_planes,          &
     &                                        tor_planes_used)
      USE vmec_input, only: nfp, ns_array
      USE stel_constants, only: twopi

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(in)              :: total_tor_planes
      INTEGER, INTENT(in)              :: tor_planes_used

!  local arguments
      INTEGER                          :: error

!  Start of executable code.
      ALLOCATE(this%magnetic_cache)

!  Set the differental volume elements.
      this%magnetic_cache%ds = 1.0/(ns_array(1) - 1.0)
      this%magnetic_cache%du = twopi/ns_array(1)
      this%magnetic_cache%dv = (twopi/nfp)/total_tor_planes

!  When stellarator symmetry is used, do not need to store all the toroidal
!  planes in a single field period.
      this%magnetic_cache%kv = tor_planes_used

      ALLOCATE(this%magnetic_cache%rsuv(ns_array(1), ns_array(1),              &
     &                                  tor_planes_used))
      ALLOCATE(this%magnetic_cache%zsuv(ns_array(1), ns_array(1),              &
     &                                  tor_planes_used))
      ALLOCATE(this%magnetic_cache%jrsuv(ns_array(1), ns_array(1),             &
     &                                   tor_planes_used))
      ALLOCATE(this%magnetic_cache%jphisuv(ns_array(1), ns_array(1),           &
     &                                     tor_planes_used))
      ALLOCATE(this%magnetic_cache%jzsuv(ns_array(1), ns_array(1),             &
     &                                   tor_planes_used))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Set magnetic cache
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_set_magnetic_cache_calc(this)
      USE read_wout_mod, only: lasym, xm, xn, mnmax, xm_nyq, xn_nyq,           &
     &                         mnmax_nyq, rmnc, zmns, rmns, zmnc,              &
     &                         currumnc, currvmnc, currumns, currvmns,         &
     &                         ns

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this

!  local variables
      INTEGER                                   :: s, u, v
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosz, sinz
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: cosu, cosv
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: sinu, sinv
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: cosz_nyq, sinz_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: cosu_nyq, cosv_nyq
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: sinu_nyq, sinv_nyq
      REAL (rprec)                              :: ru, rv, zu, zv
      REAL (rprec)                              :: ju, jv

!  Start of executable code.
      ALLOCATE(cosz(mnmax), sinz(mnmax))
      ALLOCATE(cosu(mnmax, ns), sinu(mnmax, ns),                               &
     &         cosv(mnmax, this%magnetic_cache%kv),                            &
     &         sinv(mnmax, this%magnetic_cache%kv))

      ALLOCATE(cosz_nyq(mnmax_nyq), sinz_nyq(mnmax_nyq))
      ALLOCATE(cosu_nyq(mnmax_nyq, ns), sinu_nyq(mnmax_nyq, ns),               &
     &         cosv_nyq(mnmax_nyq, this%magnetic_cache%kv),                    &
     &         sinv_nyq(mnmax_nyq, this%magnetic_cache%kv))

! U grid is square with the s grid.
      DO u = 1, ns
         cosu(:,u) = COS(xm*(u - 1)*this%magnetic_cache%du)
         sinu(:,u) = SIN(xm*(u - 1)*this%magnetic_cache%du)
         cosu_nyq(:,u) = COS(xm_nyq*(u - 1)*this%magnetic_cache%du)
         sinu_nyq(:,u) = SIN(xm_nyq*(u - 1)*this%magnetic_cache%du)
      END DO
      DO v = 1, this%magnetic_cache%kv
         cosv(:,v) = COS(xn*(v - 1)*this%magnetic_cache%dv)
         sinv(:,v) = SIN(xn*(v - 1)*this%magnetic_cache%dv)
         cosv_nyq(:,v) = COS(xn_nyq*(v - 1)*this%magnetic_cache%dv)
         sinv_nyq(:,v) = SIN(xn_nyq*(v - 1)*this%magnetic_cache%dv)
      END DO

!  Set the suv grids
!
!  ru = dR/du
!  rv = dRdv   Note: dRHat/dv = PhiHat
!
!  zu = dZ/du
!  zu = dZ/dv
!
!  jr = ju*ru + jv*rv
!  jphi = jv*R
!  jz = ju*zu + jv*zv
      DO u = 1, ns

         DO v = 1, this%magnetic_cache%kv

            cosz = cosu(:,u)*cosv(:,v) + sinu(:,u)*sinv(:,v)
            sinz = sinu(:,u)*cosv(:,v) - cosu(:,u)*sinv(:,v)

            cosz_nyq = cosu_nyq(:,u)*cosv_nyq(:,v)                             &
     &               + sinu_nyq(:,u)*sinv_nyq(:,v)
            sinz_nyq = sinu_nyq(:,u)*cosv_nyq(:,v)                             &
     &               - cosu_nyq(:,u)*sinv_nyq(:,v)

            DO s = 1, ns

               this%magnetic_cache%rsuv(s,u,v) = SUM(rmnc(:,s)*cosz)
               this%magnetic_cache%zsuv(s,u,v) = SUM(zmns(:,s)*sinz)
               ru = -SUM(xm*rmnc(:,s)*sinz)
               zu = SUM(xm*zmns(:,s)*cosz)
               rv = SUM(xn*rmnc(:,s)*sinz)
               zv = -SUM(xn*zmns(:,s)*cosz)

               ju = SUM(currumnc(:,s)*cosz_nyq)
               jv = SUM(currvmnc(:,s)*cosz_nyq)

               IF (lasym) THEN
                  this%magnetic_cache%rsuv(s, u, v) =                          &
     &               this%magnetic_cache%rsuv(s, u, v) +                       &
     &               SUM(rmns(:,s)*sinz)
                  this%magnetic_cache%zsuv(s, u, v) =                          &
     &               this%magnetic_cache%zsuv(s, u, v) +                       &
     &               SUM(zmnc(:,s)*cosz)
                  ru = ru + SUM(xm*rmns(:,s)*cosz)
                  zu = zu - SUM(xm*zmnc(:,s)*sinz)
                  rv = rv - SUM(xn*rmns(:,s)*cosz)
                  zv = zv + SUM(xn*zmnc(:,s)*sinz)

                  ju = ju + SUM(currumnc(:,s)*sinz_nyq)
                  jv = jv + SUM(currvmnc(:,s)*sinz_nyq)
               END IF

               this%magnetic_cache%jrsuv(s,u,v) = ju*ru + jv*rv
               this%magnetic_cache%jphisuv(s,u,v) =                            &
     &            jv*this%magnetic_cache%rsuv(s, u, v)
               this%magnetic_cache%jzsuv(s,u,v) = ju*zu + jv*zv
            END DO
         END DO
      END DO

      DEALLOCATE(cosz, sinz)
      DEALLOCATE(cosu, sinu, cosv, sinv)
      DEALLOCATE(cosz_nyq, sinz_nyq)
      DEALLOCATE(cosu_nyq, sinu_nyq, cosv_nyq, sinv_nyq)

      END SUBROUTINE

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

!-------------------------------------------------------------------------------
!  Get plasma only magnetic signal.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_plasma_mag_signal(this, r_fun, phi_fun, z_fun,         &
     &                                    r_max, r_min, num_r,                 &
     &                                    z_max, z_min, num_z,                 &
     &                                    is_stell_sym, num_phi)
      USE read_wout_mod, only: isigng, ns
      USE vmec_input, only: nfp
      USE bivariate
      USE stel_constants, only: twopi

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: vmec_get_plasma_mag_signal
      TYPE (vmec_class), INTENT(in)              :: this
      REAL (rprec), DIMENSION(:,:,:), INTENT(in) :: r_fun
      REAL (rprec), DIMENSION(:,:,:), INTENT(in) :: phi_fun
      REAL (rprec), DIMENSION(:,:,:), INTENT(in) :: z_fun
      REAL (rprec), INTENT(in)                   :: r_max
      REAL (rprec), INTENT(in)                   :: r_min
      INTEGER, INTENT(in)                        :: num_r
      REAL (rprec), INTENT(in)                   :: z_max
      REAL (rprec), INTENT(in)                   :: z_min
      INTEGER, INTENT(in)                        :: num_z
      LOGICAL, INTENT(in)                        :: is_stell_sym
      INTEGER, INTENT(in)                        :: num_phi

!  local variables
      TYPE (bivariate_type)                      :: bivariate_object
      INTEGER                                    :: v, r, z
      REAL (rprec), DIMENSION(num_r)             :: rgrid
      REAL (rprec), DIMENSION(num_z)             :: zgrid
      REAL (rprec), DIMENSION(ns, ns)            :: responce_rsuv
      REAL (rprec), DIMENSION(ns, ns)            :: responce_phisuv
      REAL (rprec), DIMENSION(ns, ns)            :: responce_zsuv
      REAL (rprec), DIMENSION(ns, ns)            :: sumsuv
      REAL (rprec), DIMENSION(num_phi)           :: sumv

!  Start of executable code
      DO r = 1, num_r
         rgrid(r) = r_min + (r - 1)*(r_max - r_min)/(num_r - 1)
      END DO
      DO z = 1, num_z
         zgrid(z) = z_min + (z - 1)*(z_max - z_min)/(num_z - 1)
      END DO

      CALL createbivariate(bivariate_object,                                   &
     &                     SIZE(this%magnetic_cache%rsuv, 1),                  &
     &                     SIZE(this%magnetic_cache%rsuv, 2),                  &
     &                     SIZE(this%magnetic_cache%rsuv, 3), 1)

      DO v = 1, num_phi
!  Generate bivariate interpolation responce tables.
         CALL setbivariate(bivariate_object,                                   &
     &                     this%magnetic_cache%rsuv(:,:,v),                    &
     &                     this%magnetic_cache%zsuv(:,:,v),                    &
     &                     rgrid, zgrid, v, 1)

         CALL bivariate4pt(bivariate_object,                                   &
     &                     r_fun, responce_rsuv, v)
         CALL bivariate4pt(bivariate_object,                                   &
     &                     phi_fun, responce_phisuv, v)
         CALL bivariate4pt(bivariate_object,                                   &
     &                     z_fun, responce_zsuv, v)

         sumsuv = this%magnetic_cache%jrsuv(:,:,v)*responce_rsuv               &
     &          + this%magnetic_cache%jphisuv(:,:,v)*responce_phisuv           &
     &          + this%magnetic_cache%jzsuv(:,:,v)*responce_zsuv

         sumv(v) = SUM(sumsuv(2:ns - 1,:))                                     &
     &           + SUM(sumsuv(1,:) + sumsuv(ns,:))/2.0

      END DO

      CALL clear_bivar(bivariate_object)

      vmec_get_plasma_mag_signal = SUM(sumv)
      IF (is_stell_sym) THEN
         vmec_get_plasma_mag_signal = vmec_get_plasma_mag_signal -             &
     &      (sumv(1) + sumv(num_phi))/2.0
      END IF

      vmec_get_plasma_mag_signal = vmec_get_plasma_mag_signal                  &
     &                           * vmec_get_volume_int_element(this)

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get volume integration element.
!
!  For vmec this is J*ds*du*dv where J is the jacobian. The jacobian is
!  integrated into Jr, Jphi and Jz but still need to account for the sign of the
!  jacobian.
!
!  ds = 1/(ns - 1)
!  du = 2*Pi/number of flux surfaces
!  dv = 2*Pi/number of field periods/number of toroidal planes
!
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_volume_int_element(this)
      USE read_wout_mod, only: isigng

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: vmec_get_volume_int_element
      TYPE (vmec_class), INTENT(in) :: this

!  Start of executable code
      vmec_get_volume_int_element = isigng*this%magnetic_cache%ds              &
     &                            * this%magnetic_cache%du                     &
     &                            * this%magnetic_cache%dv

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get external current magnetic signal.
!-------------------------------------------------------------------------------
      FUNCTION vmec_get_ext_mag_signal(this, diag_coils,                       &
     &                                 current_scales, num_coil_groups)
      USE vmec_input, only: lfreeb
      USE mgrid_mod, only: mgrid_mode, nextcur
      USE vmec_input, only: extcur

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: vmec_get_ext_mag_signal
      TYPE (vmec_class), INTENT(in)          :: this
      REAL (rprec), DIMENSION(:), INTENT(in) :: diag_coils
      REAL (rprec), DIMENSION(:), INTENT(in) :: current_scales
      INTEGER, INTENT(in)                    :: num_coil_groups

!  local variables
      INTEGER                       :: num_currents

!  Start of executable code
      IF (lfreeb) THEN
         num_currents = MIN(nextcur, SIZE(diag_coils, 1),                      &
     &                      num_coil_groups)

         SELECT CASE (mgrid_mode)

            CASE ('S', 's', 'N', 'n')
               vmec_get_ext_mag_signal =                                       &
     &            SUM(diag_coils(1:num_currents) *                             &
     &                extcur(1:num_currents) /                                 &
     &                current_scales(1:num_currents))

            CASE DEFAULT
               vmec_get_ext_mag_signal =                                       &
     &            SUM(diag_coils(1:num_currents) *                             &
     &                extcur(1:num_currents))

         END SELECT
      ELSE
         vmec_get_ext_mag_signal = 0.0
      END IF

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
     &         vmec_pres_scale_id)
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
      vmec_is_3d_array = .false.

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
!-------------------------------------------------------------------------------
!  Converge the equilibrium.
!-------------------------------------------------------------------------------
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
      LOGICAL                          :: vmec_converge
      TYPE (vmec_class), INTENT(inout) :: this
      INTEGER, INTENT(inout)           :: num_iter

!  local variables
!  delt_local   Value of delt on entry. Store, for restoration after error
!               recovery.
!  ictrl_array  Array containing various vmec control parameters.
      REAL(rprec)                      :: delt_local
      INTEGER, DIMENSION(5)            :: ictrl_array = 0
      INTEGER                          :: recovery_index = 0
      INTEGER                          :: error, num_iter_init

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
         CALL runvmec(ictrl_array, this%vmec_file_name, .false.)

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

               IF (ASSOCIATED(this%magnetic_cache)) THEN
                  CALL vmec_set_magnetic_cache(this)
               END IF
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

!-------------------------------------------------------------------------------
!  Save the internal state of the equilibrium.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_save_state(this)
      USE xstuff, only: xc

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%xc)) THEN
         IF (SIZE(xc, 1) .eq. SIZE(this%xc, 1)) THEN
            this%xc = xc
         ELSE
            DEALLOCATE(this%xc)
            ALLOCATE(this%xc(SIZE(xc, 1)))
            this%xc = xc
         END IF
      ELSE
         ALLOCATE(this%xc(SIZE(xc, 1)))
         this%xc = xc
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Reset the internal state of the equilibrium.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_reset_state(this)
      USE xstuff, only: xc

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vmec_class), INTENT(in) :: this

!  Start of executable code
      IF (ASSOCIATED(this%xc)) THEN
         IF (SIZE(xc, 1) .eq. SIZE(this%xc, 1)) THEN
            xc = this%xc
         ELSE IF (SIZE(xc, 1) .eq. SIZE(this%xc, 1)) THEN
            xc(1:SIZE(this%xc, 1)) = this%xc
         ELSE
            xc = this%xc(1:SIZE(xc, 1))
         END IF
      END IF

      END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************

      END MODULE
