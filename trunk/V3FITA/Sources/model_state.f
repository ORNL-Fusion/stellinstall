!*******************************************************************************
!>  @file model_state.f
!>  @brief Contains module @ref model_state.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Contains parameters defining the bit positions for flags that mark changes
!>  in different parts of the model.
!*******************************************************************************
      MODULE model_state

      IMPLICIT NONE

!*******************************************************************************
!  module parameters
!*******************************************************************************
!>  Set all flags off.
      INTEGER, PARAMETER :: model_state_all_off = 0
!>  Set all flags on.
      INTEGER, PARAMETER ::                                                    &
     &   model_state_all_on = NOT(model_state_all_off)

!>  Equilibrium changed bit position.
      INTEGER, PARAMETER :: model_state_eq_flag    = 0
!>  Denisty profile changed bit position.
      INTEGER, PARAMETER :: model_state_ne_flag    = 1
!>  Temperature profile changed bit position.
      INTEGER, PARAMETER :: model_state_te_flag    = 2
!>  Shift parameter changed bit position.
      INTEGER, PARAMETER :: model_state_shift_flag = 3
!>  Soft x-ray emissivity profile changed bit position. This flag needs to
!>  always be the last flag an can be no larger than 31 minus the maximum number
!>  of possible sxrem profiles @ref v3fit_input::v3fit_max_sxrem_profiles. If
!>  the number of model states exceeds 31, then the size of
!>  @ref model::model_class%state_flags needs to be increased.
      INTEGER, PARAMETER :: model_state_sxrem_flag = 4

      END MODULE
