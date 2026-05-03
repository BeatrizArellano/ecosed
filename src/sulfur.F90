#include "fabm_driver.h"

! ---------------------------------------------------------------------
! Sulfur redox cycle
!
! Minimal representation:
!   - so4     : dissolved sulfate, SO4--
!   - sulfide : total dissolved sulfide, ΣH2S = H2S + HS- + S--
!
! Sulfide oxidation written as:
!   H2S + 2 O2 + 2 HCO3- -> SO4-- + 2 CO2 + 2 H2O
!
! Model representation:
!   - Sulfide is consumed and SO4 is produced (1:1 in S units)
!   - O2 is consumed at 2 mol O2 per mol sulfide oxidised
!   - Alkalinity decreases by 2 eq per mol sulfide oxidised
!   - No explicit DIC change: HCO3- is converted to CO2
!
! Sulfate reduction by organic matter is handled in the OM degradation
! module, because it is an organic matter remineralisation pathway.
! ---------------------------------------------------------------------
module sulfur

   use fabm_types
   use molecular_diff, only: DIFF_ION_LINEAR, m0_SO4, m1_SO4, m0_HS, m1_HS
   implicit none
   private

   type, extends(type_base_model), public :: type_sulfur

      ! --- State variables
      type(type_state_variable_id) :: id_so4
      type(type_state_variable_id) :: id_sulfide

      ! --- Optional couplings
      type(type_state_variable_id) :: id_o2
      type(type_state_variable_id) :: id_alk

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_sulfide_ox

      ! --- Parameters
      real(rk) :: k_sulfide_ox      ! Maximum sulfide oxidation rate (s-1 internally)
      real(rk) :: o2_sulfide_ox_sat ! O2 where sulfide oxidation is fully active

   contains
      procedure :: initialize
      procedure :: do
   end type type_sulfur

contains

   subroutine initialize(self, configunit)
      class(type_sulfur), intent(inout), target :: self
      integer,            intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
      real(rk), parameter :: eps     = 1.0e-12_rk

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%k_sulfide_ox, 'k_sulfide_ox', 'd-1', 'Maximum first-order sulfide oxidation rate', default=5.0_rk, scale_factor=d_per_s, minimum=0.0_rk)
      call self%get_parameter(self%o2_sulfide_ox_sat, 'o2_sulfide_ox_sat', 'mmol m-3', 'O2 concentration above which sulfide oxidation is fully active', default=0.1_rk, minimum=eps)

      ! ---------------- State variables ----------------
      call self%register_state_variable(self%id_so4, 'so4', 'mmol m-3', &
                                        'Dissolved sulfate', initial_value=28000.0_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_so4, 'is_solute', .true.)
      call self%set_variable_property(self%id_so4, 'diff_method', DIFF_ION_LINEAR)
      call self%set_variable_property(self%id_so4, 'm0', m0_SO4)
      call self%set_variable_property(self%id_so4, 'm1', m1_SO4)

      call self%register_state_variable(self%id_sulfide, 'sulfide', 'mmol m-3', &
                                        'Total dissolved sulfide', initial_value=0.0_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_sulfide, 'is_solute', .true.)
      call self%set_variable_property(self%id_sulfide, 'diff_method', DIFF_ION_LINEAR)
      call self%set_variable_property(self%id_sulfide, 'm0', m0_HS)
      call self%set_variable_property(self%id_sulfide, 'm1', m1_HS)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_o2, 'o2', 'mmol m-3', 'Dissolved oxygen', required=.true.)
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.false.)

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_sulfide_ox, 'sulfide_ox', 'mmol m-3 d-1', 'Sulfide oxidation rate')

   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)
      class(type_sulfur), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sulfide, o2
      real(rk) :: fO2, xO2
      real(rk) :: sulfide_ox
      real(rk) :: o2_cons_sulfide_ox
      real(rk) :: alk_change

      real(rk), parameter :: secs_per_day        = 86400.0_rk

      real(rk), parameter :: o2_per_sulfide_ox  = 2.0_rk
      real(rk), parameter :: alk_per_sulfide_ox = -2.0_rk

      _LOOP_BEGIN_

         _GET_(self%id_sulfide, sulfide)
         _GET_(self%id_o2, o2)

         ! ------------------------------------------------------------------
         ! Sulfide oxidation
         !
         ! H2S + 2O2 + 2HCO3- -> SO4-- + 2CO2 + 2H2O
         !
         ! Consumes dissolved sulfide and O2, produces dissolved sulfate.
         ! Alkalinity decreases by 2 equivalents per mol sulfide oxidised.
         ! No explicit DIC change: HCO3- is converted to CO2.
         !
         ! The sulfide state variable represents total dissolved sulfide,
         ! ΣH2S = H2S + HS- + S--. Its molecular diffusivity is approximated
         ! using HS-, because HS- dominates dissolved sulfide at seawater pH.
         ! ------------------------------------------------------------------

         ! Sulfide oxidation is represented as first-order in sulfide, with a capped
         ! maximum rate. Oxygen acts as an activation factor via a smoothstep
         ! function: oxidation is zero at O2 = 0, increases smoothly at low O2,
         ! and reaches full activity above a small O2 threshold.
         !
         ! This replaces the standard mass-action formulation (rate ∝ sulfide × O2),
         ! which can lead to very fast, numerically stiff reactions under oxic
         ! conditions. The present formulation preserves O2 control while avoiding
         ! excessive rates and improving numerical stability.
         xO2 = max(0.0_rk, min(1.0_rk, max(o2, 0.0_rk) / self%o2_sulfide_ox_sat))
         fO2 = xO2*xO2*xO2 * (xO2 * (6.0_rk*xO2 - 15.0_rk) + 10.0_rk)

         sulfide_ox = self%k_sulfide_ox * fO2 * max(sulfide, 0.0_rk)

         o2_cons_sulfide_ox = o2_per_sulfide_ox  * sulfide_ox
         alk_change         = alk_per_sulfide_ox * sulfide_ox

         _ADD_SOURCE_(self%id_sulfide, -sulfide_ox)
         _ADD_SOURCE_(self%id_so4,      sulfide_ox)

         _ADD_SOURCE_(self%id_o2, -o2_cons_sulfide_ox)

         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_change)

         _SET_DIAGNOSTIC_(self%id_sulfide_ox, sulfide_ox * secs_per_day)

      _LOOP_END_

   end subroutine do

end module sulfur