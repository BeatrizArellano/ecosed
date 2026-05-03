#include "fabm_driver.h"

! ---------------------------------------------------------------------
! Manganese redox cycle
!
! Minimal representation:
!   - mn2  : dissolved reduced manganese, Mn(II)
!   - mno2 : particulate oxidised manganese, Mn(IV) oxide
!
! Oxidation:
!   Mn2+ + 0.5 O2 + 2 HCO3- -> MnO2(s) + 2 CO2 + H2O
!
! Model representation:
!   - Mn2 is consumed and MnO2 is produced (1:1 in Mn units)
!   - O2 is consumed at 0.5 mol O2 per mol Mn2 oxidised
!   - Alkalinity decreases by 2 eq per mol Mn2 oxidised
!   - MnO2 is particulate and sinks
!
! Mn reduction by organic matter is handled in the OM degradation module,
! because it is an organic matter remineralisation pathway.
! ---------------------------------------------------------------------
module manganese

   use fabm_types
   use molecular_diff, only: DIFF_ION_LINEAR, m0_Mn2, m1_Mn2
   implicit none
   private

   type, extends(type_base_model), public :: type_manganese

      ! --- State variables
      type(type_state_variable_id) :: id_mn2
      type(type_state_variable_id) :: id_mno2

      ! --- Optional couplings
      type(type_state_variable_id) :: id_o2
      type(type_state_variable_id) :: id_alk

      ! --- Dependencies
      type(type_dependency_id)     :: id_porosity

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_mn_ox

      ! --- Parameters
      real(rk) :: k_mn_ox        ! Mn2 oxidation rate (s-1 internally)
      real(rk) :: o2_mn_ox_sat   ! O2 concentration where Mn oxidation is fully active

   contains
      procedure :: initialize
      procedure :: do
   end type type_manganese

contains

   subroutine initialize(self, configunit)
      class(type_manganese), intent(inout), target :: self
      integer,                intent(in)            :: configunit

      real(rk), parameter :: d_per_s     = 1.0_rk / 86400.0_rk
      real(rk), parameter :: m_d_per_m_s = 1.0_rk / 86400.0_rk
      real(rk), parameter :: eps         = 1.0e-12_rk

      real(rk) :: w_mno2       ! MnO2 sinking velocity (m s-1 internally)

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%k_mn_ox, 'k_mn_ox', 'd-1', 'Maximum first-order Mn(II) oxidation rate', default=5.0_rk, scale_factor=d_per_s, minimum=0.0_rk)
      call self%get_parameter(self%o2_mn_ox_sat, 'o2_mn_ox_sat', 'mmol m-3', 'O2 concentration above which Mn(II) oxidation is fully active', default=0.1_rk, minimum=eps)
      call self%get_parameter(w_mno2, 'w_mno2', 'm d-1', 'Sinking velocity of particulate MnO2', default=-1.0_rk, scale_factor=m_d_per_m_s)

      ! ---------------- State variables ----------------
      call self%register_state_variable(self%id_mn2, 'mn2', 'mmol m-3', 'Dissolved Mn(II)', initial_value=0.001_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_mn2, 'is_solute', .true.)
      call self%set_variable_property(self%id_mn2, 'diff_method', DIFF_ION_LINEAR)
      call self%set_variable_property(self%id_mn2, 'm0', m0_Mn2)
      call self%set_variable_property(self%id_mn2, 'm1', m1_Mn2)

      call self%register_state_variable(self%id_mno2, 'mno2', 'mmol m-3', 'Particulate Mn(IV) oxide', initial_value=0.001_rk, minimum=0.0_rk, vertical_movement=w_mno2)
      call self%set_variable_property(self%id_mno2, 'is_solute', .false.)

      ! ---------------- Dependencies ----------------
      call self%register_dependency(self%id_porosity, type_interior_standard_variable(name='porosity', units='1'), required=.false.)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_o2, 'o2', 'mmol m-3', 'Dissolved oxygen', required=.true.)
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.false.)

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_mn_ox, 'mn2_ox', 'mmol m-3 d-1', 'Mn(II) oxidation rate')

   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)
      class(type_manganese), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: mn2, o2
      real(rk) :: fO2, xO2
      real(rk) :: mn_ox
      real(rk) :: o2_cons_mn_ox
      real(rk) :: alk_change
      real(rk) :: phi, phi_s, d2p
      real(rk) :: mno2_prod      

      real(rk), parameter :: secs_per_day   = 86400.0_rk
      real(rk), parameter :: o2_per_mn_ox   = 0.5_rk
      real(rk), parameter :: alk_per_mn_ox  = -2.0_rk
      real(rk), parameter :: eps_phi        = 1.0e-7_rk

      _LOOP_BEGIN_

         _GET_(self%id_mn2, mn2)
         _GET_(self%id_o2,  o2)

         !--------------------------------------------------------------------
         !  Conversion factors for cross-phase reactions (corrected by porosity)
         !--------------------------------------------------------------------
         if (_AVAILABLE_(self%id_porosity)) then
            _GET_(self%id_porosity, phi)
         else
            phi = 1.0_rk
         end if

         phi   = max(min(phi, 1.0_rk), 0.0_rk)
         phi_s = max(1.0_rk - phi, 0.0_rk)

         if (phi > eps_phi .and. phi < 1.0_rk - eps_phi) then
            d2p = phi / phi_s
         else
            d2p = 1.0_rk
         end if

         ! ------------------------------------------------------------------
         ! Mn(II) oxidation
         !
         ! Mn2+ + 0.5 O2 + 2HCO3- -> MnO2(s) + 2 CO2 + H2O
         !
         ! Consumes dissolved Mn2 and O2, produces particulate MnO2.
         ! Alkalinity decreases by 2 equivalents per mol Mn oxidised 
         ! because 2 equivalents of HCO3⁻ are removed
         ! No explicit DIC changes.
         ! ------------------------------------------------------------------

         ! O2 activation using a smooth step function.
         ! fO2 = 0 when O2 = 0, smoothly increases for trace O2,
         ! and reaches 1 when O2 >= o2_mn_ox_sat.
         xO2 = max(0.0_rk, min(1.0_rk, max(o2, 0.0_rk) / self%o2_mn_ox_sat))
         fO2 = xO2*xO2*xO2 * (xO2 * (6.0_rk*xO2 - 15.0_rk) + 10.0_rk)

         ! Mn(II) oxidation is represented as first-order in Mn2, with a capped
         ! maximum rate. Oxygen acts as an activation factor via a smoothstep
         ! function: oxidation is zero at O2 = 0, increases smoothly at low O2,
         ! and reaches full activity above a small O2 threshold.
         !
         ! This replaces the standard mass-action formulation (rate ∝ Mn2 × O2),
         ! which can lead to very fast, numerically stiff reactions under oxic
         ! conditions. The present formulation preserves O2 control of the process
         ! while avoiding excessive rates and improving numerical stability.
         mn_ox = self%k_mn_ox * fO2 * max(mn2, 0.0_rk)

         o2_cons_mn_ox = o2_per_mn_ox * mn_ox     ! Consumed O2
         alk_change = alk_per_mn_ox * mn_ox       ! Alkalinity changes

         ! Mn2+ is in the dissolved phase and MnO2 in the solid phase
         ! Therefore multiplying by the conversion factor
         mno2_prod = d2p * mn_ox

         _ADD_SOURCE_(self%id_mn2,  -mn_ox)
         _ADD_SOURCE_(self%id_mno2,  mno2_prod)

         _ADD_SOURCE_(self%id_o2,  -o2_cons_mn_ox)
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_change)

         ! -- Diagnostics (mmol m-3 d-1)
         _SET_DIAGNOSTIC_(self%id_mn_ox, mn_ox * secs_per_day)

      _LOOP_END_

   end subroutine do

end module manganese