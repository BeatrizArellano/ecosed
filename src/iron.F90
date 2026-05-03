#include "fabm_driver.h"

! ---------------------------------------------------------------------
! Iron redox cycle
!
! Minimal representation:
!   - fe2   : dissolved reduced iron, Fe(II)
!   - fe3ox : particulate reactive Fe(III) oxide/hydroxide pool
!
! Oxidation written as:
!   Fe2+ + 0.25 O2 + 2 HCO3- + 0.5 H2O -> Fe(OH)3(s) + 2 CO2
!
! Model representation:
!   - Fe2 is consumed and Fe3ox is produced (1:1 in Fe units)
!   - O2 is consumed at 0.25 mol O2 per mol Fe2 oxidised
!   - Alkalinity decreases by 2 eq per mol Fe2 oxidised
!   - Fe3ox is particulate and sinks
!
! Fe(III) reduction by organic matter is handled in the OM degradation
! module, because it is an organic matter remineralisation pathway.
! ---------------------------------------------------------------------
module iron

   use fabm_types
   use molecular_diff, only: DIFF_ION_LINEAR, m0_Fe2, m1_Fe2
   implicit none
   private

   type, extends(type_base_model), public :: type_iron

      ! --- State variables
      type(type_state_variable_id) :: id_fe2
      type(type_state_variable_id) :: id_fe3ox

      ! --- Optional couplings
      type(type_state_variable_id) :: id_o2
      type(type_state_variable_id) :: id_alk

      ! --- Dependencies
      type(type_dependency_id) :: id_porosity

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_fe_ox

      ! --- Parameters
      real(rk) :: k_fe_ox        ! Maximum Fe2 oxidation rate (s-1 internally)
      real(rk) :: o2_fe_ox_sat   ! O2 concentration where Fe oxidation is fully active

   contains
      procedure :: initialize
      procedure :: do
   end type type_iron

contains

   subroutine initialize(self, configunit)
      class(type_iron), intent(inout), target :: self
      integer,          intent(in)            :: configunit

      real(rk), parameter :: d_per_s       = 1.0_rk / 86400.0_rk
      real(rk), parameter :: m_d_per_m_s   = 1.0_rk / 86400.0_rk
      real(rk), parameter :: eps           = 1.0e-12_rk

      real(rk) :: w_fe3ox

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%k_fe_ox, 'k_fe_ox', 'd-1', 'Maximum first-order Fe(II) oxidation rate', default=5.0_rk, scale_factor=d_per_s, minimum=0.0_rk)
      call self%get_parameter(self%o2_fe_ox_sat, 'o2_fe_ox_sat', 'mmol m-3', 'O2 concentration above which Fe(II) oxidation is fully active', default=0.1_rk, minimum=eps)
      call self%get_parameter(w_fe3ox, 'w_fe3ox', 'm d-1', 'Sinking velocity of particulate Fe(III) oxides', default=-1.0_rk, scale_factor=m_d_per_m_s)

      ! ---------------- State variables ----------------
      call self%register_state_variable(self%id_fe2, 'fe2', 'mmol m-3', 'Dissolved Fe(II)', initial_value=0.001_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_fe2, 'is_solute', .true.)
      call self%set_variable_property(self%id_fe2, 'diff_method', DIFF_ION_LINEAR)
      call self%set_variable_property(self%id_fe2, 'm0', m0_Fe2)
      call self%set_variable_property(self%id_fe2, 'm1', m1_Fe2)

      call self%register_state_variable(self%id_fe3ox, 'fe3ox', 'mmol m-3', 'Particulate reactive Fe(III) oxides/hydroxides', &
                                        initial_value=0.001_rk, minimum=0.0_rk, vertical_movement=w_fe3ox)
      call self%set_variable_property(self%id_fe3ox, 'is_solute', .false.)

      ! ---------------- Dependencies ----------------
      call self%register_dependency(self%id_porosity, type_interior_standard_variable(name='porosity', units='1'), required=.false.)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_o2, 'o2', 'mmol m-3', 'Dissolved oxygen', required=.true.)
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.false.)

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_fe_ox, 'fe2_ox', 'mmol m-3 d-1', 'Fe(II) oxidation rate')

   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)
      class(type_iron), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: fe2, o2
      real(rk) :: fO2, xO2
      real(rk) :: fe_ox
      real(rk) :: o2_cons_fe_ox
      real(rk) :: alk_change
      real(rk) :: phi, phi_s, d2p
      real(rk) :: fe3ox_prod

      real(rk), parameter :: secs_per_day   = 86400.0_rk
      real(rk), parameter :: o2_per_fe_ox   = 0.25_rk
      real(rk), parameter :: alk_per_fe_ox  = -2.0_rk
      real(rk), parameter :: eps_phi        = 1.0e-7_rk

      _LOOP_BEGIN_

         _GET_(self%id_fe2, fe2)
         _GET_(self%id_o2, o2)

         !--------------------------------------------------------------------
         ! Conversion factors for cross-phase reactions.
         ! fe2 is dissolved; fe3ox is particulate.
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
         ! Fe(II) oxidation
         !
         ! Fe2+ + 0.25 O2 + 2HCO3- + 0.5H2O -> Fe(OH)3(s) + 2CO2
         !
         ! Consumes dissolved Fe2 and O2, produces particulate Fe3ox.
         ! Alkalinity decreases by 2 equivalents per mol Fe oxidised.
         ! No explicit DIC change: HCO3- is converted to CO2.
         ! ------------------------------------------------------------------       

         ! Fe(II) oxidation is represented as first-order in Fe2, with a capped
         ! maximum rate. Oxygen acts as an activation factor via a smoothstep
         ! function: oxidation is zero at O2 = 0, increases smoothly at low O2,
         ! and reaches full activity above a small O2 threshold.
         !
         ! This replaces the standard mass-action formulation (rate ∝ Fe2 × O2),
         ! which can lead to very fast, numerically stiff reactions under oxic
         ! conditions. The present formulation preserves O2 control of the process
         ! while avoiding excessive rates and improving numerical stability.         
         xO2 = max(0.0_rk, min(1.0_rk, max(o2, 0.0_rk) / self%o2_fe_ox_sat))
         fO2 = xO2*xO2*xO2 * (xO2 * (6.0_rk*xO2 - 15.0_rk) + 10.0_rk)

         fe_ox = self%k_fe_ox * fO2 * max(fe2, 0.0_rk)

         o2_cons_fe_ox = o2_per_fe_ox * fe_ox
         alk_change    = alk_per_fe_ox * fe_ox

         ! Fe2+ is in the dissolved phase and Fe(OH)3/Fe2O3  in the solid phase
         ! Therefore multiplying by the conversion factor
         fe3ox_prod = d2p * fe_ox

         _ADD_SOURCE_(self%id_fe2,   -fe_ox)
         _ADD_SOURCE_(self%id_fe3ox,  fe3ox_prod)

         _ADD_SOURCE_(self%id_o2, -o2_cons_fe_ox)
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_change)

         _SET_DIAGNOSTIC_(self%id_fe_ox, fe_ox * secs_per_day)

      _LOOP_END_

   end subroutine do

end module iron