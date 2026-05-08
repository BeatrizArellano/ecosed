#include "fabm_driver.h"

! ---------------------------------------------------------------------
! Pyrite diagenesis
!
! Minimal representation:
!   - pyrite : particulate/refractory pyrite, FeS2
!
! Lumped pyrite formation:
!   Fe2+ + 2 HS- -> FeS2(s) + H2
!
! Model representation:
!   - Fe2 and dissolved sulfide are consumed in a 1:2 ratio
!   - Pyrite is produced as a particulate solid
!   - Alkalinity decreases by 2 eq per mol pyrite formed
!   - No explicit DIC change
!
! Pyrite oxidation by oxygen:
!   FeS2(s) + 15/4 O2 + 7/2 H2O -> Fe(OH)3(s) + 2 SO4-- + 4 H+
!
! Model representation:
!   - Pyrite is consumed
!   - O2 is consumed at 3.75 mol O2 per mol pyrite oxidised
!   - SO4 is produced at 2 mol SO4 per mol pyrite oxidised
!   - Fe(III) oxide/hydroxide is produced at 1 mol Fe per mol pyrite oxidised
!   - Alkalinity decreases by 4 eq per mol pyrite oxidised
!
! Notes on kinetics:
!   - Pyrite formation is represented with mass-action kinetics:
!       r_form = k_pyrite_form * [Fe2+] * [sulfide]^2
!     with k_pyrite_form in m6 mmol-2 d-1 before FABM scale conversion.
!
!   - Pyrite oxidation is represented with mass-action kinetics:
!       r_ox = k_pyrite_ox_o2 * [pyrite] * [O2]
!     with k_pyrite_ox_o2 in m3 mmol-1 d-1 before FABM scale conversion.
! ---------------------------------------------------------------------
module pyrite

   use fabm_types
   implicit none
   private

   type, extends(type_base_model), public :: type_pyrite

      ! --- State variables
      type(type_state_variable_id) :: id_pyrite

      ! --- Couplings
      type(type_state_variable_id) :: id_fe2
      type(type_state_variable_id) :: id_fe3ox
      type(type_state_variable_id) :: id_sulfide
      type(type_state_variable_id) :: id_o2
      type(type_state_variable_id) :: id_so4

      ! --- Optional couplings
      type(type_state_variable_id) :: id_alk

      ! --- Dependencies
      type(type_dependency_id) :: id_porosity

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_pyrite_form
      type(type_diagnostic_variable_id) :: id_pyrite_ox_o2
      type(type_diagnostic_variable_id) :: id_alk_chn_pyrite_form
      type(type_diagnostic_variable_id) :: id_alk_chn_pyrite_ox

      ! --- Parameters
      real(rk) :: k_pyrite_form       ! Lumped pyrite formation, m6 mmol-2 s-1 effectively
      real(rk) :: k_pyrite_ox_o2      ! Pyrite oxidation by O2, m3 mmol-1 s-1 effectively

   contains
      procedure :: initialize
      procedure :: do
   end type type_pyrite

contains

   subroutine initialize(self, configunit)
      class(type_pyrite), intent(inout), target :: self
      integer,             intent(in)            :: configunit

      real(rk), parameter :: d_per_s     = 1.0_rk / 86400.0_rk
      real(rk), parameter :: m_d_per_m_s = 1.0_rk / 86400.0_rk

      real(rk) :: w_pyrite

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%k_pyrite_form, 'k_pyrite_form', 'm3 mmol-1 d-1', 'Mass-action rate constant for lumped pyrite formation from Fe2 and sulfide', &
                              default=1.0e-3_rk, scale_factor=d_per_s, minimum=0.0_rk)

      call self%get_parameter(self%k_pyrite_ox_o2, 'k_pyrite_ox_o2', 'm3 mmol-1 d-1', 'Apparent mass-action pyrite oxidation rate constant by oxygen', &
                              default=1.0e-3_rk, scale_factor=d_per_s, minimum=0.0_rk)

      call self%get_parameter(w_pyrite, 'w_pyrite', 'm d-1', 'Sinking velocity of particulate pyrite', default=-100.0_rk, scale_factor=m_d_per_m_s)

      ! ---------------- State variables ----------------
      call self%register_state_variable(self%id_pyrite, 'pyrite', 'mmol m-3', 'Particulate pyrite FeS2', initial_value=0.0_rk, minimum=0.0_rk, vertical_movement=w_pyrite)
      call self%set_variable_property(self%id_pyrite, 'is_solute', .false.)

      ! ---------------- Dependencies ----------------
      call self%register_dependency(self%id_porosity, type_interior_standard_variable(name='porosity', units='1'), required=.false.)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_fe2,     'fe2',     'mmol m-3',    'Dissolved Fe(II)', required=.true.)
      call self%register_state_dependency(self%id_fe3ox,   'fe3ox',   'mmol m-3',    'Particulate reactive Fe(III) oxide/hydroxide', required=.true.)
      call self%register_state_dependency(self%id_sulfide, 'sulfide', 'mmol m-3',    'Dissolved sulfide', required=.true.)
      call self%register_state_dependency(self%id_o2,      'o2',      'mmol m-3',    'Dissolved oxygen', required=.true.)
      call self%register_state_dependency(self%id_so4,     'so4',     'mmol m-3',    'Sulfate', required=.true.)
      call self%register_state_dependency(self%id_alk,     'alk',     'mmol eq m-3', 'Total alkalinity', required=.false.)

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_pyrite_form, 'pyrite_form', 'mmol m-3 d-1', 'Pyrite formation rate from Fe2 and sulfide')
      call self%register_diagnostic_variable(self%id_pyrite_ox_o2, 'pyrite_ox_o2', 'mmol m-3 d-1', 'Pyrite oxidation rate by oxygen')
      call self%register_diagnostic_variable(self%id_alk_chn_pyrite_form, 'alk_pyrite_form', 'mmol eq m-3 d-1', 'Alkalinity change by lumped pyrite formation')
      call self%register_diagnostic_variable(self%id_alk_chn_pyrite_ox, 'alk_pyrite_ox', 'mmol eq m-3 d-1', 'Alkalinity change by pyrite oxidation')

   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)
      class(type_pyrite), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: fe2, sulfide, o2, pyrite
      real(rk) :: phi, phi_s, d2p, p2d
      real(rk) :: pyrite_form
      real(rk) :: pyrite_prod_part
      real(rk) :: pyrite_ox_o2
      real(rk) :: o2_from_pyrite_ox
      real(rk) :: so4_from_pyrite_ox
      real(rk) :: fe3ox_from_pyrite_ox
      real(rk) :: alk_change_pyrite_form
      real(rk) :: alk_change_pyrite_ox

      logical  :: is_water

      real(rk), parameter :: secs_per_day         = 86400.0_rk
      real(rk), parameter :: o2_per_pyrite_ox     = 3.75_rk
      real(rk), parameter :: so4_per_pyrite_ox    = 2.0_rk
      real(rk), parameter :: fe3ox_per_pyrite_ox  = 1.0_rk
      real(rk), parameter :: alk_per_pyrite_form  = -2.0_rk
      real(rk), parameter :: alk_per_pyrite_ox    = -4.0_rk
      real(rk), parameter :: eps_phi              = 1.0e-7_rk

      _LOOP_BEGIN_

         _GET_(self%id_fe2,     fe2)
         _GET_(self%id_sulfide, sulfide)
         _GET_(self%id_o2,      o2)
         _GET_(self%id_pyrite,  pyrite)

         !--------------------------------------------------------------------
         ! Phase conversion factors.
         !
         ! Dissolved variables are per porewater volume.
         ! Particulate variables are per solid volume.
         !
         ! d2p converts dissolved-phase rates to particulate source terms.
         ! p2d converts particulate-phase rates to dissolved source terms.
         !--------------------------------------------------------------------
         if (_AVAILABLE_(self%id_porosity)) then
            _GET_(self%id_porosity, phi)
         else
            phi = 1.0_rk
         end if

         phi   = max(min(phi, 1.0_rk), 0.0_rk)
         phi_s = max(1.0_rk - phi, 0.0_rk)

         if (phi > eps_phi .and. phi < 1.0_rk - eps_phi) then
            is_water = .false.
            d2p = phi / phi_s
            p2d = phi_s / phi
         else
            is_water = .true.
            d2p = 1.0_rk
            p2d = 1.0_rk
         end if

         !--------------------------------------------------------------------
         ! Lumped pyrite precipitation
         !
         ! Fe2+ + 2HS- -> FeS2(s) + H2
         !
         ! This reaction represents the net effect of intermediate sulfur
         ! and iron transformations involved in pyrite formation, such as:
         ! Important intermediate steps include:
         !
         !   Fe2+ + HS- -> FeS(s) + H+
         !
         !   2FeOOH + 3H2S -> 2FeS + S0 + 4 H2O
         !
         ! In this approach, intermediate sulfur species are not represented 
         ! explicitly because their cycling is poorly constrained and would 
         ! require additional chemistry (e.g. elemental sulfur and polysulfides)
         ! that adds complexity with limited benefit for the present alkalinity-
         ! centered implementation. Many intermediate alkalinity changes are
         ! also expected to largely offset each other over the full reaction
         ! pathway when FeS and subproducts are re-oxidised. 
         !
         ! Fe2 and sulfide are dissolved, so the base rate is dissolved-phase.
         ! Pyrite production is converted to particulate phase with d2p.
         !--------------------------------------------------------------------
         pyrite_form = self%k_pyrite_form * max(fe2, 0.0_rk) * max(sulfide, 0.0_rk)

         pyrite_prod_part        = d2p * pyrite_form
         alk_change_pyrite_form  = alk_per_pyrite_form * pyrite_form

         !--------------------------------------------------------------------
         ! Pyrite oxidation by oxygen
         !
         ! FeS2(s) + 15/4 O2 + 7/2 H2O -> Fe(OH)3(s) + 2 SO4-- + 4 H+
         !
         ! Pyrite is particulate and O2 is dissolved. The primary rate is
         ! computed in particulate units, then dissolved sinks/products are
         ! converted to porewater units with p2d.
         !--------------------------------------------------------------------
         pyrite_ox_o2 = self%k_pyrite_ox_o2 * max(pyrite, 0.0_rk) * max(o2, 0.0_rk)

         o2_from_pyrite_ox     = -p2d * o2_per_pyrite_ox    * pyrite_ox_o2
         so4_from_pyrite_ox    =  p2d * so4_per_pyrite_ox   * pyrite_ox_o2
         fe3ox_from_pyrite_ox  =         fe3ox_per_pyrite_ox * pyrite_ox_o2
         alk_change_pyrite_ox  =  p2d * alk_per_pyrite_ox   * pyrite_ox_o2

         !--------------------------------------------------------------------
         ! Sources/sinks
         !--------------------------------------------------------------------

         ! Pyrite precipitation
         _ADD_SOURCE_(self%id_fe2,     -pyrite_form)
         _ADD_SOURCE_(self%id_sulfide, -2.0_rk * pyrite_form)
         _ADD_SOURCE_(self%id_pyrite,   pyrite_prod_part)
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_change_pyrite_form)

         ! Pyrite oxidation by O2
         _ADD_SOURCE_(self%id_pyrite, -pyrite_ox_o2)
         _ADD_SOURCE_(self%id_o2,      o2_from_pyrite_ox)
         _ADD_SOURCE_(self%id_so4,     so4_from_pyrite_ox)
         _ADD_SOURCE_(self%id_fe3ox,   fe3ox_from_pyrite_ox)
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_change_pyrite_ox)

         ! Diagnostics
         _SET_DIAGNOSTIC_(self%id_pyrite_form,          pyrite_form * secs_per_day)
         _SET_DIAGNOSTIC_(self%id_pyrite_ox_o2,         pyrite_ox_o2 * secs_per_day)
         _SET_DIAGNOSTIC_(self%id_alk_chn_pyrite_form,  alk_change_pyrite_form * secs_per_day)
         _SET_DIAGNOSTIC_(self%id_alk_chn_pyrite_ox,    alk_change_pyrite_ox * secs_per_day)

      _LOOP_END_

   end subroutine do

end module pyrite
