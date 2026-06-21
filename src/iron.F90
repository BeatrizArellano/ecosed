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
      type(type_state_variable_id) :: id_mn2
      type(type_state_variable_id) :: id_mno2

      ! --- Dependencies
      type(type_dependency_id) :: id_porosity

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_fe_o2_ox
      type(type_diagnostic_variable_id) :: id_fe_mno2_ox
      type(type_diagnostic_variable_id) :: id_alk_fe_o2_ox
      type(type_diagnostic_variable_id) :: id_alk_fe_mno2_ox

      ! --- Parameters
      real(rk) :: k_fe_ox        ! Fe(II) oxidation rate constant (m3 mmol-1 s-1 internally)
      real(rk) :: k_fe_mno2      ! Fe(II) oxidation by MnO2 rate constant (m3 mmol-1 s-1 internally)
      logical  :: save_process_rates
      logical  :: save_alkalinity_changes

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

      real(rk) :: w_fe3ox

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%k_fe_ox, 'k_fe_ox', 'm3 mmol-1 d-1', 'Second-order mass-action rate constant for Fe(II) oxidation by oxygen', &
                              default=0.02738_rk, scale_factor=d_per_s, minimum=0.0_rk)
      call self%get_parameter(self%k_fe_mno2, 'k_fe_mno2', 'm3 mmol-1 d-1', 'Second-order rate constant for Fe(II) oxidation by MnO2', & 
                              default=2.74e-5_rk, scale_factor=d_per_s, minimum=0.0_rk)
      call self%get_parameter(w_fe3ox, 'w_fe3ox', 'm d-1', 'Sinking velocity of particulate Fe(III) oxides', default=-1.0_rk, scale_factor=m_d_per_m_s)

      call self%get_parameter(self%save_process_rates, 'process_rates', '', 'Save process-rate diagnostics', default=.false.)
      call self%get_parameter(self%save_alkalinity_changes, 'alkalinity_changes', '', 'Save alkalinity-change diagnostics', default=.false.)

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
      call self%register_state_dependency(self%id_mn2,  'mn2',  'mmol m-3', 'Dissolved Mn(II)', required=.false.)
      call self%register_state_dependency(self%id_mno2, 'mno2', 'mmol m-3', 'Particulate MnO2', required=.false.)

      ! ---------------- Diagnostics ----------------
      if (self%save_process_rates) then
         call self%register_diagnostic_variable(self%id_fe_o2_ox, 'FE_OX', 'mmol m-3 d-1', 'Fe(II) oxidation rate by oxygen')
         call self%register_diagnostic_variable(self%id_fe_mno2_ox, 'FE_MNO2_OX', 'mmol m-3 d-1', 'Fe(II) oxidation by MnO2')
      end if

      if (self%save_alkalinity_changes) then
         call self%register_diagnostic_variable(self%id_alk_fe_o2_ox, 'ALK_FE_OX', 'mmol eq m-3 d-1', 'Alkalinity change due to Fe(II) oxidation by oxygen')
         call self%register_diagnostic_variable(self%id_alk_fe_mno2_ox, 'ALK_FE_MNO2_OX', 'mmol eq m-3 d-1', 'Alkalinity change due to Fe(II) oxidation by MnO2')
      end if

   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)
      class(type_iron), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: fe2, o2
      real(rk) :: mn2, mno2
      real(rk) :: fe_ox, fe_mno2
      real(rk) :: o2_cons_fe_ox
      real(rk) :: alk_change
      real(rk) :: alk_fe_o2_change, alk_fe_mno2_change
      real(rk) :: phi, phi_s, d2p, p2d
      real(rk) :: fe3ox_prod, fe2_cons

      real(rk), parameter :: secs_per_day   = 86400.0_rk
      real(rk), parameter :: o2_per_fe_ox   = 0.25_rk
      real(rk), parameter :: alk_per_fe_ox  = -2.0_rk
      real(rk), parameter :: eps_phi        = 1.0e-7_rk

      _LOOP_BEGIN_

         _GET_(self%id_fe2, fe2)
         _GET_(self%id_o2, o2)

         alk_fe_mno2_change = 0.0_rk
         alk_fe_o2_change   = 0.0_rk
         alk_change         = 0.0_rk

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
            p2d = phi_s / phi
            d2p = phi / phi_s
         else
            p2d = 1.0_rk  
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

         ! Fe(II) oxidation is represented using second-order mass-action
         ! kinetics:
         !
         !   rate = k_fe_ox * [Fe2+] * [O2]
         !
         ! Fe2 and O2 are dissolved-phase tracers, so the rate is computed per
         ! porewater volume in sediments and per water volume in the water column.
         ! Fe3ox is particulate, so its production is converted from dissolved-
         ! phase to solid-phase concentration using d2p.
         fe_ox = self%k_fe_ox * max(fe2, 0.0_rk) * max(o2, 0.0_rk)

         o2_cons_fe_ox    = o2_per_fe_ox * fe_ox
         alk_fe_o2_change = alk_per_fe_ox * fe_ox

         ! Fe2+ is in the dissolved phase and Fe(OH)3/Fe2O3  in the solid phase
         ! Therefore multiplying by the conversion factor
         fe3ox_prod = d2p * fe_ox
         fe2_cons   = -fe_ox
         
         
         ! ------------------------------------------------------------------
         ! Fe(II) oxidation by MnO2
         !
         ! 2Fe2+ + MnO2 + 2HCO3 + 2H2O -> 2Fe(OH)3 + Mn2+ + 2C02
         !
         ! Consumes dissolved Fe2 and MnO2, produces particulate Fe3ox and Mn2+.
         ! Alkalinity decreases by 2 equivalents per mol reaction.
         ! No explicit DIC change: HCO3- is converted to CO2.
         ! Fe(II) oxidation is represented using second-order mass-action
         ! kinetics:
         !
         !   rate = k_fe_mno2 * [Fe2+] * [MnO2]
         ! ------------------------------------------------------------------      
         fe_mno2 = 0.0_rk
         if (_AVAILABLE_(self%id_mn2) .and. _AVAILABLE_(self%id_mno2)) then
            _GET_(self%id_mn2,  mn2)
            _GET_(self%id_mno2, mno2)
            ! Reaction rate interpreted on the particulate/solid-volume basis.
            fe_mno2 = self%k_fe_mno2 * max(fe2, 0.0_rk) * max(mno2, 0.0_rk)
            alk_fe_mno2_change = -2.0_rk * p2d * fe_mno2
            fe3ox_prod = fe3ox_prod + 2.0_rk * fe_mno2
            fe2_cons   = fe2_cons   -2.0_rk * p2d  * fe_mno2
         end if

         alk_change = alk_fe_o2_change + alk_fe_mno2_change

         _ADD_SOURCE_(self%id_fe2,    fe2_cons)
         _ADD_SOURCE_(self%id_fe3ox,  fe3ox_prod)

         _ADD_SOURCE_(self%id_o2, -o2_cons_fe_ox)
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_change)

         if (fe_mno2 > 0.0_rk) then
            _ADD_SOURCE_(self%id_mno2,  -fe_mno2)
            _ADD_SOURCE_(self%id_mn2,    p2d *fe_mno2)
         end if

         if (self%save_process_rates) then
            _SET_DIAGNOSTIC_(self%id_fe_o2_ox,   fe_ox   * secs_per_day)
            _SET_DIAGNOSTIC_(self%id_fe_mno2_ox, fe_mno2 * secs_per_day)
         end if

         if (self%save_alkalinity_changes) then
            _SET_DIAGNOSTIC_(self%id_alk_fe_o2_ox,   alk_fe_o2_change   * secs_per_day)
            _SET_DIAGNOSTIC_(self%id_alk_fe_mno2_ox, alk_fe_mno2_change * secs_per_day)
         end if

      _LOOP_END_

   end subroutine do

end module iron