#include "fabm_driver.h"

! ---------------------------------------------------------------------
! Calcium carbonate mineral dynamics
!
! Minimal representation:
!   - calcite : particulate CaCO3
!
! Processes:
!   1. Calcite dissolution
!        CaCO3(s) -> Ca2+ + CO3--
!
!        calcite decreases by 1
!        DIC     increases by 1
!        ALK     increases by 2
!
!   2. Abiotic calcite precipitation, sediments only
!        Ca2+ + CO3-- -> CaCO3(s)
!
!        calcite increases by 1
!        DIC     decreases by 1
!        ALK     decreases by 2
!
! Notes:
!   - Calcium is not represented explicitly.
!   - Calcite saturation state is supplied by carbonate chemistry.
! ---------------------------------------------------------------------
module calcium_carbonate

   use fabm_types
   implicit none
   private

   type, extends(type_base_model), public :: type_calcium_carbonate

      ! --- State variables
      type(type_state_variable_id) :: id_calcite

      ! --- Couplings to carbonate chemistry state variables
      type(type_state_variable_id) :: id_dic
      type(type_state_variable_id) :: id_alk

      ! --- Dependencies from carbonate chemistry diagnostics
      type(type_dependency_id) :: id_omega_ca

      ! --- Optional dependencies
      type(type_dependency_id) :: id_porosity

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_calcite_diss
      type(type_diagnostic_variable_id) :: id_calcite_precip
      type(type_diagnostic_variable_id) :: id_alk_calcite_diss
      type(type_diagnostic_variable_id) :: id_alk_calcite_precip

      ! --- Parameters
      real(rk) :: k_diss_near     ! Near-equilibrium dissolution rate constant, s-1 internally
      real(rk) :: k_diss_far      ! Strongly undersaturated dissolution rate constant, s-1 internally
      real(rk) :: k_prec          ! Abiotic precipitation rate constant, internally converted to s-1
      logical  :: save_process_rates
      logical  :: save_alkalinity_changes

   contains
      procedure :: initialize
      procedure :: do
   end type type_calcium_carbonate

contains

   subroutine initialize(self, configunit)
      class(type_calcium_carbonate), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      real(rk), parameter :: d_per_s     = 1.0_rk / 86400.0_rk
      real(rk), parameter :: yr_per_s    = 1.0_rk / (86400.0_rk * 365.0_rk)
      real(rk), parameter :: m_d_per_m_s = 1.0_rk / 86400.0_rk

      real(rk) :: w_calcite

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%k_diss_near, 'k_diss_near', 'yr-1', 'Near-equilibrium calcite dissolution rate constant', &
                              default=0.00632_rk, scale_factor=yr_per_s, minimum=0.0_rk)

      call self%get_parameter(self%k_diss_far, 'k_diss_far', 'yr-1', 'Strongly undersaturated calcite dissolution rate constant', &
                              default=20.0_rk, scale_factor=yr_per_s, minimum=0.0_rk)

      call self%get_parameter(self%k_prec, 'k_prec', 'mmol m-3 yr-1', 'Abiotic calcite precipitation coefficient', &
                              default=0.0_rk, scale_factor=yr_per_s, minimum=0.0_rk)

      call self%get_parameter(w_calcite, 'w_calcite', 'm d-1', 'Sinking velocity of calcite', &
                              default=-150.0_rk, scale_factor=m_d_per_m_s)

      call self%get_parameter(self%save_process_rates, 'process_rates', '', 'Save process-rate diagnostics', default=.false.)
      call self%get_parameter(self%save_alkalinity_changes, 'alkalinity_changes', '', 'Save alkalinity-change diagnostics', default=.false.)

      ! ---------------- State variables ----------------
      call self%register_state_variable(self%id_calcite, 'calcite', 'mmol m-3', 'Calcite (CaCO3)', &
                                        initial_value=0.0_rk, minimum=0.0_rk, vertical_movement=w_calcite)
      call self%set_variable_property(self%id_calcite, 'is_solute', .false.)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol m-3', 'Dissolved inorganic carbon', required=.true.)
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.true.)

      ! ---------------- Dependencies ----------------
      call self%register_dependency(self%id_omega_ca, 'omega_ca', '1', 'Calcite saturation state', required=.true.)
      call self%register_dependency(self%id_porosity, type_interior_standard_variable(name='porosity', units='1'), required=.false.)

      ! ---------------- Diagnostics ----------------
      if (self%save_process_rates) then
         call self%register_diagnostic_variable(self%id_calcite_diss, 'CALCITE_DISS', 'mmol m-3 d-1', 'Calcite dissolution rate')
         call self%register_diagnostic_variable(self%id_calcite_precip, 'CALCITE_PRECIP', 'mmol m-3 d-1', 'Abiotic calcite precipitation rate')
      end if

      if (self%save_alkalinity_changes) then
         call self%register_diagnostic_variable(self%id_alk_calcite_diss, 'ALK_CALCITE_DISS', 'mmol eq m-3 d-1', &
                                                'Alkalinity change due to calcite dissolution')
         call self%register_diagnostic_variable(self%id_alk_calcite_precip, 'ALK_CALCITE_PRECIP', 'mmol eq m-3 d-1', &
                                                'Alkalinity change due to calcite precipitation')
      end if

   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)
      class(type_calcium_carbonate), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: calcite
      real(rk) :: omega_ca
      real(rk) :: r_diss, r_prec
      real(rk) :: calcite_change
      real(rk) :: dic_change, alk_change

      real(rk) :: phi, phi_s, p2d

      real(rk), parameter :: omega_transition = 0.8275_rk
      real(rk), parameter :: secs_per_day     = 86400.0_rk      
      real(rk), parameter :: eps_phi          = 1.0e-7_rk

      _LOOP_BEGIN_

         _GET_(self%id_calcite, calcite)
         _GET_(self%id_omega_ca, omega_ca)

         !--------------------------------------------------------------------
         ! Conversion factors for cross-phase reactions.
         !
         ! calcite is particulate.
         ! DIC and ALK are dissolved.
         !
         ! In sediments:
         !   particulate tracers are expressed per solid volume
         !   dissolved tracers are expressed per porewater volume.
         !
         ! In the water column:
         !   all tracers are expressed per water volume, so p2d = 1.
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
         else
            p2d = 1.0_rk
         end if

         ! ------------------------------------------------------------------
         ! Calcite dissolution
         !
         !   CaCO3​(s) -> Ca2+ + CO32−​
         ! Generates 1 mol DIC and 2 eq mol Alkalinity per mol Calcite dissolved.
         !
         ! Non-linear kinetics described in Naviaux et al. (2019):
         !   omega_transition < Omega <= 1:
         !      R = [calcite] * k_diss_near * (1 - Omega)^0.11
         !   Omega <= omega_transition:
         !      R = [calcite] * k_diss_far  * (1 - Omega)^4.7
         !
         ! Dissolution rate is computed on the particulate calcite concentration basis.
         ! DIC and ALK sources are converted to dissolved-phase concentration
         ! using p2d.
         ! ------------------------------------------------------------------
         r_diss = 0.0_rk

         if (omega_ca > omega_transition .and. omega_ca <= 1.0_rk) then
            r_diss = calcite * self%k_diss_near * (1.0_rk - omega_ca)**0.11_rk

         else if (omega_ca <= omega_transition) then
            r_diss = calcite * self%k_diss_far * (1.0_rk - omega_ca)**4.7_rk
         end if

         ! ------------------------------------------------------------------
         ! Abiotic calcite precipitation (only in sediments)
         ! R_prec is computed on the particulate phase.
         !
         !   Ca2+ + CO32− ​-> CaCO3​(s)
         !
         ! Consumes 1 mol DIC and 2 eq mol Alkalinity per mol CaCO3 precipitated.
         !
         ! Reaction kinetics adapted from Zuddas & Mucci (1998) 
         !    R = k_prec * (Omega - 1)^1.76       for Omega > 1
         !
         ! 1.76 corresponds to seawater ionic strength (I ≈ 0.7)
         ! ------------------------------------------------------------------
         r_prec = 0.0_rk

         if (self%k_prec > 0.0_rk .and. omega_ca > 1.0_rk .and. phi < 1.0_rk - eps_phi) then
               r_prec = self%k_prec * (omega_ca - 1.0_rk)**1.76_rk  
         end if

         ! ------------------------------------------------------------------
         ! Source terms
         !
         ! Dissolution:
         !   calcite -1, DIC +1, ALK +2
         !
         ! Precipitation:
         !   calcite +1, DIC -1, ALK -2
         ! ------------------------------------------------------------------
         calcite_change = -r_diss + r_prec
         dic_change     =  p2d * (r_diss - r_prec)
         alk_change     = 2.0_rk * p2d * (r_diss - r_prec)

         _ADD_SOURCE_(self%id_calcite, calcite_change)
         _ADD_SOURCE_(self%id_dic,     dic_change)
         _ADD_SOURCE_(self%id_alk,     alk_change)

         if (self%save_process_rates) then
            _SET_DIAGNOSTIC_(self%id_calcite_diss,   r_diss * secs_per_day)
            _SET_DIAGNOSTIC_(self%id_calcite_precip, r_prec * secs_per_day)
         end if

         if (self%save_alkalinity_changes) then
            _SET_DIAGNOSTIC_(self%id_alk_calcite_diss,    2.0_rk * p2d * r_diss * secs_per_day)
            _SET_DIAGNOSTIC_(self%id_alk_calcite_precip, -2.0_rk * p2d * r_prec * secs_per_day)
         end if

      _LOOP_END_

   end subroutine do

end module calcium_carbonate