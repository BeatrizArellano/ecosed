#include "fabm_driver.h"

!-------------------------------------------------------------------------------------------------------
! Simple water-column detritus module.
!
! Purpose:
!   Minimal particulate organic matter closure for validating phytoplankton / chlorophyll dynamics.
!
! Currency:
!   det is stored in nitrogen units [mmol N m-3].
!
! Coupling logic:
!   - pelagic_ecosystem provides diagnostic pom_prod_n [mmol N m-3 s-1]
!   - this module stores that production as detritus
!   - detritus remineralises linearly to NH4 and PO4 using fixed N:P stoichiometry
!   - nitrogen.F90 then handles NH4 -> NO3 through nitrification
!-------------------------------------------------------------------------------------------------------
module detritus

   use fabm_types
   implicit none
   private

   type, extends(type_base_model), public :: type_detritus
      ! --- State variables
      type(type_state_variable_id) :: id_det

      ! --- Couplings to dissolved nutrient module
      type(type_state_variable_id) :: id_nh4
      type(type_state_variable_id) :: id_po4
      !--- Optional couplings
      type(type_state_variable_id) :: id_o2
      type(type_state_variable_id) :: id_dic
      type(type_state_variable_id) :: id_alk

      ! --- Dependency provided by pelagic_ecosystem
      type(type_dependency_id) :: id_pom_prod_n

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_total_det
      type(type_diagnostic_variable_id) :: id_remin_det
      type(type_diagnostic_variable_id) :: id_alk_prod_aer

      ! --- Parameters
      real(rk) :: k_remin_det      ! Detritus remineralisation rate [s-1 internally]
      real(rk) :: c_to_n_det       ! C:N ratio [-]
      real(rk) :: n_to_p_det       ! Detrital molar N:P ratio [-]
      real(rk) :: o2_per_c
      real(rk) :: atten_det        ! Detritus-specific PAR attenuation [m2 mmol N-1]
      real(rk) :: w_det            ! Detritus settling velocity [m s-1 internally]
      logical  :: settling         ! Enable detritus settling and seabed loss

      logical  :: save_process_rates
      logical  :: save_alkalinity_changes
      
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
   end type type_detritus

contains

   subroutine initialize(self, configunit)
      class(type_detritus), intent(inout), target :: self
      integer,              intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
      real(rk), parameter :: eps     = 1.0e-12_rk

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%k_remin_det, 'k_remin_det', 'd-1', 'First-order detritus remineralisation rate', default=0.05_rk, &
                              scale_factor=d_per_s, minimum=0.0_rk)

      call self%get_parameter(self%c_to_n_det, 'c_to_n_det', '-', 'Detrital molar C:N ratio', default=6.625_rk, minimum=eps)
      call self%get_parameter(self%n_to_p_det, 'n_to_p_det', '-', 'Detrital molar N:P ratio', default=16.0_rk, minimum=eps)

      call self%get_parameter(self%o2_per_c, 'o2_per_c', 'mol O2 mol C-1', 'Effective O2 consumed per mol organic C remineralised aerobically', default=1.3_rk, minimum=0.0_rk)

      call self%get_parameter(self%w_det, 'w_det', 'm d-1', 'Vertical velocity of detritus (<0 sinking)', default=-5.0_rk, &
                              maximum=0.0_rk, scale_factor=d_per_s)
      call self%get_parameter(self%settling, 'settling', '', 'Enable detritus settling and loss to abstract sediment', default=.false.)

      call self%get_parameter(self%atten_det, 'atten_det', 'm2 mmol-1', 'Specific light extinction of detritus', default=0.03_rk, minimum=0.0_rk)

      call self%get_parameter(self%save_process_rates, 'process_rates', '', 'Save process-rate diagnostics', default=.false.)
      call self%get_parameter(self%save_alkalinity_changes, 'alkalinity_changes', '', 'Save alkalinity-change diagnostics', default=.false.)

      ! ---------------- State variables ----------------
      call self%register_state_variable(self%id_det, 'det', 'mmol N m-3', 'Detritus', initial_value=0.001_rk, minimum=0.0_rk, vertical_movement=self%w_det)
      call self%set_variable_property(self%id_det, 'is_solute', .false.)

      ! Detritus contributes to total nitrogen conservation.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_det)

      ! Detritus contributes to PAR attenuation.
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
                                          self%id_det, scale_factor=self%atten_det)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mmol m-3', 'Dissolved ammonium', required=.true.)
      call self%register_state_dependency(self%id_po4, 'po4', 'mmol m-3', 'Dissolved phosphate', required=.true.)

      call self%register_state_dependency(self%id_o2,  'o2',  'mmol O2 m-3', 'Dissolved oxygen', required=.false.)
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C m-3',  'Dissolved inorganic carbon', required=.false.)
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.false.)

      ! ---------------- Dependencies ----------------
      ! This is provided by pelagic_ecosystem as a diagnostic in N units.
      call self%register_dependency(self%id_pom_prod_n, 'pom_prod_n', 'mmol N m-3 s-1', 'Production of particulate organic matter from pelagic biology',required=.false.)

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_total_det, 'total_det', 'mmol N m-3', 'Total detritus')

      if (self%save_process_rates) then
         call self%register_diagnostic_variable(self%id_remin_det, 'REMIN_DET', 'mmol N m-3 d-1', 'Detritus remineralisation rate')
      end if
      if (self%save_alkalinity_changes) then
         call self%register_diagnostic_variable(self%id_alk_prod_aer, 'ALK_PROD_AER', 'mmol eq m-3 d-1', 'Alkalinity change due to aerobic remineralisation')
      end if

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: det
      real(rk) :: pom_prod_n
      real(rk) :: remin_n, remin_c, remin_p
      real(rk) :: o2_cons_aer
      real(rk) :: alk_prod_aer

      real(rk), parameter :: secs_per_day = 86400.0_rk

      _LOOP_BEGIN_

         _GET_(self%id_det, det)
         pom_prod_n = 0.0_rk
         if (_AVAILABLE_(self%id_pom_prod_n)) then
            _GET_(self%id_pom_prod_n, pom_prod_n)
         end if

         ! Linear remineralisation in N units.
         remin_n = self%k_remin_det * max(det, 0.0_rk)
         remin_c = self%c_to_n_det * remin_n
         remin_p = remin_n / self%n_to_p_det  
         
         ! O2 consumption associated with aerobic remineralisation.
         ! O2 consumption is computed from the amount of organic C remineralised,
         ! using o2_per_c.
         o2_cons_aer = self%o2_per_c * remin_c
         
         ! Alkalinity generation during aerobic remineralisation.
         ! Equation R1 in Middelburg et al. (2020)
         ! (CH2O)(NH3)n/c(H3PO4)p/c + O2 → CO2 + n/c NH3 + p/c H3PO4 + H2O 	Alk change: n/c-p/c per mol C or 1-p/n per mol N   
         alk_prod_aer = remin_n * (1.0_rk - 1.0_rk/self%n_to_p_det)

         ! Detritus receives particulate production from pelagic_ecosystem and loses material by remineralisation.
         _ADD_SOURCE_(self%id_det, pom_prod_n - remin_n)

         ! Remineralised detrital N and P return to dissolved inorganic nutrients.
         _ADD_SOURCE_(self%id_nh4, remin_n)
         _ADD_SOURCE_(self%id_po4, remin_p)

         ! Aerobic remineralisation stoichiometry
         if (_AVAILABLE_(self%id_o2)) _ADD_SOURCE_(self%id_o2, -o2_cons_aer)
         if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic,  remin_c)
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk,  alk_prod_aer)

         ! Diagnostics.
         _SET_DIAGNOSTIC_(self%id_total_det, det)
         if (self%save_process_rates) then
            _SET_DIAGNOSTIC_(self%id_remin_det, remin_n * secs_per_day)
         end if

         if (self%save_alkalinity_changes) then
            _SET_DIAGNOSTIC_(self%id_alk_prod_aer, alk_prod_aer * secs_per_day)
         end if

      _LOOP_END_

   end subroutine do


   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: det
      real(rk) :: fsettle

      if (.not. self%settling) return

      _BOTTOM_LOOP_BEGIN_

         _GET_(self%id_det, det)

         ! Settling velocity is negative for downward movement.
         ! Bottom flux is positive into water, so settling loss is negative.
         fsettle = max(0.0_rk, -self%w_det) * max(det, 0.0_rk)

         _ADD_BOTTOM_FLUX_(self%id_det, -fsettle)

      _BOTTOM_LOOP_END_

   end subroutine do_bottom

end module detritus