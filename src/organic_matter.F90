#include "fabm_driver.h"

! Organic matter module.
! POM pools are stored in N units and partitioned into labile, semi-labile, and refractory fractions.
! Aerobic remineralisation regenerates inorganic nitrogen (currently returned directly to NO3),
! produces DIC according to pool-specific C:N ratios, and consumes O2 according to an effective
! O2:C remineralisation stoichiometry.
module organic_matter

   use fabm_types
   implicit none
   private

   type, extends(type_base_model), public :: type_organic_matter
      ! --- State variables
      type(type_state_variable_id) :: id_pom_l, id_pom_s, id_pom_r

      ! --- Couplings
      type(type_state_variable_id) :: id_no3

      !--- Optional Coupling
      type(type_state_variable_id) :: id_dic
      type(type_state_variable_id) :: id_o2

      ! --- Dependencies
      type(type_dependency_id)     :: id_temp
      type(type_dependency_id)     :: id_porosity
      type(type_dependency_id)     :: id_pom_prod_n

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_rem
      type(type_diagnostic_variable_id) :: id_total_pom
      type(type_diagnostic_variable_id) :: id_o2_cons_rem

      ! --- Parameters
      real(rk) :: frac_lab, frac_semi, frac_ref
      real(rk) :: k_remin_pom_l, k_remin_pom_s, k_remin_pom_r  ! Remineralisation rates
      real(rk) :: o2_remin_thr                                 ! O2 threshold for aerobic remineralisation (mmol m-3)      
      real(rk) :: c_to_n_pom_l
      real(rk) :: c_to_n_pom_s
      real(rk) :: c_to_n_pom_r
      real(rk) :: o2_per_c
      real(rk) :: kc

   contains
      procedure :: initialize
      procedure :: do
   end type type_organic_matter

contains

   subroutine initialize(self, configunit)
      class(type_organic_matter), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
      real(rk) :: w_pom

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%frac_lab,  'frac_lab',  '-', 'Fraction of newly produced labile organic matter', default=0.5_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%frac_semi, 'frac_semi', '-', 'Fraction of newly produced semi-labile organic matter', default=0.4_rk, minimum=0.0_rk, maximum=1.0_rk)

      if (self%frac_lab + self%frac_semi .gt. 1.0_rk) then
         write(*,*) 'yaml_settings: frac_lab + frac_semi must be <= 1'
         stop
      end if
      self%frac_ref = max(0.0_rk, 1.0_rk - self%frac_lab - self%frac_semi)

      call self%get_parameter(self%k_remin_pom_l,'k_remin_pom_l','d-1','Labile OM remineralisation rate (oxic)', default=0.1_rk, scale_factor=d_per_s)
      call self%get_parameter(self%k_remin_pom_s,'k_remin_pom_s','d-1','Semi-labile OM remineralisation rate (oxic)', default=0.01_rk, scale_factor=d_per_s)
      call self%get_parameter(self%k_remin_pom_r,'k_remin_pom_r','d-1','Refractory OM remineralisation rate (oxic)', default=0.001_rk, scale_factor=d_per_s)

      call self%get_parameter(self%o2_remin_thr,'o2_remin_thr','mmol m-3', 'O2 concentration where aerobic remineralisation is half active', default=6.0_rk)

      call self%get_parameter(self%c_to_n_pom_l, 'c_to_n_l', '-', 'Molar C:N ratio of labile POM', default=6.625_rk)
      call self%get_parameter(self%c_to_n_pom_s, 'c_to_n_s', '-', 'Molar C:N ratio of semi-labile POM', default=8.0_rk)
      call self%get_parameter(self%c_to_n_pom_r, 'c_to_n_r', '-', 'Molar C:N ratio of refractory POM', default=10.0_rk)

      call self%get_parameter(self%o2_per_c,  'o2_per_c',  'Effective O2 consumed per mol organic C remineralised aerobically', default=1.3_rk)
      call self%get_parameter(self%kc, 'kc', 'm2 mmol-1', 'specific light extinction of POM', default=0.03_rk)
      call self%get_parameter(w_pom,'w_pom','m d-1','vertical velocity of POM (<0 sinking)', default=-10.0_rk, scale_factor=d_per_s)

      ! ---------------- State variables ----------------
      call self%register_state_variable(self%id_pom_l, 'pom_l', 'mmol m-3', 'POM labile',      0.001_rk, minimum=0.0_rk, vertical_movement=w_pom)
      call self%register_state_variable(self%id_pom_s, 'pom_s', 'mmol m-3', 'POM semi-labile', 0.001_rk, minimum=0.0_rk, vertical_movement=w_pom)
      call self%register_state_variable(self%id_pom_r, 'pom_r', 'mmol m-3', 'POM refractory',  0.001_rk, minimum=0.0_rk, vertical_movement=w_pom)
      call self%set_variable_property(self%id_pom_l, 'is_solute', .false.)
      call self%set_variable_property(self%id_pom_s, 'is_solute', .false.)
      call self%set_variable_property(self%id_pom_r, 'is_solute', .false.)

      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_pom_l)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_pom_s)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_pom_r)

      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_pom_l, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_pom_s, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_pom_r, scale_factor=self%kc)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_no3, 'no3', 'mmol m-3', 'Dissolved nitrate', required=.true.)
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol m-3', 'Total dissolved inorganic carbon', required=.false.)
      call self%register_state_dependency(self%id_o2,  'o2',  'mmol m-3', 'Dissolved oxygen', required=.false.)

      ! ---------------- Dependencies ----------------
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_porosity, type_interior_standard_variable(name='porosity', units='1'), required=.false.)
      call self%register_dependency(self%id_pom_prod_n, 'pom_prod_n', 'mmol m-3 s-1', 'Production of particulate organic Matter from pelagic biology')

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_rem, 'REM_aerobic', 'mmol m-3 d-1', 'Total aerobic remineralisation')
      call self%register_diagnostic_variable(self%id_total_pom, 'total_pom', 'mmol m-3', 'Total particulate organic matter')
      call self%register_diagnostic_variable(self%id_o2_cons_rem, 'O2_CONS_REM', 'mmol m-3 d-1', 'O2 consumption by aerobic remineralisation')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_organic_matter), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: pom_l, pom_s, pom_r
      real(rk) :: temp, pom_prod_n
      real(rk) :: phi, phi_s          ! Porosity
      real(rk) :: p2d, d2p
      real(rk) :: o2, fo2, fT
      real(rk) :: rem_pom_l, rem_pom_s, rem_pom_r, rem_total
      real(rk) :: o2_cons_rem, dic_prod_rem      

      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: o2_width = 0.2_rk            ! Determines the smoothness of the transition from no remineralisation to remineralisation. Larger means smoother, close to 0 is a fast change from 0 to 1
      real(rk), parameter :: eps_phi = 1.0e-6_rk

      _LOOP_BEGIN_

         _GET_(self%id_pom_l, pom_l)
         _GET_(self%id_pom_s, pom_s)
         _GET_(self%id_pom_r, pom_r)
         _GET_(self%id_temp, temp)
         _GET_(self%id_pom_prod_n, pom_prod_n)

         ! Porosity (important if sediments are enabled)
         if (_AVAILABLE_(self%id_porosity)) then
            _GET_(self%id_porosity, phi)
         else 
            phi = 1.0_rk              ! 1 (assuming water)
         end if

         phi   = max(min(phi, 1.0_rk), 0.0_rk)
         phi_s = max(1.0_rk - phi, 0.0_rk)
         ! POM pools are stored on a solid-phase volume basis in sediments.
         ! Dissolved tracers are stored on a porewater-volume basis.
         ! Cross-phase reaction terms are converted using:
         !   p2d = (1-phi)/phi   for particulate -> dissolved
         !   d2p = phi/(1-phi)   for dissolved -> particulate
         if (phi > eps_phi .and. phi < 1.0_rk - eps_phi) then
            p2d = phi_s / phi
            d2p = phi / phi_s
         else
            p2d = 1.0_rk            ! No phase conversion in water
            d2p = 1.0_rk
         end if

         ! Eppley's (1972) exponential scaling of metabolic rates with temperature.
         fT = 1.066_rk ** temp
         ! Smooth Oxic limitation factor for remineralisation
         if (_AVAILABLE_(self%id_o2)) then
            _GET_(self%id_o2, o2)
            fo2 = 0.5_rk * (1.0_rk + tanh((max(o2,0.0_rk) - self%o2_remin_thr) / o2_width))
         else
            fo2 = 1.0_rk
         end if

         ! Aerobic remineralisation (first-order decay process)
         rem_pom_l = self%k_remin_pom_l * fT * fo2 * pom_l
         rem_pom_s = self%k_remin_pom_s * fT * fo2 * pom_s
         rem_pom_r = self%k_remin_pom_r * fT * fo2 * pom_r
         rem_total = rem_pom_l + rem_pom_s + rem_pom_r

         dic_prod_rem = self%c_to_n_pom_l * rem_pom_l   &
                      + self%c_to_n_pom_s * rem_pom_s   &
                      + self%c_to_n_pom_r * rem_pom_r
         
         o2_cons_rem = self%o2_per_c * dic_prod_rem   

         ! POM budgets: production from pelagic biology, loss by remineralisation
         _ADD_SOURCE_(self%id_pom_l, self%frac_lab  * pom_prod_n - rem_pom_l)
         _ADD_SOURCE_(self%id_pom_s, self%frac_semi * pom_prod_n - rem_pom_s)
         _ADD_SOURCE_(self%id_pom_r, self%frac_ref  * pom_prod_n - rem_pom_r)

         ! Remineralisation products 
         _ADD_SOURCE_(self%id_no3, p2d * rem_total)

         if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic, p2d * dic_prod_rem)
         if (_AVAILABLE_(self%id_o2))  _ADD_SOURCE_(self%id_o2, -p2d * o2_cons_rem)

         _SET_DIAGNOSTIC_(self%id_rem, rem_total * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_total_pom, pom_l + pom_s + pom_r)
         _SET_DIAGNOSTIC_(self%id_o2_cons_rem, o2_cons_rem * secs_pr_day)

      _LOOP_END_
   end subroutine do

end module organic_matter