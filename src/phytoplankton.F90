#include "fabm_driver.h"

!-------------------------------------------------------------------------------------------------------
! Phytoplankton functional group.
!
! Phytoplankton growth is controlled by light, temperature, and inorganic nitrogen availability.
! Nitrogen uptake is represented using a smooth two-substrate limitation formulation with uptake from
! both NH4 and NO3, allowing preferential use of NH4 via the half-saturation terms.
!
! The resulting phytoplankton production is partitioned diagnostically into regenerated production
! (NH4-based) and new production (NO3-based).
!
! Phytoplankton mortality contributes to particulate organic matter production.
! Chlorophyll is diagnosed from phytoplankton carbon using a fixed Chl:C ratio.
!
! Multiple instances of this module may be used to represent different phytoplankton functional groups.
!-------------------------------------------------------------------------------------------------------
module phytoplankton

   use fabm_types
   use pelagic_common, only: pom_production_n, total_chlorophyll

   implicit none
   private

   type, extends(type_base_model), public :: type_phytoplankton

      ! --- State variable
      type(type_state_variable_id) :: id_biomass    ! Phytoplankton functional group biomass concentration

      ! --- Couplings
      type(type_state_variable_id) :: id_no3
      type(type_state_variable_id) :: id_nh4

      ! --- Optional couplings
      type(type_state_variable_id) :: id_o2
      type(type_state_variable_id) :: id_dic
      type(type_state_variable_id) :: id_alk

      ! --- Environmental dependencies
      type(type_dependency_id) :: id_par
      type(type_dependency_id) :: id_temp
      type(type_dependency_id) :: id_porosity

      ! --- Internal diagnostics used for aggregate variables
      type(type_diagnostic_variable_id) :: id_pom_prod_n
      type(type_diagnostic_variable_id) :: id_chl

      ! --- Optional diagnostics
      type(type_diagnostic_variable_id) :: id_TPP
      type(type_diagnostic_variable_id) :: id_o2_prod
      type(type_diagnostic_variable_id) :: id_new_prod
      type(type_diagnostic_variable_id) :: id_reg_prod
      type(type_diagnostic_variable_id) :: id_alk_pp

      ! --- Parameters
      real(rk) :: mu_max
      real(rk) :: k_nh4
      real(rk) :: k_no3
      real(rk) :: m_phy

      real(rk) :: chl_per_c
      real(rk) :: alpha_phy

      real(rk) :: c_to_n_phy
      real(rk) :: o2_per_c

      ! --- Optical properties
      real(rk) :: kc

      ! --- Diagnostic switches
      logical :: save_process_rates
      logical :: save_alkalinity_changes

   contains
      procedure :: initialize
      procedure :: do
   end type type_phytoplankton

contains

   subroutine initialize(self, configunit)

      class(type_phytoplankton), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
      real(rk) :: w_phy

      !-------------------------------------------------------------------------------------------------
      ! Parameters
      !-------------------------------------------------------------------------------------------------

      ! Rates are supplied in d-1 in the configuration and converted internally to s-1.
      call self%get_parameter(self%mu_max, 'mu_max', 'd-1', 'Maximum specific phytoplankton growth rate', &
                              default=1.0_rk, scale_factor=d_per_s)

      call self%get_parameter(self%k_nh4, 'k_nh4', 'mmol m-3', 'Half-saturation constant for ammonium uptake', &
                              default=0.2_rk)

      call self%get_parameter(self%k_no3, 'k_no3', 'mmol m-3', 'Half-saturation constant for nitrate uptake', &
                              default=0.5_rk)

      call self%get_parameter(self%alpha_phy, 'alpha_phy', 'm2 W-1 d-1', 'Initial slope of phytoplankton P-I curve', &
                              default=0.03_rk, scale_factor=d_per_s)

      call self%get_parameter(self%chl_per_c, 'chl_per_c', 'mgChl mmolC-1', 'Fixed Chl:C ratio', &
                              default=0.3_rk)

      call self%get_parameter(self%m_phy, 'm_phy', 'd-1', 'Linear phytoplankton mortality', &
                              default=0.05_rk, scale_factor=d_per_s)

      call self%get_parameter(self%c_to_n_phy, 'c_to_n_phy', '-', 'Phytoplankton molar C:N ratio', default=6.625_rk)

      call self%get_parameter(self%o2_per_c, 'o2_per_c', 'mol O2 mol C-1', 'Effective O2 produced per C fixed', &
                              default=1.3_rk)

      call self%get_parameter(self%kc, 'kc', 'm2 mmol-1', 'Specific light extinction of phytoplankton', &
                              default=0.03_rk)

      call self%get_parameter(w_phy, 'w_phy', 'm d-1', 'Vertical velocity of phytoplankton (<0 sinking)', &
                              default=0.0_rk, scale_factor=d_per_s)

      call self%get_parameter(self%save_process_rates, 'process_rates', '', 'Save process-rate diagnostics', &
                              default=.false.)

      call self%get_parameter(self%save_alkalinity_changes, 'alkalinity_changes', '', 'Save alkalinity-change diagnostics', &
                              default=.false.)

      !-------------------------------------------------------------------------------------------------
      ! State variable
      !-------------------------------------------------------------------------------------------------

      call self%register_state_variable(self%id_biomass, 'biomass', 'mmol N m-3', 'Biomass', &
                                        1.0e-12_rk, minimum=0.0_rk, vertical_movement=w_phy)

      call self%set_variable_property(self%id_biomass, 'is_solute', .false.)

      ! Phytoplankton contributes to total nitrogen.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_biomass)

      ! Phytoplankton contributes to light attenuation.
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
                                           self%id_biomass, scale_factor=self%kc)

      !-------------------------------------------------------------------------------------------------
      ! Couplings
      !-------------------------------------------------------------------------------------------------

      call self%register_state_dependency(self%id_no3, 'no3', 'mmol m-3', 'Dissolved nitrate', required=.true.)
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mmol m-3', 'Dissolved ammonium', required=.true.)

      call self%register_state_dependency(self%id_o2, 'o2', 'mmol m-3', 'Dissolved oxygen', required=.false.)
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C m-3', 'Total dissolved inorganic carbon', required=.false.)
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.false.)

      !-------------------------------------------------------------------------------------------------
      ! Environmental dependencies
      !-------------------------------------------------------------------------------------------------

      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      call self%register_dependency(self%id_porosity, type_interior_standard_variable(name='porosity', units='1'), required=.false.)

      !-------------------------------------------------------------------------------------------------
      ! Internal diagnostics contributing to shared aggregate variables
      !-------------------------------------------------------------------------------------------------

      call self%register_diagnostic_variable(self%id_pom_prod_n, 'pom_prod_n', &
                                             'mmol N m-3 s-1', 'Phytoplankton contribution to POM production', &
                                             output=output_none)
      call self%add_to_aggregate_variable(pom_production_n, self%id_pom_prod_n)

      call self%register_diagnostic_variable(self%id_chl, 'chl', 'mg m-3', &
                                             'Chlorophyll-a', output=output_none)
      call self%add_to_aggregate_variable(total_chlorophyll, self%id_chl)

      !-------------------------------------------------------------------------------------------------
      ! Optional process diagnostics
      !-------------------------------------------------------------------------------------------------

      if (self%save_process_rates) then

         call self%register_diagnostic_variable(self%id_TPP, 'TPP', 'mmol N m-3 d-1', &
                                                'Total primary production (N-based)')

         call self%register_diagnostic_variable(self%id_o2_prod, 'O2_PROD', 'mmol m-3 d-1', &
                                                'O2 production from primary production')

         call self%register_diagnostic_variable(self%id_new_prod, 'NEW_PP', 'mmol N m-3 d-1', &
                                                'New production (NO3-supported)')

         call self%register_diagnostic_variable(self%id_reg_prod, 'REG_PP', 'mmol N m-3 d-1', &
                                                'Regenerated production (NH4-supported)')

      end if

      if (self%save_alkalinity_changes) then

         call self%register_diagnostic_variable(self%id_alk_pp, 'ALK_PP', 'mmol eq m-3 d-1', &
                                                'Alkalinity change due to primary production')

      end if

   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)

      class(type_phytoplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: biomass
      real(rk) :: no3, nh4
      real(rk) :: par, temp, phi

      real(rk) :: fN, fT
      real(rk) :: mu_T, alphaI
      real(rk) :: lim_nh4, lim_no3, denom_n
      real(rk) :: tpp, new_prod, reg_prod
      real(rk) :: cphy, jden

      real(rk) :: mort_p
      real(rk) :: phy_source
      real(rk) :: phy_kill

      real(rk) :: dno3, dnh4
      real(rk) :: dic_change
      real(rk) :: alk_change
      real(rk) :: o2_change

      real(rk) :: pom_prod_n
      real(rk) :: chl_diag

      logical :: is_sediment

      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: k_sed_kill_phy = 200.0_rk / secs_pr_day
      real(rk), parameter :: n_to_p = 16.0_rk
      real(rk), parameter :: eps = 1.0e-12_rk

      _LOOP_BEGIN_

         !----------------------------------------------------------------------------------------------
         ! State and environmental conditions
         !----------------------------------------------------------------------------------------------

         _GET_(self%id_biomass, biomass)
         _GET_(self%id_no3, no3)
         _GET_(self%id_nh4, nh4)
         _GET_(self%id_par, par)
         _GET_(self%id_temp, temp)

         ! Porosity to distinguish water-column and sediment layers.
         if (_AVAILABLE_(self%id_porosity)) then
            _GET_(self%id_porosity, phi)
         else
            phi = 1.0_rk
         end if

         is_sediment = (phi < 1.0_rk - eps)

         !----------------------------------------------------------------------------------------------
         ! Water column
         !----------------------------------------------------------------------------------------------

         if (.not. is_sediment) then
            !-------------------------------------------------------------------------------------------
            ! Nitrogen limitation 
            ! Nitrogen limitation is computed with a smooth two-substrate formulation using both NO3 and NH4.
            ! The resulting terms lim_no3 and lim_nh4 represent the contributions of nitrate and ammonium to
            ! total nitrogen limitation, with preferential NH4 use emerging from the lower NH4 half-saturation.
            ! Their sum fN is the total nitrogen limitation factor (0-1).
            !-------------------------------------------------------------------------------------------

            denom_n = self%k_no3 * self%k_nh4 &
                    + self%k_nh4 * max(no3, 0.0_rk) &
                    + self%k_no3 * max(nh4, 0.0_rk)

            if (denom_n > eps) then
               lim_no3 = max(no3, 0.0_rk) * self%k_nh4 / denom_n
               lim_nh4 = max(nh4, 0.0_rk) * self%k_no3 / denom_n
            else
               lim_no3 = 0.0_rk
               lim_nh4 = 0.0_rk
            end if

            fN = lim_no3 + lim_nh4

            !-------------------------------------------------------------------------------------------
            ! Temperature dependence
            !-------------------------------------------------------------------------------------------
            ! Eppley's (1972) exponential scaling of metabolic rates with temperature.
            ! Sets the temperature-dependent maximum growth rate.
            fT = 1.066_rk ** temp
            mu_T = self%mu_max * fT         ! Temperature-adjusted maximum growth rate

            !-------------------------------------------------------------------------------------------
            ! Light limitation
            !-------------------------------------------------------------------------------------------
            ! Light-limited potential phytoplankton growth rate.
            alphaI = self%alpha_phy * max(par, 0.0_rk)

            !-------------------------------------------------------------------------------------------
            ! Primary production
            !-------------------------------------------------------------------------------------------
            ! Primary production follows a smooth saturating light-response formulation,
            ! transitioning between light-limited (alphaI) and temperature-limited (mu_T) growth.
            ! The resulting N-based production is further reduced by total nitrogen limitation fN.
            jden = sqrt(mu_T * mu_T + alphaI * alphaI + eps)
            tpp = (mu_T * alphaI / jden) * fN * biomass

            ! New and regenerated production
            ! Partition total phytoplankton production into NO3-supported (new) and NH4-supported
            ! (regenerated) production using the relative NO3 and NH4 limitation terms.
            new_prod = tpp * lim_no3 / (fN + eps)
            reg_prod = tpp * lim_nh4 / (fN + eps)

            !-------------------------------------------------------------------------------------------
            ! Biomass and chlorophyll
            !-------------------------------------------------------------------------------------------
            ! Convert N-based phytoplankton biomass to carbon.
            cphy = self%c_to_n_phy * max(biomass, 0.0_rk)
            ! Compute chlorophyll from phytoplankton carbon biomass
            chl_diag = self%chl_per_c * cphy

            ! Linear phytoplankton mortality: constant fractional loss to POM
            mort_p = self%m_phy * biomass

            ! Chlorophyll decreases in proportion to phytoplankton biomass loss
            ! Grazing is not included here.
            ! It will be applied to id_phy by the zooplankton module.
            phy_source = tpp - mort_p

            ! Mortality contributes particulate organic matter.
            pom_prod_n = mort_p

            !-------------------------------------------------------------------------------------------
            ! Dissolved nutrient and carbonate-system effects
            !-------------------------------------------------------------------------------------------
            ! Changes in inorganic nitrogen pools:
            ! - regenerated production consumes NH4
            ! - new production consumes NO3
            dno3 = -new_prod
            dnh4 = -reg_prod

            ! Changes in DIC
            ! DIC decreases through phytoplankton carbon fixation and increases through zooplankton respiration.
            dic_change = -self%c_to_n_phy * tpp

            ! Changes in Alkalinity
            ! New production removes NO3- increasing alkalinity, while consumption of NH4 decreases it. 
            ! Equations R12 and R13 in Middelburg et al. (2020)
            ! PP (NO3): CO2 + n/c HNO3 + p/c H3PO4 + (1+n)H2O → (CH2O)(NH3)n(H3PO4)p + (1+2n)O2  Alk change:	p/c + n/c per mol C or p/n+1 per mol N
            ! PP (NH4): CO2 + n/c NH3 + p/c H3PO4 + H2O → (CH2O)(NH3)n(H3PO4)p + O2 Alk change: p/c-n/c per mol C or p/n-1 per mol N 
            alk_change = (1.0_rk + 1.0_rk / n_to_p) * new_prod &
                       + (-1.0_rk + 1.0_rk / n_to_p) * reg_prod

            ! Oxygen production associated with carbon fixation.
            o2_change = self%o2_per_c * self%c_to_n_phy * tpp

         else

            !-------------------------------------------------------------------------------------------
            ! Sediment layers
            !
            ! Pelagic phytoplankton entering sediment cells are rapidly removed and transferred
            ! to particulate organic matter.
            !-------------------------------------------------------------------------------------------

            phy_kill = k_sed_kill_phy * biomass

            tpp        = 0.0_rk
            new_prod   = 0.0_rk
            reg_prod   = 0.0_rk
            mort_p     = 0.0_rk

            dno3       = 0.0_rk
            dnh4       = 0.0_rk
            dic_change = 0.0_rk
            alk_change = 0.0_rk
            o2_change  = 0.0_rk

            chl_diag   = 0.0_rk

            pom_prod_n = phy_kill
            phy_source = -phy_kill

         end if

         !----------------------------------------------------------------------------------------------
         ! Tendencies
         !----------------------------------------------------------------------------------------------

         _ADD_SOURCE_(self%id_biomass, phy_source)

         _ADD_SOURCE_(self%id_no3, dno3)
         _ADD_SOURCE_(self%id_nh4, dnh4)

         if (_AVAILABLE_(self%id_dic)) then
            _ADD_SOURCE_(self%id_dic, dic_change)
         end if

         if (_AVAILABLE_(self%id_alk)) then
            _ADD_SOURCE_(self%id_alk, alk_change)
         end if

         if (_AVAILABLE_(self%id_o2)) then
            _ADD_SOURCE_(self%id_o2, o2_change)
         end if

         !----------------------------------------------------------------------------------------------
         ! Aggregate contributions
         !----------------------------------------------------------------------------------------------

         _SET_DIAGNOSTIC_(self%id_pom_prod_n, pom_prod_n)
         _SET_DIAGNOSTIC_(self%id_chl, chl_diag)

         !----------------------------------------------------------------------------------------------
         ! Optional diagnostics
         !----------------------------------------------------------------------------------------------

         if (self%save_process_rates) then

            _SET_DIAGNOSTIC_(self%id_TPP,      tpp      * secs_pr_day)
            _SET_DIAGNOSTIC_(self%id_new_prod, new_prod * secs_pr_day)
            _SET_DIAGNOSTIC_(self%id_reg_prod, reg_prod * secs_pr_day)
            _SET_DIAGNOSTIC_(self%id_o2_prod,  o2_change * secs_pr_day)

         end if

         if (self%save_alkalinity_changes) then
            _SET_DIAGNOSTIC_(self%id_alk_pp, alk_change * secs_pr_day)
         end if

      _LOOP_END_

   end subroutine do

end module phytoplankton