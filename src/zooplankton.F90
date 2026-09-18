#include "fabm_driver.h"

!-------------------------------------------------------------------------------------------------------
! Zooplankton functional group.
!
! Zooplankton grazes on a coupled phytoplankton biomass state using a Holling-II formulation.
! Grazed material is partitioned into sloppy feeding losses, egestion, zooplankton growth,
! nitrogen excretion, and carbon respiration.
!
! Particulate losses contribute to the shared pelagic POM-production aggregate.
! A fixed Redfield C:N ratio is used to convert ingested prey nitrogen to carbon.
!
! This initial implementation supports one phytoplankton prey and preserves the grazing and
! zooplankton-growth equations of the original monolithic pelagic ecosystem.
!-------------------------------------------------------------------------------------------------------
module zooplankton

   use fabm_types
   use pelagic_common, only: redfield_c_to_n, pom_production_n

   implicit none
   private

   type, extends(type_base_model), public :: type_zooplankton

      ! --- State variable
      type(type_state_variable_id) :: id_biomass

      ! --- Coupled prey state
      type(type_state_variable_id) :: id_prey

      ! --- Couplings
      type(type_state_variable_id) :: id_nh4

      ! --- Optional couplings
      type(type_state_variable_id) :: id_o2
      type(type_state_variable_id) :: id_dic
      type(type_state_variable_id) :: id_alk

      ! --- Environmental dependencies
      type(type_dependency_id) :: id_porosity

      ! --- Internal diagnostic used for aggregate POM production
      type(type_diagnostic_variable_id) :: id_pom_prod_n

      ! --- Optional diagnostics
      type(type_diagnostic_variable_id) :: id_GRAZE
      type(type_diagnostic_variable_id) :: id_o2_cons
      type(type_diagnostic_variable_id) :: id_alk_resp

      ! --- Parameters
      real(rk) :: g_max
      real(rk) :: k_g
      real(rk) :: m_zoo2

      real(rk) :: sloppy_feed
      real(rk) :: beta_n
      real(rk) :: beta_c
      real(rk) :: k_zoo

      real(rk) :: o2_per_c

      ! --- Diagnostic switches
      logical :: save_process_rates
      logical :: save_alkalinity_changes

   contains
      procedure :: initialize
      procedure :: do
   end type type_zooplankton

contains

   subroutine initialize(self, configunit)

      class(type_zooplankton), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk

      !-------------------------------------------------------------------------------------------------
      ! Parameters
      !-------------------------------------------------------------------------------------------------

      call self%get_parameter(self%g_max, 'g_max', 'd-1', 'Maximum specific grazing rate', &
                              default=0.5_rk, scale_factor=d_per_s)

      call self%get_parameter(self%k_g, 'k_g', 'mmol m-3', 'Holling-II half-saturation prey biomass', &
                              default=0.5_rk)

      call self%get_parameter(self%m_zoo2, 'm_zoo2', '(mmol m-3)-1 d-1', 'Quadratic zooplankton mortality', &
                              default=0.2_rk, scale_factor=d_per_s)

      call self%get_parameter(self%sloppy_feed, 'sloppy_feed', '-', &
                              'Fraction of grazed material lost as sloppy feeding to POM', &
                              default=0.20_rk, minimum=0.0_rk, maximum=1.0_rk)

      call self%get_parameter(self%beta_n, 'beta_n', '-', 'Zooplankton nitrogen assimilation efficiency', &
                              default=0.77_rk, minimum=0.0_rk, maximum=1.0_rk)

      call self%get_parameter(self%beta_c, 'beta_c', '-', 'Zooplankton carbon assimilation efficiency', &
                              default=0.64_rk, minimum=0.0_rk, maximum=1.0_rk)

      call self%get_parameter(self%k_zoo, 'k_zoo', '-', 'Zooplankton net carbon growth efficiency', &
                              default=0.80_rk, minimum=0.0_rk, maximum=1.0_rk)

      call self%get_parameter(self%o2_per_c, 'o2_per_c', 'mol O2 mol C-1', &
                              'Effective O2 consumed per C respired', default=1.3_rk)

      call self%get_parameter(self%save_process_rates, 'process_rates', '', &
                              'Save process-rate diagnostics', default=.false.)

      call self%get_parameter(self%save_alkalinity_changes, 'alkalinity_changes', '', &
                              'Save alkalinity-change diagnostics', default=.false.)

      !-------------------------------------------------------------------------------------------------
      ! State variable
      !-------------------------------------------------------------------------------------------------

      call self%register_state_variable(self%id_biomass, 'biomass', 'mmol N m-3', 'Biomass', &
                                        1.0e-12_rk, minimum=0.0_rk)

      call self%set_variable_property(self%id_biomass, 'is_solute', .false.)

      ! Zooplankton contributes to total nitrogen.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_biomass)

      !-------------------------------------------------------------------------------------------------
      ! Couplings
      !-------------------------------------------------------------------------------------------------

      ! Grazing target. The zooplankton module both reads this state and applies the grazing loss to it.
      call self%register_state_dependency(self%id_prey, 'prey_biomass', 'mmol N m-3', &
                                          'Phytoplankton prey biomass', required=.true.)

      ! Nitrogen excretion is returned to ammonium.
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mmol m-3', &
                                          'Dissolved ammonium', required=.true.)

      call self%register_state_dependency(self%id_o2, 'o2', 'mmol m-3', &
                                          'Dissolved oxygen', required=.false.)

      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C m-3', &
                                          'Total dissolved inorganic carbon', required=.false.)

      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', &
                                          'Total alkalinity', required=.false.)

      !-------------------------------------------------------------------------------------------------
      ! Environmental dependencies
      !-------------------------------------------------------------------------------------------------

      call self%register_dependency(self%id_porosity, &
                                    type_interior_standard_variable(name='porosity', units='1'), &
                                    required=.false.)

      !-------------------------------------------------------------------------------------------------
      ! Internal diagnostic contributing to shared aggregate POM production
      !-------------------------------------------------------------------------------------------------

      call self%register_diagnostic_variable(self%id_pom_prod_n, 'pom_prod_n', &
                                             'mmol N m-3 s-1', &
                                             'Zooplankton contribution to POM production', &
                                             output=output_none)

      call self%add_to_aggregate_variable(pom_production_n, self%id_pom_prod_n)

      !-------------------------------------------------------------------------------------------------
      ! Optional process diagnostics
      !-------------------------------------------------------------------------------------------------

      if (self%save_process_rates) then

         call self%register_diagnostic_variable(self%id_GRAZE, 'GRAZE', 'mmol N m-3 d-1', &
                                                'Grazing rate (N-based)')

         call self%register_diagnostic_variable(self%id_o2_cons, 'O2_CONS', 'mmol m-3 d-1', &
                                                'O2 consumption from zooplankton respiration')

      end if

      if (self%save_alkalinity_changes) then

         call self%register_diagnostic_variable(self%id_alk_resp, 'ALK_RESP', 'mmol eq m-3 d-1', &
                                                'Alkalinity change due to zooplankton excretion')

      end if

   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)

      class(type_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: biomass
      real(rk) :: prey
      real(rk) :: phi

      real(rk) :: graze
      real(rk) :: ing_n, ing_c
      real(rk) :: sloppy_n
      real(rk) :: assim_n, assim_c
      real(rk) :: egestion_n

      real(rk) :: grow_n_pot
      real(rk) :: grow_c_pot
      real(rk) :: zoo_growth

      real(rk) :: zoo_excr_n
      real(rk) :: zoo_resp_c
      real(rk) :: mortality

      real(rk) :: biomass_source
      real(rk) :: prey_source
      real(rk) :: biomass_kill

      real(rk) :: nh4_change
      real(rk) :: dic_change
      real(rk) :: alk_change
      real(rk) :: o2_change

      real(rk) :: pom_prod_n
      real(rk) :: o2_cons

      logical :: is_sediment

      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: k_sed_kill_zoo = 200.0_rk / secs_pr_day
      real(rk), parameter :: n_to_p = 16.0_rk
      real(rk), parameter :: eps = 1.0e-12_rk

      _LOOP_BEGIN_

         !----------------------------------------------------------------------------------------------
         ! State and environmental conditions
         !----------------------------------------------------------------------------------------------

         _GET_(self%id_biomass, biomass)
         _GET_(self%id_prey, prey)

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

            ! Zooplankton consumes phytoplankton with the same Holling-II response
            ! used in the original monolithic pelagic ecosystem.
            graze = self%g_max * biomass * prey / (self%k_g + max(prey, 0.0_rk))

            ! Grazed prey expressed in nitrogen and carbon units.
            ing_n = graze
            ing_c = redfield_c_to_n * graze

            ! A fraction of ingestion is lost directly as particulate matter.
            sloppy_n = self%sloppy_feed * ing_n

            ! Remaining ingested material is available for assimilation.
            assim_n = max(0.0_rk, (1.0_rk - self%sloppy_feed) * ing_n)
            assim_c = max(0.0_rk, (1.0_rk - self%sloppy_feed) * ing_c)

            ! Non-assimilated nitrogen is egested as particulate matter.
            egestion_n = (1.0_rk - self%beta_n) * assim_n

            ! Potential growth supported independently by assimilated N and C.
            grow_n_pot = self%beta_n * assim_n
            grow_c_pot = self%beta_c * self%k_zoo * assim_c / redfield_c_to_n

            ! Actual zooplankton growth is limited by the more restrictive element.
            zoo_growth = min(grow_n_pot, grow_c_pot)

            ! Assimilated material not incorporated into zooplankton biomass is returned
            ! to dissolved pools as NH4 and DIC.
            zoo_excr_n = max(0.0_rk, self%beta_n * assim_n - zoo_growth)
            zoo_resp_c = max(0.0_rk, self%beta_c * assim_c - redfield_c_to_n * zoo_growth)

            ! Quadratic mortality contributes to POM.
            mortality = self%m_zoo2 * biomass * biomass

            biomass_source = zoo_growth - mortality
            prey_source = -graze

            pom_prod_n = sloppy_n + egestion_n + mortality

            nh4_change = zoo_excr_n
            dic_change = zoo_resp_c

            ! Respiration/excretion alkalinity term retained from the original model.
            alk_change = zoo_excr_n - zoo_excr_n / n_to_p

            o2_cons = self%o2_per_c * zoo_resp_c
            o2_change = -o2_cons

         else

            !-------------------------------------------------------------------------------------------
            ! Sediment layers
            !
            ! Pelagic zooplankton entering sediment cells are rapidly removed and transferred
            ! to particulate organic matter.
            !-------------------------------------------------------------------------------------------

            biomass_kill = k_sed_kill_zoo * biomass

            graze         = 0.0_rk
            zoo_growth    = 0.0_rk
            zoo_excr_n    = 0.0_rk
            zoo_resp_c    = 0.0_rk
            mortality     = 0.0_rk

            biomass_source = -biomass_kill
            prey_source    = 0.0_rk

            nh4_change = 0.0_rk
            dic_change = 0.0_rk
            alk_change = 0.0_rk
            o2_cons    = 0.0_rk
            o2_change  = 0.0_rk

            pom_prod_n = biomass_kill

         end if

         !----------------------------------------------------------------------------------------------
         ! Tendencies
         !----------------------------------------------------------------------------------------------

         _ADD_SOURCE_(self%id_biomass, biomass_source)
         _ADD_SOURCE_(self%id_prey, prey_source)
         _ADD_SOURCE_(self%id_nh4, nh4_change)

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
         ! Aggregate contribution
         !----------------------------------------------------------------------------------------------

         _SET_DIAGNOSTIC_(self%id_pom_prod_n, pom_prod_n)

         !----------------------------------------------------------------------------------------------
         ! Optional diagnostics
         !----------------------------------------------------------------------------------------------

         if (self%save_process_rates) then
            _SET_DIAGNOSTIC_(self%id_GRAZE,   graze   * secs_pr_day)
            _SET_DIAGNOSTIC_(self%id_o2_cons, o2_cons * secs_pr_day)
         end if

         if (self%save_alkalinity_changes) then
            _SET_DIAGNOSTIC_(self%id_alk_resp, alk_change * secs_pr_day)
         end if

      _LOOP_END_

   end subroutine do

end module zooplankton
