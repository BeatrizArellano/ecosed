#include "fabm_driver.h"
!-------------------------------------------------------------------------------------------------------
! Phytoplankton growth is controlled by nutrient (Michaelis–Menten), light, and temperature limitation,
! producing biomass (GPP) that is consumed by zooplankton via a Holling-II grazing formulation.
! Grazed material is partitioned into sloppy feeding losses (direct to POM), egestion (faecal pellets),
! and assimilated material, which is split between zooplankton growth and respiration.
! Excess assimilated nitrogen is excreted to DIN, while excess carbon is respired to DIC.
! Mortality and feeding losses contribute to production of POM/detritus. 
!--------------------------------------------------------------------------------------------------------
module pelagic_ecosystem

   use fabm_types
   implicit none
   private

   type, extends(type_base_model), public :: type_pelagic_ecosystem
      ! --- State variable identifiers
      type(type_state_variable_id) :: id_phy, id_chl, id_zoo

      ! --- Couplings 
      type(type_state_variable_id) :: id_no3

      ! --- Optional Couplings 
      type(type_state_variable_id) :: id_dic
      type(type_state_variable_id) :: id_o2

      ! --- Environmental dependencies
      type(type_dependency_id)          :: id_par      
      type(type_dependency_id)          :: id_temp
      type(type_dependency_id)          :: id_porosity

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_GPP, id_GRAZE
      type(type_diagnostic_variable_id) :: id_o2_prod
      type(type_diagnostic_variable_id) :: id_o2_cons
      type(type_diagnostic_variable_id) :: id_pom_prod_n  ! Produced POM
      type(type_diagnostic_variable_id) :: id_chl_diag    ! When photoacclimation is disabled, chl is diagnostic

      ! --- Parameters
      real(rk) :: mu_max, k_n
      real(rk) :: g_max, k_g
      real(rk) :: m_phy, m_zoo2

      logical  :: photoacclimation
      real(rk) :: chl_per_c
      real(rk) :: alpha_phy
      real(rk) :: alpha_chl      
      real(rk) :: theta_chl_max
      ! Grazing
      real(rk) :: sloppy_feed     ! messy feeding / sloppy loss fraction
      real(rk) :: beta_n          ! N assimilation efficiency
      real(rk) :: beta_c          ! C assimilation efficiency
      real(rk) :: k_zoo           ! Fraction of assimilated carbon allocated to zooplankton growth; remainder is lost via respiration.

      real(rk) :: c_to_n_phy      ! mol C per mol N in phytoplankton
      real(rk) :: c_to_n_zoo      ! mol C per mol N in zooplankton
      real(rk) :: o2_per_c        ! mol O2 per C produced/respired

      real(rk) :: kc              ! specific light extinction [m2 mmol-1] applied to phy + POM pools

   contains
      procedure :: initialize
      procedure :: do
   end type type_pelagic_ecosystem

contains

   subroutine initialize(self, configunit)
      class(type_pelagic_ecosystem), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
      real(rk) :: w_phy

      ! ---------------- Parameters ----------------
      ! Rates are provided in d-1 in config file and scaled to s-1 here.
      call self%get_parameter(self%mu_max,'mu_max','d-1','max specific phyto growth rate', default=1.0_rk, scale_factor=d_per_s)
      call self%get_parameter(self%k_n,   'k_n',   'mmol m-3','DIN half-saturation for growth', default=0.5_rk)

      call self%get_parameter(self%photoacclimation, 'photoacclimation', '-','Enable chlorophyll photoacclimation', default=.false.)
      call self%get_parameter(self%alpha_chl,'alpha_chl','mmolC mgChl-1 (W m-2)-1 d-1', 'Chlorophyll-specific initial slope of P-I curve', default=0.7_rk, scale_factor=d_per_s)
      call self%get_parameter(self%theta_chl_max,'theta_chl_max','mgChl mmolC-1', 'Maximum chlorophyll-to-carbon ratio', default=0.6_rk)
      call self%get_parameter(self%alpha_phy, 'alpha_phy', 'm2 W-1 d-1', 'Initial slope of phytoplankton P-I curve when photoacclimation is disabled', default=0.03_rk, scale_factor=d_per_s)
      call self%get_parameter(self%chl_per_c,'chl_per_c','mgChl mmolC-1', 'Fixed Chl:C ratio when photoacclimation is disabled', default=0.3_rk)

      call self%get_parameter(self%g_max, 'g_max', 'd-1','max specific grazing rate', default=0.5_rk, scale_factor=d_per_s)
      call self%get_parameter(self%k_g,   'k_g',   'mmol m-3','Holling-II half-sat prey for grazing', default=0.5_rk)

      call self%get_parameter(self%m_phy, 'm_phy','d-1','linear phyto mortality to POM_S', default=0.05_rk, scale_factor=d_per_s)
      call self%get_parameter(self%m_zoo2,'m_zoo2','(mmol m-3)-1 d-1','quadratic zoo mortality to POM_S', default=0.2_rk, scale_factor=d_per_s)

      call self%get_parameter(self%sloppy_feed, 'sloppy_feed', '-', 'Fraction of grazed material lost as sloppy/messy feeding to POM', default=0.20_rk, minimum=0.0_rk, maximum=1.0_rk)

      call self%get_parameter(self%beta_n, 'beta_n', '-', 'Zooplankton Nitrogen assimilation efficiency', default=0.77_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%beta_c, 'beta_c', '-', 'Zooplankton Carbon assimilation efficiency',  default=0.64_rk, minimum=0.0_rk, maximum=1.0_rk)

      call self%get_parameter(self%k_zoo, 'k_zoo', '-', 'Zooplankton net carbon growth efficiency', default=0.80_rk, minimum=0.0_rk, maximum=1.0_rk)

      call self%get_parameter(self%c_to_n_phy, 'c_to_n_phy', '-', 'Phytoplankton molar C:N ratio', default=6.625_rk)
      call self%get_parameter(self%c_to_n_zoo, 'c_to_n_zoo', '-', 'Zooplankton molar C:N ratio', default=5.625_rk)
      call self%get_parameter(self%o2_per_c,  'o2_per_c',  'mol O2 mol C-1', 'Effective O2 produced/consumed per C fixed/respired', default=1.3_rk)

      call self%get_parameter(self%kc,'kc','m2 mmol-1','specific light extinction of Phytoplankton', default=0.03_rk)

      ! sinking speeds for phytoplankton (m d-1 in config, scaled to m s-1)
      call self%get_parameter(w_phy,'w_phy','m d-1','vertical velocity of Phy (<0 sinking)', default=0.0_rk, scale_factor=d_per_s)

      ! ---------------- State variables ----------------      
      call self%register_state_variable(self%id_phy,   'phy',   'mmol m-3', 'phytoplankton',   1.0e-12_rk, minimum=0.0_rk, vertical_movement=w_phy)
      call self%set_variable_property(self%id_phy, 'is_solute', .false.)
      if (self%photoacclimation) then
         call self%register_state_variable(self%id_chl,   'chl',   'mg m-3', 'chlorophyll',       1.0e-12_rk, minimum=0.0_rk, vertical_movement=w_phy)
         call self%set_variable_property(self%id_chl, 'is_solute', .false.)
      end if
      call self%register_state_variable(self%id_zoo,   'zoo',   'mmol m-3', 'zooplankton',     1.0e-12_rk, minimum=0.0_rk)
      call self%set_variable_property(self%id_zoo, 'is_solute', .false.)

      ! aggregate nitrogen (useful for assessing Conservation)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_phy)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_zoo)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_no3, 'no3', 'mmol m-3', 'Dissolved Nitrate', required=.true.)
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol m-3', 'Total dissolved inorganic carbon', required=.false.)
      call self%register_state_dependency(self%id_o2,  'o2',  'mmol m-3', 'dissolved oxygen',          required=.false.)

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_GPP,   'GPP',   'mmol m-3 d-1', 'Gross Primary Production (N-based)')
      call self%register_diagnostic_variable(self%id_GRAZE, 'GRAZE', 'mmol m-3 d-1', 'Grazing rate (N-based)')
      call self%register_diagnostic_variable(self%id_o2_prod, 'O2_PROD', 'mmol m-3 d-1', 'O2 production from primary production')
      call self%register_diagnostic_variable(self%id_o2_cons, 'O2_CONS', 'mmol m-3 d-1', 'Total O2 consumption')
      call self%register_diagnostic_variable(self%id_pom_prod_n, 'pom_prod_n', 'mmol m-3 s-1', 'Production rate of particulate organic matter from pelagic biology')
      if (.not. self%photoacclimation) then
         call self%register_diagnostic_variable(self%id_chl_diag, 'chl', 'mg m-3','chlorophyll diagnosed from phytoplankton biomass')
      end if

      ! ---------------- Environmental dependencies ----------------
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_porosity, type_interior_standard_variable(name='porosity', units='1'), required=.false.)

      ! light attenuation feedback (PHY + POM)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_phy,   scale_factor=self%kc)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_pelagic_ecosystem), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: phy, chl, zoo
      real(rk)            :: chl_diag
      real(rk)            :: no3
      real(rk)            :: par, temp, phi

      real(rk)            :: fN, fT, mu_T, alphaI
      real(rk)            :: theta, cphy, fch, jden
      real(rk)            :: chl_prod, chl_loss
      real(rk)            :: phy_source
      real(rk)            :: gpp, graze
      real(rk)            :: mort_p, mort_z
      real(rk)            :: dno3
      real(rk)            :: dic_change, o2_change
      real(rk)            :: ing_n, ing_c
      real(rk)            :: sloppy_n
      real(rk)            :: assim_n, assim_c
      real(rk)            :: egestion_n
      real(rk)            :: grow_n_pot, grow_c_pot
      real(rk)            :: zoo_growth
      real(rk)            :: zoo_excr_n
      real(rk)            :: zoo_resp_c
      real(rk)            :: pom_prod_n
      real(rk)            :: o2_prod, o2_cons_resp
      real(rk)            :: phy_kill, zoo_kill, chl_kill
      logical             :: is_sediment

      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: k_sed_kill_phy = 200.0_rk * (1.0_rk/secs_pr_day)        ! Rapid removal of phytoplankton in sediments
      real(rk), parameter :: k_sed_kill_zoo = 200.0_rk * (1.0_rk/secs_pr_day)        ! Rapid removal of zooplankton in sediments
      real(rk), parameter :: eps = 1.0e-12_rk

      _LOOP_BEGIN_

         _GET_(self%id_no3, no3)
         _GET_(self%id_phy, phy)
         if (self%photoacclimation) then
            _GET_(self%id_chl, chl)
         else
            chl = 0.0_rk
         end if
         _GET_(self%id_zoo, zoo)
         ! Environment
         _GET_(self%id_par, par)
         _GET_(self%id_temp, temp)

         ! Using porosity to differentiate water from sediments
         if (_AVAILABLE_(self%id_porosity)) then
            _GET_(self%id_porosity, phi)
         else 
            phi = 1.0_rk                        ! Porosity is 1 for water
         end if
         ! phi=1 in water layers and phi<1 in sediment layers.
         is_sediment = (phi < 1.0_rk - eps)

         ! ==============================================================
         ! In sediment layers: kill pelagic biomass rapidly and transfer
         ! phy + zoo nitrogen to particulate organic matter.
         ! Chlorophyll is removed too, but does not contribute to POM.
         ! ==============================================================
         if (is_sediment) then

            phy_kill = k_sed_kill_phy * phy
            zoo_kill = k_sed_kill_zoo * zoo

            if (self%photoacclimation) then
               chl_kill = k_sed_kill_phy * chl
            else
               chl_kill = 0.0_rk
               chl_diag = 0.0_rk
            end if

            ! No active pelagic biology in sediments
            gpp         = 0.0_rk
            graze       = 0.0_rk
            mort_p      = 0.0_rk
            mort_z      = 0.0_rk
            zoo_growth  = 0.0_rk
            zoo_excr_n  = 0.0_rk
            zoo_resp_c  = 0.0_rk
            dno3        = 0.0_rk
            dic_change  = 0.0_rk
            o2_prod     = 0.0_rk
            o2_cons_resp= 0.0_rk
            o2_change   = 0.0_rk
            chl_prod    = 0.0_rk
            chl_loss    = 0.0_rk

            ! Dead pelagic biomass becomes POM (N-based)
            pom_prod_n = phy_kill + zoo_kill

            ! Apply tracer tendencies
            _ADD_SOURCE_(self%id_phy, -phy_kill)
            _ADD_SOURCE_(self%id_zoo, -zoo_kill)

            if (self%photoacclimation) then
               _ADD_SOURCE_(self%id_chl, -chl_kill)
            end if

         else
         ! In the water column   
            !-------------------------------------------------------------------
            !                 Phytoplankton and Chlorophyll
            !-------------------------------------------------------------------

            ! --- Nutrient limitation ---
            ! Michaelis-Menten uptake: growth increases with DIN and saturates at high values.
            ! fN controls the strength of nutrient limitation (0–1 scaling).
            fN = max(no3,0.0_rk) / (self%k_n + max(no3,0.0_rk))

            ! --- Temperature dependence ---
            ! Eppley's (1972) exponential scaling of metabolic rates with temperature.
            ! Sets the temperature-dependent maximum growth rate.
            fT = 1.066_rk ** temp         
            mu_T = self%mu_max * fT    ! Temperature-adjusted maximum growth rate

            ! Convert N-based phytoplankton biomass to carbon.
            cphy = self%c_to_n_phy * max(phy,0.0_rk)  ! Carbon in phytoplankton

            ! --- Light limitation ---
            chl_loss = 0.0_rk
            chl_prod = 0.0_rk
            chl_diag = 0.0_rk
            ! In photoacclimation mode, chlorophyll is prognostic and defines Chl:C.
            ! Otherwise chlorophyll is diagnosed from a fixed Chl:C ratio.
            if (self%photoacclimation) then
               ! Compute current Chl:C ratio from chlorophyll
               if (cphy > eps) then
                  theta = chl / cphy                     ! Chl:C ratio
               else
                  theta = self%theta_chl_max             ! Maximum Chl:C ratio
               end if
               theta = min(self%theta_chl_max, max(theta, eps))  
               ! Chlorophyll-dependent light harvesting; proportional to Chl:C ratio.
               ! alphaI represents the light-limited potential growth rate.
               alphaI = self%alpha_chl * theta * max(par,0.0_rk)
            else
               ! Simple non-chlorophyll light limitation
               alphaI = self%alpha_phy * max(par,0.0_rk)
               ! No independent chlorophyll production in this mode
               chl_diag = self%chl_per_c * cphy
               chl_prod = 0.0_rk
            end if

            ! --- Primary production ---
            ! Primary production (N-based) follows a smooth saturating light-response formulation,
            ! transitioning between light-limited (alphaI) and temperature-limited (mu_T) growth.
            ! Resulting N-based production is further reduced by nutrient limitation.
            jden = sqrt(mu_T*mu_T + alphaI*alphaI + eps)
            gpp = (mu_T * alphaI / jden) * fN * phy

            ! --- Chlorophyll photoacclimation ---
            ! Chlorophyll synthesis tied to carbon fixation, enhanced under low light.
            ! Regulates Chl:C ratio toward a light-dependent optimum, limited by theta_chl_max.
            if (self%photoacclimation) then
               fch = mu_T / jden                                           ! Light-control term for chlorophyll synthesis
               ! Chl production enhanced under low-light while consrained by nutrient availability
               chl_prod = (self%theta_chl_max * fch * fN / theta) * (gpp * self%c_to_n_phy) 
            end if

            !-------------------------------------------------------------------
            !             Zooplankton grazing
            !-------------------------------------------------------------------
            
            ! Zooplankton consumes phytoplankton with a saturating (Holling II) response
            ! increases with prey and saturates at high phy (g_max * zoo)
            graze = self%g_max * zoo * phy / (self%k_g + max(phy, 0.0_rk))

            ! --- Ingestion ---
            ! Grazed phytoplankton is converted into ingested nitrogen and carbon
            ! Grazed material expressed in N and C units
            ing_n = graze
            ing_c = self%c_to_n_phy * graze

            ! --- Sloppy feeding ---
            ! A fraction of ingestion is lost directly as particles (messy feeding)
            sloppy_n = self%sloppy_feed * ing_n         

            ! --- Assimilation ---
            ! Remaining ingested material is available for assimilation by zooplankton
            assim_n = max(0.0_rk, (1.0_rk - self%sloppy_feed) * ing_n)
            assim_c = max(0.0_rk, (1.0_rk - self%sloppy_feed) * ing_c)

            ! Non-assimilated nitrogen from ingestion is lost to particulate matter (faecal pellets)
            egestion_n = (1.0_rk - self%beta_n) * assim_n

            ! --- Growth limitation ---
            ! Potential growth from nitrogen and carbon constraints
            grow_n_pot = self%beta_n * assim_n
            grow_c_pot = self%beta_c * self%k_zoo * assim_c / self%c_to_n_zoo

            ! Actual zooplankton growth is limited by the most restrictive element (C or N)
            zoo_growth = min(grow_n_pot, grow_c_pot)

            !-------------------------------------------------------------------
            !             Biomass loss
            !-------------------------------------------------------------------

            ! ---  Mortality ----
            ! Mortality losses to particulate organic matter
            ! Linear phytoplankton mortality: constant fractional loss to POM
            mort_p = self%m_phy * phy     
            ! Quadratic zooplankton mortality: density-dependent loss to POM
            mort_z = self%m_zoo2 * zoo * zoo

            ! --- Chlorophyll loss ---
            ! Chlorophyll decreases in proportion to phytoplankton biomass losses
            phy_source = gpp - graze - mort_p 
            if (self%photoacclimation) then
               if (phy > eps) then
                  chl_loss = (chl/phy) * (graze + mort_p)
               else
                  chl_loss = 0.0_rk
               end if
            end if

            ! --- Excretion from assimilated biomass ---
            ! Material assimilated but not used for growth is returned to dissolved pools
            ! Excess assimilated N is excreted as DIN; excess assimilated C is respired as DIC
            zoo_excr_n = max(0.0_rk, self%beta_n * assim_n - zoo_growth)
            zoo_resp_c = max(0.0_rk, self%beta_c * assim_c - self%c_to_n_zoo * zoo_growth)

            ! Total particulate matter is produced from feeding losses, egestion, and mortality
            pom_prod_n = sloppy_n + egestion_n + mort_p + mort_z

            !--------------------------------------------------------------------------
            !               Tendencies (all in s-1) 
            !--------------------------------------------------------------------------
            ! Changes in NO3
            dno3 = -gpp + zoo_excr_n 
            ! Changes in DIC
            dic_change = -self%c_to_n_phy * gpp + zoo_resp_c           
            ! Changes in O2
            o2_prod      = self%o2_per_c * self%c_to_n_phy * gpp           
            o2_cons_resp = self%o2_per_c * zoo_resp_c
            o2_change    = o2_prod - o2_cons_resp

            ! NO3         
            _ADD_SOURCE_(self%id_no3, dno3)

            ! PHY
            _ADD_SOURCE_(self%id_phy, phy_source)
            if (self%photoacclimation) then
               _ADD_SOURCE_(self%id_chl, chl_prod - chl_loss)
            end if

            ! ZOO
            _ADD_SOURCE_(self%id_zoo, zoo_growth - mort_z)

            ! ---------------- Optional couplings ----------------
            ! Carbon and oxygen are affected by primary production and zooplankton respiration
            if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic, dic_change)
            if (_AVAILABLE_(self%id_o2))  _ADD_SOURCE_(self%id_o2,  o2_change)
         end if

         ! ---------------- Diagnostics (convert to rates per day) ----------------
         _SET_DIAGNOSTIC_(self%id_GPP,       gpp*secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_GRAZE,     graze*secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_o2_prod,   o2_prod * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_o2_cons,   o2_cons_resp * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_pom_prod_n, pom_prod_n)
         if (.not. self%photoacclimation) then
            _SET_DIAGNOSTIC_(self%id_chl_diag, chl_diag)
         end if

      _LOOP_END_
   end subroutine do

end module pelagic_ecosystem