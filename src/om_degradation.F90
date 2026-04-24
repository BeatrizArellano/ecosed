#include "fabm_driver.h"

! Organic matter degradation.
! POM pools are stored in N units and partitioned into labile, semi-labile, and refractory fractions.
module om_degradation

   use fabm_types
   implicit none
   private

   type, extends(type_base_model), public :: type_om_degradation
      ! --- State variables
      type(type_state_variable_id) :: id_pom_l, id_pom_s, id_pom_r

      ! --- Couplings
      type(type_state_variable_id) :: id_no3, id_nh4
      type(type_state_variable_id) :: id_po4
      type(type_state_variable_id) :: id_o2

      !--- Optional Coupling
      type(type_state_variable_id) :: id_dic, id_alk      

      ! --- Dependencies
      type(type_dependency_id)     :: id_temp
      type(type_dependency_id)     :: id_porosity
      type(type_dependency_id)     :: id_pom_prod_n

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_rem
      type(type_diagnostic_variable_id) :: id_total_pom
      type(type_diagnostic_variable_id) :: id_o2_cons_rem
      type(type_diagnostic_variable_id) :: id_n_loss_denit

      ! --- Parameters
      real(rk) :: frac_lab, frac_semi, frac_ref
      real(rk) :: k_remin_pom_l, k_remin_pom_s, k_remin_pom_r  ! Remineralisation rates
      real(rk) :: c_to_n_pom_l
      real(rk) :: c_to_n_pom_s
      real(rk) :: c_to_n_pom_r
      real(rk) :: n_to_p_pom_l
      real(rk) :: n_to_p_pom_s
      real(rk) :: n_to_p_pom_r
      real(rk) :: o2_per_c
      real(rk) :: k_o2_aer
      real(rk) :: ki_o2_redox
      real(rk) :: o2_thr_redox
      real(rk) :: k_no3_denit 
      real(rk) :: atten

   contains
      procedure :: initialize
      procedure :: do
   end type type_om_degradation

contains

   subroutine initialize(self, configunit)
      class(type_om_degradation), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
      real(rk), parameter :: eps     = 1.0e-12_rk
      real(rk) :: w_pom

      ! ---------------- Parameters ----------------
      call self%get_parameter(self%frac_lab,  'frac_lab',  '-', 'Fraction of newly produced labile organic matter', default=0.5_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%frac_semi, 'frac_semi', '-', 'Fraction of newly produced semi-labile organic matter', default=0.4_rk, minimum=0.0_rk, maximum=1.0_rk)

      if (self%frac_lab + self%frac_semi .gt. 1.0_rk) then
         write(*,*) 'yaml_settings: frac_lab + frac_semi must be <= 1'
         stop
      end if
      self%frac_ref = max(0.0_rk, 1.0_rk - self%frac_lab - self%frac_semi)

      call self%get_parameter(self%k_remin_pom_l,'k_remin_pom_l','d-1','Labile OM remineralisation rate', default=0.1_rk, scale_factor=d_per_s)
      call self%get_parameter(self%k_remin_pom_s,'k_remin_pom_s','d-1','Semi-labile OM remineralisation rate', default=0.01_rk, scale_factor=d_per_s)
      call self%get_parameter(self%k_remin_pom_r,'k_remin_pom_r','d-1','Refractory OM remineralisation rate', default=0.001_rk, scale_factor=d_per_s)

      call self%get_parameter(self%c_to_n_pom_l, 'c_to_n_l', '-', 'Molar C:N ratio of labile POM', default=6.625_rk)
      call self%get_parameter(self%c_to_n_pom_s, 'c_to_n_s', '-', 'Molar C:N ratio of semi-labile POM', default=8.0_rk)
      call self%get_parameter(self%c_to_n_pom_r, 'c_to_n_r', '-', 'Molar C:N ratio of refractory POM', default=10.0_rk)
      call self%get_parameter(self%n_to_p_pom_l, 'n_to_p_l', '-', 'Molar N:P ratio of labile POM', default=16.0_rk)
      call self%get_parameter(self%n_to_p_pom_s, 'n_to_p_s', '-', 'Molar N:P ratio of semi-labile POM', default=16.0_rk)
      call self%get_parameter(self%n_to_p_pom_r, 'n_to_p_r', '-', 'Molar N:P ratio of refractory POM', default=16.0_rk)

      call self%get_parameter(self%k_o2_aer, 'k_o2_aer', 'mmol m-3', 'Half-saturation constant for O2 limitation in aerobic remineralisation', default=3.0_rk, minimum=eps)
      call self%get_parameter(self%ki_o2_redox, 'ki_o2_redox', 'mmol m-3', 'Half-saturation constant for O2 inhibition in suboxic remineralisation pathways', default=10.0_rk, minimum=eps)
      call self%get_parameter(self%o2_thr_redox, 'o2_thr_redox', 'mmol m-3', 'O2 threshold below which suboxic/anoxic remineralisation pathways are allowed', default=10.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%k_no3_denit, 'k_no3_denit', 'mmol m-3', 'Half-saturation constant for NO3 limitation in denitrification', default=30.0_rk, minimum=eps)

      call self%get_parameter(self%o2_per_c, 'o2_per_c', 'mol O2 mol C-1', 'Effective O2 consumed per mol organic C remineralised aerobically', default=1.3_rk)
      call self%get_parameter(self%atten, 'atten', 'm2 mmol-1', 'specific light extinction of POM', default=0.03_rk)
      call self%get_parameter(w_pom,'w_pom','m d-1','vertical velocity of POM (<0 sinking)', default=-10.0_rk, maximum=0.0_rk, scale_factor=d_per_s)

      ! ---------------- State variables ----------------
      call self%register_state_variable(self%id_pom_l, 'pom_l', 'mmol N m-3', 'POM labile',      0.001_rk, minimum=0.0_rk, vertical_movement=w_pom)
      call self%register_state_variable(self%id_pom_s, 'pom_s', 'mmol N m-3', 'POM semi-labile', 0.001_rk, minimum=0.0_rk, vertical_movement=w_pom)
      call self%register_state_variable(self%id_pom_r, 'pom_r', 'mmol N m-3', 'POM refractory',  0.001_rk, minimum=0.0_rk, vertical_movement=w_pom)
      call self%set_variable_property(self%id_pom_l, 'is_solute', .false.)
      call self%set_variable_property(self%id_pom_s, 'is_solute', .false.)
      call self%set_variable_property(self%id_pom_r, 'is_solute', .false.)

      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_pom_l)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_pom_s)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_pom_r)

      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_pom_l, scale_factor=self%atten)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_pom_s, scale_factor=self%atten)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_pom_r, scale_factor=self%atten)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_no3, 'no3', 'mmol m-3', 'Dissolved nitrate', required=.true.)
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mmol m-3', 'Dissolved ammonium', required=.true.)
      call self%register_state_dependency(self%id_po4, 'po4', 'mmol P m-3', 'Dissolved phosphate', required=.true.)
      call self%register_state_dependency(self%id_o2,  'o2',  'mmol m-3', 'Dissolved oxygen', required=.true.)
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol m-3', 'Total dissolved inorganic carbon', required=.false.)
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.false.) 

      ! ---------------- Dependencies ----------------
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_porosity, type_interior_standard_variable(name='porosity', units='1'), required=.false.)
      call self%register_dependency(self%id_pom_prod_n, 'pom_prod_n', 'mmol m-3 s-1', 'Production of particulate organic Matter from pelagic biology')

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_rem, 'REM', 'mmol m-3 d-1', 'Total remineralisation rate of particulate matter.')
      call self%register_diagnostic_variable(self%id_total_pom, 'total_pom', 'mmol m-3', 'Total particulate organic matter')
      call self%register_diagnostic_variable(self%id_o2_cons_rem, 'O2_CONS_REM', 'mmol m-3 d-1', 'O2 consumption by aerobic remineralisation')
      call self%register_diagnostic_variable(self%id_n_loss_denit, 'NLOSS_DENIT', 'mmol m-3 d-1', 'N loss to N2 by denitrification')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_om_degradation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: pom_l, pom_s, pom_r
      real(rk) :: temp, pom_prod_n
      real(rk) :: phi, phi_s                     ! Porosity
      real(rk) :: p2d, d2p
      real(rk) :: o2, no3, fT
      logical  :: is_water
      ! Pathway limiting/inhibition functions
      real(rk) :: fi_o2_redox
      real(rk) :: fsum, faer, fdenit

      ! Pathway-specific remineralisation rates
      real(rk) :: rem_pom_l_aer, rem_pom_s_aer, rem_pom_r_aer
      real(rk) :: rem_pom_l_denit, rem_pom_s_denit, rem_pom_r_denit
      real(rk) :: rem_pom_l_tot, rem_pom_s_tot, rem_pom_r_tot
      real(rk) :: rem_aer_total, rem_denit_total, rem_total

      ! Stoichiometric products and sinks
      real(rk) :: dic_prod_rem
      real(rk) :: o2_cons_rem
      real(rk) :: no3_cons_denit
      real(rk) :: n_loss_denit
      real(rk) :: nh4_prod
      real(rk) :: po4_prod_rem
      real(rk) :: alk_prod_aer, alk_prod_denit, alk_prod_tot

      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: eps_phi     = 1.0e-7_rk      

      _LOOP_BEGIN_

         _GET_(self%id_pom_l, pom_l)
         _GET_(self%id_pom_s, pom_s)
         _GET_(self%id_pom_r, pom_r)
         _GET_(self%id_temp, temp)
         _GET_(self%id_pom_prod_n, pom_prod_n)
         _GET_(self%id_o2, o2)
         _GET_(self%id_no3, no3)

         !--------------------------------------------------------------------
         !  Conversion factors for cross-phase reactions (corrected by porosity)
         !--------------------------------------------------------------------
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
            is_water = .false.       ! It is a sediment layer
            p2d = phi_s / phi
            d2p = phi / phi_s            
         else
            is_water = .true.        ! It is a water-column layer
            p2d = 1.0_rk             ! No phase conversion in water
            d2p = 1.0_rk
         end if
         !------------------------------------------------------------------------

         ! Account for temperature effects on degradation rates in the water column only.    
         if (is_water) then
            ! Eppley's (1972) exponential scaling of metabolic rates with temperature.
            fT = 1.066_rk ** temp
         else
            fT = 1.0_rk
         end if

         !-------------------------------------------------------------------------
         ! Limitation and inhibition functions
         !-------------------------------------------------------------------------

         ! O2 limitation for aerobic remineralisation.
         faer = o2 / (self%k_o2_aer + o2)

         ! O2 inhibition of suboxic/anoxic remineralisation pathways.
         ! ki_o2_redox controls the smooth inhibition below the redox threshold.
         ! o2_thr_redox is an explicit threshold: above it, suboxic/anoxic pathways are inhibited.
         if (o2 <= self%o2_thr_redox) then
            fi_o2_redox = self%ki_o2_redox / (self%ki_o2_redox + o2)
         else
            fi_o2_redox = 0.0_rk
         end if

         ! NO3 limitation for denitrification and inhibition by O2.
         ! Denitrification requires both low-O2 conditions and NO3 availability.
         fdenit = fi_o2_redox * no3 / (self%k_no3_denit + no3)

         ! Combined activity (demand) of all remineralisation pathways under current redox conditions.
         fsum = faer + fdenit

         ! If the combined activity exceeds the maximum potential, rescale pathway factors so that
         ! total remineralisation does not exceed the intrinsic first-order degradation potential (k * POM).
         ! NOTE: Always rescaling (for any fsum) forces fsum=1, so total remineralisation
         !       remains fixed at k * POM and oxidants only control pathway partitioning, but not its magnitude.
         !       Here, rescaling is applied only when fsum > 1, allowing oxidant limitation to
         !       determine pathway-specific rates when the maximum potential is not exceeded.
         !       This avoids artificially increasing rates when oxidant concentrations are low.
         if (fsum > 1.0_rk) then
            ! After rescaling, the limiting functions represent fractions of total remineralisation.
            faer   = faer   / fsum
            fdenit = fdenit / fsum
         end if         

         !-------------------------------------------------------------------------
         !   Aerobic remineralisation (first-order decay process).
         !-------------------------------------------------------------------------
         ! Aerobic remineralisation rate for each POM pool, scaled by temperature and O2 availability.
         rem_pom_l_aer   = self%k_remin_pom_l * fT * faer * pom_l
         rem_pom_s_aer   = self%k_remin_pom_s * fT * faer * pom_s
         rem_pom_r_aer   = self%k_remin_pom_r * fT * faer * pom_r

         ! O2 consumption associated with aerobic remineralisation.
         ! O2 consumption is computed from the amount of organic C remineralised, using o2_per_c.
         ! Theoretically, 1 mol of O2 is consumed per mol C remineralised, but observations suggest the ratio can be higher ~1.2.
         ! Only the aerobic pathway consumes O2.
         o2_cons_rem = self%o2_per_c * (self%c_to_n_pom_l * rem_pom_l_aer &
                                      + self%c_to_n_pom_s * rem_pom_s_aer &
                                      + self%c_to_n_pom_r * rem_pom_r_aer)
         
         ! Convert O2 demand from solid-phase to porewater-volume units,
         ! because remineralisation is computed from POM concentrations per solid volume,
         ! whereas O2 is expressed per porewater volume.
         o2_cons_rem = p2d * o2_cons_rem

         ! Alkalinity generation during aerobic remineralisation.
         ! Equation R1 in Middelburg et al. (2020)
         ! (CH2O)(NH3)n/c(H3PO4)p/c + O2 → CO2 + n/c NH3 + p/c H3PO4 + H2O 	Alk change: n/c-p/c per mol C or 1-p/n per mol N         
         alk_prod_aer = rem_pom_l_aer * (1.0_rk - 1.0_rk/self%n_to_p_pom_l)  &
                      + rem_pom_s_aer * (1.0_rk - 1.0_rk/self%n_to_p_pom_s)  &
                      + rem_pom_r_aer * (1.0_rk - 1.0_rk/self%n_to_p_pom_r)

         ! ------------------------------------------------------------------------
         !    Denitrification 
         ! ------------------------------------------------------------------------
         ! Denitrification rate for each POM pool, limited by NO3 and inhibited by O2.
         rem_pom_l_denit = self%k_remin_pom_l * fdenit * pom_l
         rem_pom_s_denit = self%k_remin_pom_s * fdenit * pom_s
         rem_pom_r_denit = self%k_remin_pom_r * fdenit * pom_r

         !rem_pom_l_denit = 0.0_rk
         !rem_pom_s_denit = 0.0_rk
         !rem_pom_r_denit = 0.0_rk

         ! NO3 consumption during denitrification 
         ! CH2O(NH3)x(H3PO4)y + 0.8 NO3- -> 0.2 CO2 + 0.4 N2 + 0.8 HCO3- + x NH3 + y H3PO4 + 0.6 H2O
         ! The factor 0.8 gives the NO3 demand per mol organic C remineralised by denitrification.
         ! Because POM is stored in N units, NO3 consumption is scaled by each pool C:N ratio.
         no3_cons_denit = 0.8_rk * (self%c_to_n_pom_l * rem_pom_l_denit &
                                  + self%c_to_n_pom_s * rem_pom_s_denit &
                                  + self%c_to_n_pom_r * rem_pom_r_denit)

         ! Convert NO3 consumption from solid-phase to porewater-volume units,
         ! because denitrification is computed from POM concentrations per solid volume,
         ! whereas NO3 is expressed per porewater volume.
         no3_cons_denit = p2d * no3_cons_denit

         ! Nitrogen loss to N2 is equal to the nitrate-N consumed by denitrification.
         ! For every 0.8 mol of NO3, 0.4 mol of N2 are produced, equivalent to a loss of 0.8 mol N.
         n_loss_denit = no3_cons_denit 

         ! Alkalinity generation by denitrification
         ! Equation R2 in Middelburg et al. (2020)
         ! (CH2O)(NH3)n/c(H3PO4)p/c + 0.8 HNO3 → CO2 + n/c NH3 + p/c H3PO4 + 0.4N2 + 1.4H2O 	Alk: 0.8+n-p per mol C or c/n(0.8)+1-p/n per mol N
         alk_prod_denit = rem_pom_l_denit * (0.8_rk*self%c_to_n_pom_l + 1.0_rk - 1.0_rk/self%n_to_p_pom_l)  &
                        + rem_pom_s_denit * (0.8_rk*self%c_to_n_pom_s + 1.0_rk - 1.0_rk/self%n_to_p_pom_s)  &
                        + rem_pom_r_denit * (0.8_rk*self%c_to_n_pom_r + 1.0_rk - 1.0_rk/self%n_to_p_pom_r) 
         
         ! ------------------------------------------------------------------------
         !    Total remineralisation rates and products
         ! ------------------------------------------------------------------------
         ! Total remineralisation rates summed across POM pools and pathways.
         rem_aer_total   = rem_pom_l_aer   + rem_pom_s_aer   + rem_pom_r_aer
         rem_denit_total = rem_pom_l_denit + rem_pom_s_denit + rem_pom_r_denit
         rem_total       = rem_aer_total + rem_denit_total

         ! Total remineralisation rate for each POM pool.
         rem_pom_l_tot = rem_pom_l_aer + rem_pom_l_denit
         rem_pom_s_tot = rem_pom_s_aer + rem_pom_s_denit
         rem_pom_r_tot = rem_pom_r_aer + rem_pom_r_denit

         ! NH4 production during remineralisation.
         ! For each mol N of remineralised organic matter, 1 mol of NH4 is produced.
         ! Convert NH4 production from solid-phase to porewater-volume units.
         nh4_prod = p2d * rem_total

         ! PO4 production
         po4_prod_rem = p2d * (rem_pom_l_tot / self%n_to_p_pom_l  &
                         + rem_pom_s_tot / self%n_to_p_pom_s  &
                         + rem_pom_r_tot / self%n_to_p_pom_r)

         ! DIC production is computed from total remineralised organic matter using the C:N ratio of each POM pool.
         ! 1 mol of DIC is produced per mol of C remineralised in OM. 
         dic_prod_rem = self%c_to_n_pom_l * (rem_pom_l_aer  + rem_pom_l_denit) &
                      + self%c_to_n_pom_s * (rem_pom_s_aer  + rem_pom_s_denit) &
                      + self%c_to_n_pom_r * (rem_pom_r_aer  + rem_pom_r_denit)    

         ! Convert DIC production from solid-phase to porewater-volume units.           
         dic_prod_rem = p2d * dic_prod_rem        

         ! Alkalinity
         ! Converts PO4 production from solid-phase volume units to porewater concentration using p2d 
         alk_prod_tot = p2d * (alk_prod_aer + alk_prod_denit)
         

         !-------------------------------------------------------------------------
         !               Update tendencies and diagnostics
         !-------------------------------------------------------------------------
         ! POM budgets: production from pelagic biology, loss by remineralisation.
         _ADD_SOURCE_(self%id_pom_l, self%frac_lab  * pom_prod_n - rem_pom_l_tot)
         _ADD_SOURCE_(self%id_pom_s, self%frac_semi * pom_prod_n - rem_pom_s_tot)
         _ADD_SOURCE_(self%id_pom_r, self%frac_ref  * pom_prod_n - rem_pom_r_tot)

         ! Remineralisation products and oxidant consumption.
         _ADD_SOURCE_(self%id_nh4,  nh4_prod)
         _ADD_SOURCE_(self%id_no3, -no3_cons_denit)
         _ADD_SOURCE_(self%id_o2,  -o2_cons_rem)
         _ADD_SOURCE_(self%id_po4, po4_prod_rem) 

         if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic, dic_prod_rem)         
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_prod_tot)        

         _SET_DIAGNOSTIC_(self%id_rem, rem_total * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_total_pom, pom_l + pom_s + pom_r)
         _SET_DIAGNOSTIC_(self%id_o2_cons_rem, o2_cons_rem * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_n_loss_denit, n_loss_denit * secs_pr_day)

      _LOOP_END_
   end subroutine do

end module om_degradation