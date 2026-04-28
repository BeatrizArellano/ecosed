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
      type(type_state_variable_id) :: id_o2
      type(type_state_variable_id) :: id_no3, id_nh4
      type(type_state_variable_id) :: id_po4      
      type(type_state_variable_id) :: id_mn2, id_mno2
      type(type_state_variable_id) :: id_fe2, id_fe3ox
      type(type_state_variable_id) :: id_so4, id_sulfide

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
      type(type_diagnostic_variable_id) :: id_rem_aer
      type(type_diagnostic_variable_id) :: id_rem_denit
      type(type_diagnostic_variable_id) :: id_rem_mn
      type(type_diagnostic_variable_id) :: id_mno2_cons_mn
      type(type_diagnostic_variable_id) :: id_mn2_prod_mn
      type(type_diagnostic_variable_id) :: id_alk_prod_aer
      type(type_diagnostic_variable_id) :: id_alk_prod_denit
      type(type_diagnostic_variable_id) :: id_alk_prod_mn
      type(type_diagnostic_variable_id) :: id_rem_fe
      type(type_diagnostic_variable_id) :: id_fe3ox_cons_fe
      type(type_diagnostic_variable_id) :: id_fe2_prod_fe
      type(type_diagnostic_variable_id) :: id_alk_prod_fe
      type(type_diagnostic_variable_id) :: id_rem_so4
      type(type_diagnostic_variable_id) :: id_so4_cons
      type(type_diagnostic_variable_id) :: id_sulfide_prod
      type(type_diagnostic_variable_id) :: id_alk_prod_so4

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
      real(rk) :: ki_no3_redox
      real(rk) :: k_mno2_red
      real(rk) :: ki_mno2_redox
      real(rk) :: k_fe3ox_red
      real(rk) :: ki_fe3ox_redox
      real(rk) :: k_so4_red
      real(rk) :: ki_so4_redox   
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
      call self%get_parameter(self%k_mno2_red, 'k_mno2_red', 'mmol m-3', 'Half-saturation constant for MnO2 limitation in Mn oxide reduction', default=42.4_rk, minimum=eps)
      call self%get_parameter(self%k_fe3ox_red, 'k_fe3ox_red', 'mmol m-3', 'Half-saturation constant for Fe(III) oxide limitation in Fe oxide reduction', default=100.0_rk, minimum=eps)
      call self%get_parameter(self%k_so4_red, 'k_so4_red', 'mmol m-3', 'Half-saturation constant for sulfate limitation in sulfate reduction', default=1000.0_rk, minimum=eps)

      call self%get_parameter(self%ki_no3_redox, 'ki_no3_redox', 'mmol m-3', 'NO3 inhibition constant for lower-energy redox pathways', default=5.0_rk, minimum=eps)
      call self%get_parameter(self%ki_mno2_redox, 'ki_mno2_redox', 'mmol m-3', 'MnO2 inhibition constant for lower-energy redox pathways', default=42.4_rk, minimum=eps)
      call self%get_parameter(self%ki_fe3ox_redox, 'ki_fe3ox_redox', 'mmol m-3', 'Fe(III) oxide inhibition constant for lower-energy redox pathways', default=100.0_rk, minimum=eps)
      call self%get_parameter(self%ki_so4_redox, 'ki_so4_redox', 'mmol m-3', 'SO4 inhibition constant for methanogenesis or lower-energy pathways', default=1000.0_rk, minimum=eps)

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
      call self%register_state_dependency(self%id_mn2,  'mn2',  'mmol m-3', 'Dissolved Mn(II)', required=.false.)
      call self%register_state_dependency(self%id_mno2, 'mno2', 'mmol m-3', 'Particulate Mn(IV) oxide', required=.false.)
      call self%register_state_dependency(self%id_fe2, 'fe2', 'mmol m-3', 'Dissolved Fe(II)', required=.false.)
      call self%register_state_dependency(self%id_fe3ox, 'fe3ox', 'mmol m-3', 'Particulate Fe(III) oxides/hydroxides', required=.false.)
      call self%register_state_dependency(self%id_so4, 'so4', 'mmol m-3', 'Dissolved sulfate', required=.false.)
      call self%register_state_dependency(self%id_sulfide, 'sulfide', 'mmol m-3', 'Total dissolved sulfide', required=.false.)

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
      call self%register_diagnostic_variable(self%id_rem_aer,   'REM_AER',   'mmol N m-3 d-1', 'Aerobic remineralisation rate')
      call self%register_diagnostic_variable(self%id_rem_denit, 'REM_DENIT', 'mmol N m-3 d-1', 'Denitrification remineralisation rate')
      call self%register_diagnostic_variable(self%id_rem_mn,    'REM_MN',    'mmol N m-3 d-1', 'Mn oxide reduction remineralisation rate')

      call self%register_diagnostic_variable(self%id_mno2_cons_mn, 'MNO2_CONS_MN', 'mmol m-3 d-1', 'MnO2 consumption by Mn oxide reduction')
      call self%register_diagnostic_variable(self%id_mn2_prod_mn,  'MN2_PROD_MN',  'mmol m-3 d-1', 'Mn2 production by Mn oxide reduction')

      call self%register_diagnostic_variable(self%id_rem_fe, 'REM_FE', 'mmol N m-3 d-1', 'Fe oxide reduction remineralisation rate')
      call self%register_diagnostic_variable(self%id_fe3ox_cons_fe, 'FE3OX_CONS_FE', 'mmol m-3 d-1', 'Fe(III) oxide consumption by Fe oxide reduction')
      call self%register_diagnostic_variable(self%id_fe2_prod_fe, 'FE2_PROD_FE', 'mmol m-3 d-1', 'Fe2 production by Fe oxide reduction')
      call self%register_diagnostic_variable(self%id_alk_prod_fe, 'ALK_PROD_FE', 'mmol eq m-3 d-1', 'Alkalinity production by Fe oxide reduction')

      call self%register_diagnostic_variable(self%id_rem_so4, 'REM_SO4', 'mmol N m-3 d-1', 'Sulfate reduction remineralisation rate')
      call self%register_diagnostic_variable(self%id_so4_cons, 'SO4_CONS', 'mmol m-3 d-1', 'Sulfate consumption')
      call self%register_diagnostic_variable(self%id_sulfide_prod, 'SULFIDE_PROD', 'mmol m-3 d-1', 'Sulfide production')
      call self%register_diagnostic_variable(self%id_alk_prod_so4, 'ALK_PROD_SO4', 'mmol eq m-3 d-1', 'Alkalinity production by sulfate reduction')

      call self%register_diagnostic_variable(self%id_alk_prod_aer,   'ALK_PROD_AER',   'mmol eq m-3 d-1', 'Alkalinity production by aerobic remineralisation')
      call self%register_diagnostic_variable(self%id_alk_prod_denit, 'ALK_PROD_DENIT', 'mmol eq m-3 d-1', 'Alkalinity production by denitrification')
      call self%register_diagnostic_variable(self%id_alk_prod_mn,    'ALK_PROD_MN',    'mmol eq m-3 d-1', 'Alkalinity production by Mn oxide reduction')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_om_degradation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: pom_l, pom_s, pom_r
      real(rk) :: temp, pom_prod_n
      real(rk) :: phi, phi_s                     ! Porosity
      real(rk) :: p2d, d2p
      real(rk) :: o2, no3, fT
      real(rk) :: mn2, mno2
      real(rk) :: fe2, fe3ox
      real(rk) :: so4, sulfide

      logical  :: is_water
      ! Pathway limiting/inhibition functions
      real(rk) :: fi_o2_redox
      real(rk) :: fi_no3_redox
      real(rk) :: fi_mno2_redox, fi_fe3ox_redox, fi_so4_redox
      real(rk) :: fsum, faer, fdenit, fmn, ffe, fso4

      ! Pathway-specific remineralisation rates
      real(rk) :: rem_pom_l_aer, rem_pom_s_aer, rem_pom_r_aer
      real(rk) :: rem_pom_l_denit, rem_pom_s_denit, rem_pom_r_denit
      real(rk) :: rem_pom_l_mn, rem_pom_s_mn, rem_pom_r_mn
      real(rk) :: rem_pom_l_fe, rem_pom_s_fe, rem_pom_r_fe
      real(rk) :: rem_pom_l_so4, rem_pom_s_so4, rem_pom_r_so4
      real(rk) :: rem_pom_l_tot, rem_pom_s_tot, rem_pom_r_tot
      real(rk) :: rem_aer_total, rem_denit_total, rem_mn_total, &
                  rem_fe_total, rem_so4_total, rem_total

      ! Stoichiometric products and sinks
      real(rk) :: dic_prod_rem
      real(rk) :: o2_cons_rem
      real(rk) :: no3_cons_denit
      real(rk) :: n_loss_denit
      real(rk) :: nh4_prod
      real(rk) :: po4_prod_rem
      real(rk) :: c_remin_mn, mno2_cons_mn, mn2_prod_mn
      real(rk) :: c_remin_fe, fe3ox_cons_fe, fe2_prod_fe
      real(rk) :: c_remin_so4, so4_cons, sulfide_prod
      real(rk) :: alk_prod_aer, alk_prod_denit, alk_prod_mn, alk_prod_fe, alk_prod_so4, alk_prod_tot

      logical  :: use_mn, use_fe, use_so4

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

         use_mn = _AVAILABLE_(self%id_mn2) .and. _AVAILABLE_(self%id_mno2)
         if (use_mn) then
            _GET_(self%id_mn2,  mn2)
            _GET_(self%id_mno2, mno2)
         else
            mn2  = 0.0_rk
            mno2 = 0.0_rk
         end if

         use_fe = _AVAILABLE_(self%id_fe2) .and. _AVAILABLE_(self%id_fe3ox)
         if (use_fe) then
            _GET_(self%id_fe2,    fe2)
            _GET_(self%id_fe3ox,  fe3ox)
         else
            fe2   = 0.0_rk
            fe3ox = 0.0_rk
         end if

         use_so4 = _AVAILABLE_(self%id_so4) .and. _AVAILABLE_(self%id_sulfide)
         if (use_so4) then
            _GET_(self%id_so4, so4)
            _GET_(self%id_sulfide, sulfide)
         else
            so4 = 0.0_rk
            sulfide = 0.0_rk
         end if

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

         ! -------------------------------------------------------------------------
         ! Sequential redox inhibition factors.
         !
         ! The redox ladder is represented by allowing lower-energy pathways only
         ! when higher-energy electron acceptors are depleted.
         !
         ! Sequence:
         !   O2 -> NO3 -> MnO2 -> Fe(III) oxides -> SO4 -> CH4
         ! -------------------------------------------------------------------------

         ! Inhibition by nitrate: suppress Mn, Fe and sulfate reduction while NO3 is available.
         fi_no3_redox = self%ki_no3_redox / (self%ki_no3_redox + no3)

         ! Inhibition by MnO2: suppress Fe and sulfate reduction while MnO2 is available.
         if (use_mn) then
            fi_mno2_redox = self%ki_mno2_redox / (self%ki_mno2_redox + mno2)
         else
            fi_mno2_redox = 1.0_rk
         end if

         ! Inhibition by Fe(III) oxides: suppress sulfate reduction while Fe oxides are available.
         if (use_fe) then
            fi_fe3ox_redox = self%ki_fe3ox_redox / (self%ki_fe3ox_redox + fe3ox)
         else
            fi_fe3ox_redox = 1.0_rk
         end if

         ! Optional sulfate inhibition, useful later if you add methanogenesis.
         if (use_so4) then
            fi_so4_redox = self%ki_so4_redox / (self%ki_so4_redox + so4)
         else
            fi_so4_redox = 1.0_rk
         end if

         ! NO3 limitation for denitrification and inhibition by O2.
         ! Denitrification requires both low-O2 conditions and NO3 availability.
         fdenit = fi_o2_redox * no3 / (self%k_no3_denit + no3)

         ! MnO2 reduction is suppressed by O2 and NO3.
         if (use_mn) then
            fmn = fi_o2_redox * fi_no3_redox * mno2 / (self%k_mno2_red + mno2)
         else
            fmn = 0.0_rk
         end if

         ! Fe(III) oxide reduction is suppressed by O2, NO3 and MnO2.
         if (use_fe) then
            ffe = fi_o2_redox * fi_no3_redox * fi_mno2_redox * fe3ox / (self%k_fe3ox_red + fe3ox)
         else
            ffe = 0.0_rk
         end if

         ! Sulfate reduction is suppressed by O2, NO3, MnO2 and Fe(III) oxides.
         if (use_so4) then
            fso4 = fi_o2_redox * fi_no3_redox * fi_mno2_redox * fi_fe3ox_redox &
                 * so4 / (self%k_so4_red + so4)
         else
            fso4 = 0.0_rk
         end if


         ! Combined activity (demand) of all remineralisation pathways under current redox conditions.
         fsum = faer + fdenit + fmn + ffe + fso4

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
            fmn    = fmn    / fsum
            ffe    = ffe    / fsum
            fso4   = fso4   / fsum
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
         !    Mn-Oxide Reduction
         ! ------------------------------------------------------------------------
         ! Organic matter degradation coupled to Mn(IV) oxide reduction.
         ! Equation R3 in Middelburg et al. (2020):
         ! (CH2O)(NH3)n/c(H3PO4)p/c + 2 MnO2 + 4H+ 
         !    -> CO2 + n/c NH3 + p/c H3PO4 + 2 Mn2+ + 3H2O
         !
         ! Per mol organic C remineralised:
         !   MnO2 consumption = 2 per mol C
         !   Mn2 production   = 2 per mol C
         !   DIC production   = 1 per mol C
         !   TA production    = 4 + n/c - p/c per mol C
         rem_pom_l_mn = 0.0_rk; rem_pom_s_mn = 0.0_rk; rem_pom_r_mn = 0.0_rk
         rem_mn_total = 0.0_rk; mno2_cons_mn = 0.0_rk; mn2_prod_mn  = 0.0_rk; alk_prod_mn  = 0.0_rk

         ! MnO2 reduction rate for each POM pool, limited by MnO2 and inhibited by O2 and NO3.
         rem_pom_l_mn = self%k_remin_pom_l * fmn * pom_l                ! mmol N m⁻3 s⁻1
         rem_pom_s_mn = self%k_remin_pom_s * fmn * pom_s
         rem_pom_r_mn = self%k_remin_pom_r * fmn * pom_r

         rem_mn_total = rem_pom_l_mn + rem_pom_s_mn + rem_pom_r_mn      ! Total degradation rate via MnO2 reduction

         ! POM is stored in N units, so C remineralisation is C:N * rem_pom_*_mn. (Converting to mmol C m⁻3 s-1)
         c_remin_mn = self%c_to_n_pom_l * rem_pom_l_mn &
                    + self%c_to_n_pom_s * rem_pom_s_mn &
                    + self%c_to_n_pom_r * rem_pom_r_mn

         ! MnO2 is particulate and stored on the solid fraction, like POM.
         ! Therefore this sink remains in solid-volume units.
         mno2_cons_mn = 2.0_rk * c_remin_mn

         ! Mn2 is dissolved and enters porewater, so convert with p2d.
         mn2_prod_mn = p2d * 2.0_rk * c_remin_mn

         ! Alkalinity generation by MnO2 reduction.
         ! Per mol N remineralised:
         !   4*C:N + 1 - 1/(N:P)
         alk_prod_mn = rem_pom_l_mn * (4.0_rk*self%c_to_n_pom_l + 1.0_rk - 1.0_rk/self%n_to_p_pom_l) &
                     + rem_pom_s_mn * (4.0_rk*self%c_to_n_pom_s + 1.0_rk - 1.0_rk/self%n_to_p_pom_s) &
                     + rem_pom_r_mn * (4.0_rk*self%c_to_n_pom_r + 1.0_rk - 1.0_rk/self%n_to_p_pom_r)  
                     
         ! ------------------------------------------------------------------------
         !    Fe-Oxide Reduction
         ! ------------------------------------------------------------------------
         ! Organic matter degradation via Fe(III) oxide reduction.
         !
         ! Written using Fe(OH)3 as the form of the reactive Fe(III)
         ! oxide/hydroxide pool:
         !
         ! (CH2O)(NH3)n/c(H3PO4)p/c + 4 Fe(OH)3 + 7 CO2
         !    -> 8 HCO3- + n/c NH3 + p/c H3PO4 + 4 Fe2+ + 3 H2O
         !
         ! Per mol organic C remineralised:
         !   Fe3ox consumption = 4 per mol C
         !   Fe2 production    = 4 per mol C
         !   DIC production    = 1 per mol C
         !   TA production     = 8 + n/c - p/c per mol C
         rem_pom_l_fe = 0.0_rk; rem_pom_s_fe = 0.0_rk; rem_pom_r_fe = 0.0_rk
         rem_fe_total = 0.0_rk; fe3ox_cons_fe = 0.0_rk; fe2_prod_fe = 0.0_rk
         alk_prod_fe  = 0.0_rk

         rem_pom_l_fe = self%k_remin_pom_l * ffe * pom_l
         rem_pom_s_fe = self%k_remin_pom_s * ffe * pom_s
         rem_pom_r_fe = self%k_remin_pom_r * ffe * pom_r

         rem_fe_total = rem_pom_l_fe + rem_pom_s_fe + rem_pom_r_fe

         ! POM is stored in N units, so C remineralisation is C:N * rem_pom_*_fe.
         c_remin_fe = self%c_to_n_pom_l * rem_pom_l_fe &
                    + self%c_to_n_pom_s * rem_pom_s_fe &
                    + self%c_to_n_pom_r * rem_pom_r_fe

         ! Fe3ox is particulate and stored on the solid fraction, like POM.
         fe3ox_cons_fe = 4.0_rk * c_remin_fe

         ! Fe2 is dissolved and enters porewater, so convert with p2d.
         fe2_prod_fe = p2d * 4.0_rk * c_remin_fe

         ! Alkalinity generation by Fe(III) oxide reduction.
         ! Per mol N remineralised:
         !   8*C:N + 1 - 1/(N:P)
         alk_prod_fe = rem_pom_l_fe * (8.0_rk*self%c_to_n_pom_l + 1.0_rk - 1.0_rk/self%n_to_p_pom_l) &
                     + rem_pom_s_fe * (8.0_rk*self%c_to_n_pom_s + 1.0_rk - 1.0_rk/self%n_to_p_pom_s) &
                     + rem_pom_r_fe * (8.0_rk*self%c_to_n_pom_r + 1.0_rk - 1.0_rk/self%n_to_p_pom_r)


         ! ------------------------------------------------------------------------
         !    Sulfate Reduction
         ! ------------------------------------------------------------------------
         ! (CH2O)(NH3)n/c(H3PO4)p/c + 0.5 SO4--
         !    -> CO2 + n/c NH3 + p/c H3PO4 + 0.5 H2S + H2O
         !
         ! Per mol C remineralised:
         !   SO4 consumption = 0.5 per mol C
         !   H2S production  = 0.5 per mol C
         !   DIC production  = 1 per mol C
         !   TA production   = 1 + n/c - p/c  per mol C
         !
         rem_pom_l_so4 = 0.0_rk; rem_pom_s_so4 = 0.0_rk; rem_pom_r_so4 = 0.0_rk
         rem_so4_total = 0.0_rk; so4_cons      = 0.0_rk; sulfide_prod  = 0.0_rk
         alk_prod_so4  = 0.0_rk

         if (use_so4) then

            rem_pom_l_so4 = self%k_remin_pom_l * fso4 * pom_l
            rem_pom_s_so4 = self%k_remin_pom_s * fso4 * pom_s
            rem_pom_r_so4 = self%k_remin_pom_r * fso4 * pom_r

            rem_so4_total = rem_pom_l_so4 + rem_pom_s_so4 + rem_pom_r_so4

            c_remin_so4 = self%c_to_n_pom_l * rem_pom_l_so4 &
                        + self%c_to_n_pom_s * rem_pom_s_so4 &
                        + self%c_to_n_pom_r * rem_pom_r_so4

             ! Convert SO4 consumption from solid-phase to porewater-volume units,
            so4_cons = p2d * 0.5_rk * c_remin_so4

            ! sulfide is dissolved -> convert with p2d
            sulfide_prod = p2d * 0.5_rk * c_remin_so4

            ! Alkalinity (per mol N):
            ! C:N + 1 - 1/(N:P)
            alk_prod_so4 = rem_pom_l_so4 * (self%c_to_n_pom_l + 1.0_rk - 1.0_rk/self%n_to_p_pom_l) &
                         + rem_pom_s_so4 * (self%c_to_n_pom_s + 1.0_rk - 1.0_rk/self%n_to_p_pom_s) &
                         + rem_pom_r_so4 * (self%c_to_n_pom_r + 1.0_rk - 1.0_rk/self%n_to_p_pom_r)
         end if

         
         ! ------------------------------------------------------------------------
         !    Total remineralisation rates and products
         ! ------------------------------------------------------------------------
         ! Total remineralisation rates summed across POM pools and pathways.
         rem_aer_total   = rem_pom_l_aer   + rem_pom_s_aer   + rem_pom_r_aer
         rem_denit_total = rem_pom_l_denit + rem_pom_s_denit + rem_pom_r_denit
         rem_total = rem_aer_total + rem_denit_total + rem_mn_total + rem_fe_total + rem_so4_total
         
         ! Total remineralisation rate for each POM pool.
         rem_pom_l_tot = rem_pom_l_aer + rem_pom_l_denit + rem_pom_l_mn + rem_pom_l_fe + rem_pom_l_so4
         rem_pom_s_tot = rem_pom_s_aer + rem_pom_s_denit + rem_pom_s_mn + rem_pom_s_fe + rem_pom_s_so4
         rem_pom_r_tot = rem_pom_r_aer + rem_pom_r_denit + rem_pom_r_mn + rem_pom_r_fe + rem_pom_r_so4

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
         dic_prod_rem = self%c_to_n_pom_l * rem_pom_l_tot &
                      + self%c_to_n_pom_s * rem_pom_s_tot &
                      + self%c_to_n_pom_r * rem_pom_r_tot   

         ! Convert DIC production from solid-phase to porewater-volume units.           
         dic_prod_rem = p2d * dic_prod_rem        

         ! Alkalinity
         ! Convert alkalinity production from solid-phase volume units to porewater-volume units.
         alk_prod_tot = p2d * (alk_prod_aer + alk_prod_denit + alk_prod_mn + alk_prod_fe + alk_prod_so4)
         

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
         if (use_mn) then
            _ADD_SOURCE_(self%id_mno2, -mno2_cons_mn)
            _ADD_SOURCE_(self%id_mn2,   mn2_prod_mn)
         end if
         if (use_fe) then
            _ADD_SOURCE_(self%id_fe3ox, -fe3ox_cons_fe)
            _ADD_SOURCE_(self%id_fe2,    fe2_prod_fe)
         end if
         if (use_so4) then
            _ADD_SOURCE_(self%id_so4,     -so4_cons)
            _ADD_SOURCE_(self%id_sulfide,  sulfide_prod)
         end if

         if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic, dic_prod_rem)         
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_prod_tot)        

         _SET_DIAGNOSTIC_(self%id_rem, rem_total * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_total_pom, pom_l + pom_s + pom_r)
         _SET_DIAGNOSTIC_(self%id_o2_cons_rem, o2_cons_rem * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_n_loss_denit, n_loss_denit * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_rem_aer,   rem_aer_total   * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_rem_denit, rem_denit_total * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_rem_mn,    rem_mn_total    * secs_pr_day)

         _SET_DIAGNOSTIC_(self%id_mno2_cons_mn, mno2_cons_mn * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_mn2_prod_mn,  mn2_prod_mn  * secs_pr_day)

         _SET_DIAGNOSTIC_(self%id_rem_fe, rem_fe_total * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_fe3ox_cons_fe, fe3ox_cons_fe * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_fe2_prod_fe,   fe2_prod_fe   * secs_pr_day)

         _SET_DIAGNOSTIC_(self%id_rem_so4, rem_so4_total * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_so4_cons, so4_cons * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_sulfide_prod, sulfide_prod * secs_pr_day)         

         _SET_DIAGNOSTIC_(self%id_alk_prod_aer,   p2d * alk_prod_aer   * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_alk_prod_denit, p2d * alk_prod_denit * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_alk_prod_mn,    p2d * alk_prod_mn    * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_alk_prod_fe,    p2d * alk_prod_fe    * secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_alk_prod_so4,   p2d * alk_prod_so4 * secs_pr_day)

      _LOOP_END_
   end subroutine do

end module om_degradation