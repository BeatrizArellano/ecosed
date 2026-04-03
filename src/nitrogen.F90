#include "fabm_driver.h"

! ---------------------------------------------------------------------
! Nitrigen dissolved species and Nitrification
!
! Nitrification: oxidation of ammonium to nitrate
! Overall simplified reaction:
!   NH4+ + 2O2 -> NO3- + 2H+ + H2O
!
! Model representation:
!   - NH4 is consumed and NO3 is produced (1:1 in N units)
!   - O2 is consumed at a fixed ratio (2 mol O2 per mol NH4 oxidised)
!   - No direct DIC production (inorganic process)
!   - Produces acidity (2 H+), thus reduces alkalinity
!
! Rate formulation:
!   - First-order in NH4
!   - Enhanced by temperature (Eppley scaling)
!   - Suppressed under low O2 (smooth threshold)
!   - Inhibited by light (PAR-dependent)
! ---------------------------------------------------------------------
module nitrogen

   use fabm_types
   use molecular_diff, only: DIFF_ION_LINEAR, m0_NO3, m1_NO3, m0_NH4, m1_NH4, m0_PO4, m1_PO4
   implicit none
   private



   type, extends(type_base_model), public :: type_nitrogen
      ! --- State variables
      type(type_state_variable_id) :: id_no3
      type(type_state_variable_id) :: id_nh4
      type(type_state_variable_id) :: id_po4

      !--- Optional Coupling
      type(type_state_variable_id) :: id_o2      
      type(type_state_variable_id) :: id_alk

      ! --- Dependencies
      type(type_dependency_id)     :: id_temp
      type(type_dependency_id)     :: id_par

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_nit

      !--- Parameters
      real(rk) :: k_nit             ! First-order nitrification rate (s-1 internally)
      real(rk) :: k_par_nit         ! PAR where nitrification is reduced by 50%

   contains
      procedure :: initialize
      procedure :: do
   end type type_nitrogen

contains

   subroutine initialize(self, configunit)
      class(type_nitrogen), intent(inout), target :: self
      integer,              intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk

      call self%get_parameter(self%k_nit, 'k_nit', 'd-1', 'First-order nitrification rate',default=0.05_rk, scale_factor=d_per_s, minimum=0.0_rk)
      call self%get_parameter(self%k_par_nit, 'k_par_nit', 'W m-2', 'PAR at which nitrification is reduced by 50%', default=20.0_rk, minimum=0.0_rk)

      ! ---------------- State variables ----------------
      ! --- NO3 ---
      call self%register_state_variable(self%id_no3, 'no3', 'mmol m-3', 'Dissolved Nitrate', initial_value=5.0_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_no3, 'is_solute', .true.)
      call self%set_variable_property(self%id_no3, 'diff_method', DIFF_ION_LINEAR)
      call self%set_variable_property(self%id_no3, 'm0', m0_NO3)
      call self%set_variable_property(self%id_no3, 'm1', m1_NO3)

      ! --- NH4 ---
      call self%register_state_variable(self%id_nh4, 'nh4', 'mmol m-3', 'Dissolved Ammonium', initial_value=5.0_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_nh4, 'is_solute', .true.)
      call self%set_variable_property(self%id_nh4, 'diff_method', DIFF_ION_LINEAR)
      call self%set_variable_property(self%id_nh4, 'm0', m0_NH4)
      call self%set_variable_property(self%id_nh4, 'm1', m1_NH4)

      ! --- PO4 ---
      ! Phosphate is defined here as it varies with N depending on stoichiometry
      call self%register_state_variable(self%id_po4, 'po4', 'mmol m-3', 'Dissolved Phosphate', initial_value=0.3125_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_po4, 'is_solute', .true.)
      call self%set_variable_property(self%id_po4, 'diff_method', DIFF_ION_LINEAR)
      call self%set_variable_property(self%id_po4, 'm0', m0_PO4)
      call self%set_variable_property(self%id_po4, 'm1', m1_PO4)
      
      ! --- Contribution to total Nitrogen
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_no3)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_nh4)

      ! ---------------- Couplings ----------------
      call self%register_state_dependency(self%id_o2,  'o2',  'mmol m-3', 'Dissolved oxygen', required=.false.)
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.false.)     

      ! ---------------- Dependencies ----------------
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

      ! ---------------- Diagnostics ----------------
      call self%register_diagnostic_variable(self%id_nit, 'NIT', 'mmol m-3 d-1', 'Nitrification rate')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_nitrogen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: nh4, temp, par, o2
      real(rk) :: fT, fI, fO2
      real(rk) :: nit, o2_cons_nit, alk_change

      real(rk), parameter :: secs_pr_day  = 86400.0_rk
      real(rk), parameter :: o2_nit_thr   = 2.0_rk     ! O2 Threshold to reduce nitrification rapidly
      real(rk), parameter :: o2_nit_width = 0.2_rk
      real(rk), parameter :: o2_per_n_nit = 2.0_rk     ! 2 mol O2 per mol of NH4

      _LOOP_BEGIN_

         _GET_(self%id_nh4, nh4)
         _GET_(self%id_temp, temp)         
         _GET_(self%id_par, par)         

         !------------------------------------------------------------------------
         !   Nitrification (NH4 -> NO3)
         !   NH4+ + 2 O2 -> NO3- + 2H+ + H2O
         ! Consumes NH4 and O2, produces NO3 (no DIC change, alkalinity decreases).
         ! First-order in NH4, enhanced by temperature, reduced by low O2 and light.
         !------------------------------------------------------------------------
         
         ! Smooth O2 limitation: nitrification declines rapidly below the threshold.
         if (_AVAILABLE_(self%id_o2)) then
            _GET_(self%id_o2, o2)
            fO2 = 0.5_rk * (1.0_rk + tanh((max(o2,0.0_rk) - o2_nit_thr) / o2_nit_width))
         else
            fO2 = 1.0_rk
         end if

         ! Eppley-type temperature scaling
         fT = 1.066_rk ** temp

         ! Light inhibits nitrification; k_par_nit sets the 50% inhibition level.
         fI = 1.0_rk / (1.0_rk + max(par, 0.0_rk) / self%k_par_nit)

         ! Nitrification rate (mmol N m-3 s-1).
         nit = self%k_nit * fT * fO2 * fI * max(nh4, 0.0_rk)

         ! Stoichiometric O2 consumption.
         o2_cons_nit = o2_per_n_nit * nit

         ! Alkalinity consumption by nitrification. No change in DIC here. 
         ! Eq. R7 in Middelburg et al. (2020): NH4+ + 2 O2 -> NO3- + 2H+ + H2O
         ! 2 equivalents of alkalinity are consumed per mol NH4 oxidised.
         alk_change = -2.0_rk * nit

         ! Update tracer tendencies
         _ADD_SOURCE_(self%id_nh4, -nit)
         _ADD_SOURCE_(self%id_no3,  nit)

         if (_AVAILABLE_(self%id_o2)) _ADD_SOURCE_(self%id_o2, -o2_cons_nit)
         if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk, alk_change)

         ! Diagnostic variables
         _SET_DIAGNOSTIC_(self%id_nit, nit * secs_pr_day)

      _LOOP_END_

   end subroutine do


end module nitrogen