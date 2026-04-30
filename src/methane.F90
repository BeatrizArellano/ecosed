#include "fabm_driver.h"

! ---------------------------------------------------------------------
! Methane redox cycle
!
! Minimal representation:
!   - ch4 : dissolved methane, CH4
!
! Aerobic methane oxidation:
!   CH4 + 2 O2 -> CO2 + 2 H2O
! Sulfate anaerobic oxidation of methane:
!   CH4 + SO4-- -> HCO3- + HS- + H2O
!
! Model representation:
!   - Methane is consumed
!   - O2 is consumed at 2 mol O2 per mol CH4 oxidised
!   - DIC is produced at 1 mol C per mol CH4 oxidised
!   - 2 mol alkalinity are produced per mol CH4 oxidised through sulfate anaerobic oxidation
!
! Methanogenesis from organic matter is handled in the OM degradation
! module, because it is an organic matter remineralisation pathway.
! ---------------------------------------------------------------------
module methane

    use fabm_types
    use molecular_diff, only: DIFF_ARRHENIUS, Ea_CH4, A0_CH4
    implicit none
    private 
    type, extends(type_base_model), public :: type_methane

        ! --- State variables
        type(type_state_variable_id) :: id_ch4      
        ! --- Optional/required couplings
        type(type_state_variable_id) :: id_o2
        type(type_state_variable_id) :: id_so4
        type(type_state_variable_id) :: id_sulfide
        type(type_state_variable_id) :: id_dic
        type(type_state_variable_id) :: id_alk      
        ! --- Diagnostics
        type(type_diagnostic_variable_id) :: id_ch4_ox      
        type(type_diagnostic_variable_id) :: id_ch4_an
        ! --- Parameters
        real(rk) :: k_ch4_ox      ! Methane oxidation rate (s-1 internally)
        real(rk) :: k_ch4_so4     ! Sulfate-dependent anaerobic CH4 oxidation rate (s-1 internally)
        real(rk) :: k_so4_ox      ! SO4 half-saturation for AOM (mmol m-3)
        
    contains
        procedure :: initialize
        procedure :: do
    end type type_methane

contains

    subroutine initialize(self, configunit)
        class(type_methane), intent(inout), target :: self
        integer,             intent(in)            :: configunit 

        real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
        real(rk), parameter :: eps     = 1.0e-12_rk 

        ! ---------------- Parameters ----------------
        call self%get_parameter(self%k_ch4_ox, 'k_ch4_ox', 'd-1', 'First-order aerobic methane oxidation rate', &
                                default=0.1_rk, scale_factor=d_per_s, minimum=0.0_rk)
        call self%get_parameter(self%k_ch4_so4, 'k_ch4_so4', 'd-1', 'First-order sulfate-driven methane oxidation rate', &
                                default=0.05_rk, scale_factor=d_per_s, minimum=0.0_rk)  
        call self%get_parameter(self%k_so4_ox, 'k_so4_ox', 'mmol m-3', 'Sulfate half-saturation concentration for anaerobic methane oxidation', &
                                default=1600.0_rk, minimum=eps)
 
        ! ---------------- State variables ----------------
        call self%register_state_variable(self%id_ch4, 'ch4', 'mmol m-3', 'Dissolved methane', initial_value=0.0_rk, minimum=0.0_rk, no_river_dilution=.true.)
        call self%set_variable_property(self%id_ch4, 'is_solute', .true.)
        call self%set_variable_property(self%id_ch4, 'diff_method', DIFF_ARRHENIUS)
        call self%set_variable_property(self%id_ch4, 'A0', A0_CH4)
        call self%set_variable_property(self%id_ch4, 'Ea', Ea_CH4)   
        ! ---------------- Couplings ----------------
        call self%register_state_dependency(self%id_o2,  'o2',  'mmol m-3',   'Dissolved oxygen', required=.true.)
        call self%register_state_dependency(self%id_so4,     'so4',     'mmol m-3', 'Dissolved sulfate', required=.false.)
        call self%register_state_dependency(self%id_sulfide, 'sulfide', 'mmol m-3', 'Total dissolved sulfide', required=.false.)
        call self%register_state_dependency(self%id_dic, 'dic', 'mmol C m-3', 'Dissolved inorganic carbon', required=.false.)
        call self%register_state_dependency(self%id_alk, 'alk', 'mmol eq m-3', 'Total alkalinity', required=.false.) 
        ! ---------------- Diagnostics ----------------
        call self%register_diagnostic_variable(self%id_ch4_ox, 'ch4_ox', 'mmol m-3 d-1', 'Aerobic methane oxidation rate')   
        call self%register_diagnostic_variable(self%id_ch4_an, 'ch4_an', 'mmol m-3 d-1', 'Anaerobic methane oxidation rate via SO4')

    end subroutine initialize 

    subroutine do(self, _ARGUMENTS_DO_)
        class(type_methane), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        real(rk) :: ch4, o2, so4
        real(rk) :: f_oxic, f_anoxic, f_so4
        real(rk) :: ch4_ox, ch4_an

        real(rk), parameter :: secs_per_day    = 86400.0_rk
        real(rk), parameter :: o2_ch4_ox_thr   = 1.5_rk
        real(rk), parameter :: o2_ch4_ox_width = 0.2_rk

        real(rk), parameter :: o2_per_ch4_ox  = 2.0_rk
        real(rk), parameter :: so4_per_ch4_an = 1.0_rk
        real(rk), parameter :: dic_per_ch4_ox = 1.0_rk
        real(rk), parameter :: dic_per_ch4_an = 1.0_rk
        real(rk), parameter :: alk_per_ch4_an = 2.0_rk

        _LOOP_BEGIN_

            _GET_(self%id_ch4, ch4)
            _GET_(self%id_o2,  o2)

            ! ---------------------------------------------------------------
            ! Redox control of methane oxidation
            !
            ! Aerobic methane oxidation is active under oxic conditions:
            !   CH4 + 2 O2 -> CO2 + 2 H2O
            !
            ! Sulfate-dependent anaerobic oxidation of methane (AOM) is
            ! carried out by anaerobic microbial consortia and is therefore
            ! restricted here to low-O2 conditions:
            !   CH4 + SO4-- -> HCO3- + HS- + H2O
            !
            ! The smooth tanh switch avoids a discontinuous transition around
            ! the O2 threshold. f_oxic approaches 1 under oxic conditions and
            ! 0 under anoxic conditions; f_anoxic is its complement.
            ! ---------------------------------------------------------------
            f_oxic = 0.5_rk * (1.0_rk + tanh((max(o2, 0.0_rk) - o2_ch4_ox_thr) / o2_ch4_ox_width))
            f_anoxic = 1.0_rk - f_oxic

            ! Aerobic methane oxidation: first-order in CH4, enabled by O2.
            ch4_ox = self%k_ch4_ox * f_oxic * max(ch4, 0.0_rk)

            ! Sulfate-AOM: first-order in CH4, enabled under anoxia and
            ! limited by sulfate availability.
            if (_AVAILABLE_(self%id_so4) .and. _AVAILABLE_(self%id_sulfide)) then
                _GET_(self%id_so4, so4)
                f_so4 = max(so4, 0.0_rk) / (self%k_so4_ox + max(so4, 0.0_rk))
                ! Anaerobic oxidation rate of CH4
                ch4_an = self%k_ch4_so4 * f_anoxic * f_so4 * max(ch4, 0.0_rk)
            else
                ch4_an = 0.0_rk
            end if

            ! ---------------------------------------------------------------
            ! Stoichiometry
            ! Aerobic oxidation:
            !   - consumes 1 CH4
            !   - consumes 2 O2
            !   - produces 1 DIC
            !   - no direct alkalinity change
            !
            ! Sulfate-AOM:
            !   - consumes 1 CH4
            !   - consumes 1 SO4
            !   - produces 1 sulfide
            !   - produces 1 DIC
            !   - increases alkalinity by 2 equivalents
            ! ---------------------------------------------------------------

            _ADD_SOURCE_(self%id_ch4, -ch4_ox - ch4_an)
            _ADD_SOURCE_(self%id_o2,  -o2_per_ch4_ox * ch4_ox)

            if (_AVAILABLE_(self%id_so4))     _ADD_SOURCE_(self%id_so4,     -so4_per_ch4_an * ch4_an)
            if (_AVAILABLE_(self%id_sulfide)) _ADD_SOURCE_(self%id_sulfide,  ch4_an)
            if (_AVAILABLE_(self%id_dic))     _ADD_SOURCE_(self%id_dic,      dic_per_ch4_ox * ch4_ox + dic_per_ch4_an * ch4_an)
            if (_AVAILABLE_(self%id_alk))     _ADD_SOURCE_(self%id_alk,      alk_per_ch4_an * ch4_an)

            _SET_DIAGNOSTIC_(self%id_ch4_ox, ch4_ox * secs_per_day)
            _SET_DIAGNOSTIC_(self%id_ch4_an, ch4_an * secs_per_day)

        _LOOP_END_
    end subroutine do

end module methane