#include "fabm_driver.h"

!---------------------------------------------------------------------------------------------------
! Carbonate chemistry module
!
! This module represents seawater carbonate chemistry with prognostic dissolved inorganic carbon
! (DIC) and total alkalinity (ALK). From the local tracer and environmental state, it diagnoses pH,
! hydrogen ion concentration, carbonate speciation (CO2*, HCO3-, CO3--), and calcite saturation
! state. The carbonate system is solved from DIC and alkalinity using routines from
! carbonate_system, including salinity-derived background seawater totals, pressure-corrected
! equilibrium constants, and a safeguarded Newton-Raphson solver for [H+].
!
! Optional couplings to phosphate, silicate, ammonium and sulfide allow their contributions to
! alkalinity to be included in the speciation calculation. At the surface, the module also computes
! air-sea CO2 exchange from diagnosed aqueous CO2, atmospheric CO2, wind speed, pressure and
! optional sea-ice cover, and applies the resulting flux to DIC.
!
! Prognostic state variables:
!   - DIC : dissolved inorganic carbon      [mmol C  m-3]
!   - ALK : total alkalinity                [mmol eq m-3]
!---------------------------------------------------------------------------------------------------
module carbonate_chemistry

   use fabm_types
   use molecular_diff,   only: DIFF_ION_LINEAR, m0_HCO3, m1_HCO3
   use carbonate_system, only: compute_background_totals_from_salinity, compute_equilibrium_constants, &
                               solve_hplus, eps_h, calcon_ref
   implicit none
   private   

   type, extends(type_base_model), public :: type_carbonate_chemistry

      ! --- State variables
      type(type_state_variable_id) :: id_dic
      type(type_state_variable_id) :: id_alk

      ! --- Optional tracer couplings
      type(type_state_variable_id) :: id_po4
      type(type_state_variable_id) :: id_si
      type(type_state_variable_id) :: id_nh4
      type(type_state_variable_id) :: id_h2s

      ! --- Main dependencies
      type(type_dependency_id) :: id_temp
      type(type_dependency_id) :: id_sal
      type(type_dependency_id) :: id_density
      type(type_dependency_id) :: id_pressure
      type(type_dependency_id) :: id_hplus_prev      ! To retrieve calculated H+ values
      type(type_dependency_id) :: id_co2_wat         ! To retrieve calculated CO2 values from speciation

      ! --- Surface dependencies for air-sea exchange
      type(type_horizontal_dependency_id) :: id_wind
      type(type_horizontal_dependency_id) :: id_patm
      type(type_horizontal_dependency_id) :: id_co2_air
      type(type_horizontal_dependency_id) :: id_ice

      ! --- Diagnostics
      type(type_diagnostic_variable_id) :: id_ph
      type(type_diagnostic_variable_id) :: id_hplus
      type(type_diagnostic_variable_id) :: id_co2
      type(type_diagnostic_variable_id) :: id_hco3
      type(type_diagnostic_variable_id) :: id_co3
      type(type_diagnostic_variable_id) :: id_omega_ca

      ! --- Surface diagnostics
      type(type_surface_diagnostic_variable_id) :: id_co2_flux
      type(type_surface_diagnostic_variable_id) :: id_pco2_sea
      type(type_surface_diagnostic_variable_id) :: id_delta_pco2
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type type_carbonate_chemistry

contains

   subroutine initialize(self, configunit)
      class(type_carbonate_chemistry), intent(inout), target :: self
      integer,                         intent(in)            :: configunit

      call self%register_implemented_routines((/source_do, source_do_surface/))

      ! ---------------- State variables ----------------
      ! --- DIC ---
      call self%register_state_variable(self%id_dic, 'dic', 'mmol C m-3', 'Dissolved inorganic carbon', initial_value=2100.0_rk, minimum=0.0_rk, &
                                        standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%set_variable_property(self%id_dic, 'is_solute', .true.)
      call self%set_variable_property(self%id_dic, 'diff_method', DIFF_ION_LINEAR)
      ! Using bicarbonate ion as the representative ion to compute molecular diffusivity
      call self%set_variable_property(self%id_dic, 'm0', m0_HCO3)
      call self%set_variable_property(self%id_dic, 'm1', m1_HCO3)

      ! --- Alkalinity ---
      call self%register_state_variable(self%id_alk,'alk','mmol eq m-3','Total alkalinity', initial_value = 2300.0_rk, minimum=0.0_rk, &
                                        standard_variable=standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%set_variable_property(self%id_alk, 'is_solute', .true.)
      call self%set_variable_property(self%id_alk, 'diff_method', DIFF_ION_LINEAR)
      ! Using bicarbonate ion as the representative ion to compute molecular diffusivity
      call self%set_variable_property(self%id_alk, 'm0', m0_HCO3)
      call self%set_variable_property(self%id_alk, 'm1', m1_HCO3)

      ! ---------------- Optional tracer couplings ----------------

      call self%register_state_dependency(self%id_po4, 'po4', 'mmol m-3', 'Phosphate', required=.false.)
      call self%register_state_dependency(self%id_si , 'si' , 'mmol m-3', 'Silicate', required=.false.)
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mmol m-3', 'Ammonium', required=.false.)
      call self%register_state_dependency(self%id_h2s, 'h2s', 'mmol m-3', 'Hydrogen sulfide', required=.false.)

      ! ---------------- Env dependencies ----------------
      call self%register_dependency(self%id_temp,    standard_variables%temperature)
      call self%register_dependency(self%id_sal,     standard_variables%practical_salinity)
      call self%register_dependency(self%id_density, standard_variables%density)
      call self%register_dependency(self%id_pressure,standard_variables%pressure)

      ! ---------------- Surface dependencies ----------------
      call self%register_dependency(self%id_wind,    standard_variables%wind_speed)
      call self%register_dependency(self%id_patm,    standard_variables%surface_air_pressure)
      call self%register_dependency(self%id_co2_air, standard_variables%mole_fraction_of_carbon_dioxide_in_air)
      call self%register_dependency(self%id_ice,     standard_variables%ice_area_fraction, required=.false.)

      ! ---------------- Diagnostics ----------------

      call self%register_diagnostic_variable(self%id_ph,      'ph',       '1',            'pH')
      call self%register_diagnostic_variable(self%id_hplus,   'hplus',    'mmol m-3',     'Hydrogen ion concentration')
      call self%register_diagnostic_variable(self%id_co2,     'co2star',  'mmol C m-3',   'Dissolved CO2* concentration (CO2+H2CO3)')
      call self%register_diagnostic_variable(self%id_hco3,    'hco3',     'mmol C m-3',   'Bicarbonate concentration')
      call self%register_diagnostic_variable(self%id_co3,     'co3',      'mmol C m-3',   'Carbonate concentration')
      call self%register_diagnostic_variable(self%id_omega_ca,'omega_ca', '1',            'Calcite saturation state')

      call self%register_surface_diagnostic_variable(self%id_co2_flux, 'co2_flux', 'mmol C m-2 d-1', 'Air-sea CO2 flux')
      call self%register_surface_diagnostic_variable(self%id_pco2_sea, 'pco2_sea',     'uatm',       'Surface ocean pCO2')
      call self%register_surface_diagnostic_variable(self%id_delta_pco2, 'delta_pco2', 'uatm',       'Delta pCO2')

      ! Dependency on hplus (diagnostic variable) to retrieve its previous value
      call self%register_dependency(self%id_hplus_prev, 'hplus', 'mmol m-3', 'Previous hydrogen ion concentration')
      call self%register_dependency(self%id_co2_wat,    'co2star', 'mmol C m-3',   'Dissolved CO2* concentration')


   end subroutine initialize


   subroutine do(self, _ARGUMENTS_DO_)
      class(type_carbonate_chemistry), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: dic, alk, temp, sal, rho, pres
      real(rk) :: po4, si, nh4, h2s
      real(rk) :: dic_kg, alk_kg, po4_kg, si_kg, nh4_kg, h2s_kg

      real(rk) :: tb, ts, tf
      real(rk) :: k1, k2, kb, kw, ks, kf
      real(rk) :: kp1, kp2, kp3, ksi, knh4, kh2s
      real(rk) :: ksp_ca, ca

      real(rk) :: h, ph
      real(rk) :: h_prev
      real(rk) :: co2_kg, hco3_kg, co3_kg
      real(rk) :: omega_ca
      real(rk) :: denom

      _LOOP_BEGIN_

         ! ---------------- Read state and environment ----------------
         _GET_(self%id_dic, dic)           ! mmol m-3
         _GET_(self%id_alk, alk)           ! mmol eq m-3
         _GET_(self%id_temp, temp)         ! °C
         _GET_(self%id_sal, sal)           ! PSU
         _GET_(self%id_density, rho)       ! kg m-3
         _GET_(self%id_pressure, pres)     ! dbar

         rho = max(rho, 1.0_rk)

         po4 = 0.0_rk
         si  = 0.0_rk
         nh4 = 0.0_rk
         h2s = 0.0_rk
         ! If coupled, then retrieve actual concentrations for these tracers [mmol m-3]
         if (_AVAILABLE_(self%id_po4)) _GET_(self%id_po4, po4)
         if (_AVAILABLE_(self%id_si )) _GET_(self%id_si , si )
         if (_AVAILABLE_(self%id_nh4)) _GET_(self%id_nh4, nh4)
         if (_AVAILABLE_(self%id_h2s)) _GET_(self%id_h2s, h2s)

         ! --- Convert from mmol m-3 to mol kg-1
         dic_kg = dic / (1000.0_rk * rho)
         alk_kg = alk / (1000.0_rk * rho)
         po4_kg = po4 / (1000.0_rk * rho)
         si_kg  = si  / (1000.0_rk * rho)
         nh4_kg = nh4 / (1000.0_rk * rho)
         h2s_kg = h2s / (1000.0_rk * rho)

         ! Retrieve previous value for H+         
         _GET_(self%id_hplus_prev, h_prev)       ! mmol m-3 from previous timestep
         h_prev = h_prev / (1000.0_rk * rho)     ! convert to mol kg-1

         ! ---------------- Background seawater totals ----------------
         call compute_background_totals_from_salinity(sal, tb, ts, tf)

         ! ---------------- Equilibrium constants ----------------
         call compute_equilibrium_constants(temp, sal, pres, ts, tf, &
                                            k1, k2, kb, kw, ks, kf,  &
                                            kp1, kp2, kp3, ksi, knh4, kh2s, ksp_ca)

         ! ------------ Solve H+ to compute pH -------------------------
         h = solve_hplus(dic_kg, alk_kg, po4_kg, si_kg, nh4_kg, h2s_kg, &
                         tb, ts, tf, &
                         k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s, h_prev)

         ph = -log10(max(h, eps_h))

         ! ---------------- Carbonate speciation ----------------
         denom = h*h + k1*h + k1*k2
         co2_kg  = dic_kg * h*h     / denom
         hco3_kg = dic_kg * k1*h    / denom
         co3_kg  = dic_kg * k1*k2   / denom

         ! ---- Calcite saturation state
         ! Approximate [Ca2+] from reference seawater value scaled by salinity
         ca = calcon_ref * sal / 35.0_rk
         omega_ca = ca * co3_kg / (ksp_ca + eps_h)

         ! --- Save diagnostics converting from mol/kg to mmol m-3 ----------------
         _SET_DIAGNOSTIC_(self%id_ph, ph)
         _SET_DIAGNOSTIC_(self%id_hplus, 1000.0_rk * rho * h)
         _SET_DIAGNOSTIC_(self%id_co2,   1000.0_rk * rho * co2_kg)
         _SET_DIAGNOSTIC_(self%id_hco3,  1000.0_rk * rho * hco3_kg)
         _SET_DIAGNOSTIC_(self%id_co3,   1000.0_rk * rho * co3_kg)
         _SET_DIAGNOSTIC_(self%id_omega_ca, omega_ca)

      _LOOP_END_

      ! call test_carbonate_solver()
      ! stop 'carbonate solver test finished'

   end subroutine do


   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class(type_carbonate_chemistry), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: temp, sal, rho
      real(rk) :: wind, patm, co2_air, ice
      real(rk) :: co2_air_molmol
      real(rk) :: co2_sfc, co2_kg

      real(rk) :: tk, temp100
      real(rk) :: ln_k0, k0
      real(rk) :: pco2_sea

      real(rk) :: temp1, temp2, temp3, temp4
      real(rk) :: sch_co2, kgas
      real(rk) :: vapor_sw, patm_atm
      real(rk) :: pco2_air_atm, delta_pco2_atm
      real(rk) :: fugcoeff, fco2_air_atm
      real(rk) :: B, delta, x2
      real(rk) :: flux_co2
      real(rk) :: ice_frac

      real(rk), parameter :: pa_to_atm    = 1.0_rk / 101325.0_rk
      real(rk), parameter :: cm_hr_to_m_s = 0.01_rk / 3600.0_rk
      real(rk), parameter :: secs_per_day = 86400.0_rk

      _SURFACE_LOOP_BEGIN_

         ! ---------------- Read environment ----------------
         _GET_(self%id_temp, temp)               ! degC
         _GET_(self%id_sal, sal)                 ! PSU
         _GET_(self%id_density, rho)             ! kg m-3

         _GET_SURFACE_(self%id_wind, wind)       ! m s-1
         _GET_SURFACE_(self%id_patm, patm)       ! Pa
         _GET_SURFACE_(self%id_co2_air, co2_air) ! ppm

         co2_air_molmol = co2_air * 1.0e-6_rk    ! Convert atmospheric CO2 mole fraction from ppm to mol mol-1

         ice = 0.0_rk
         if (_AVAILABLE_HORIZONTAL_(self%id_ice)) _GET_SURFACE_(self%id_ice, ice)

         rho = max(rho, 1.0_rk)

         ! ---------------- Retrieve diagnosed CO2* value ----------------
         ! co2_sfc is stored as mmol C m-3, convert to mol kg-1
         _GET_(self%id_co2_wat, co2_sfc)
         co2_kg = co2_sfc / (1000.0_rk * rho)

         ! ------------------------------------------------------------
         ! K0: CO2 solubility (Weiss 1974), Adapted from CO2Sys
         ! ------------------------------------------------------------
         tk      = temp + 273.15_rk
         temp100 = tk / 100.0_rk

         ln_k0 = -60.2409_rk + 93.4517_rk / temp100 + 23.3585_rk * log(temp100) + &
                  sal * (0.023517_rk - 0.023656_rk * temp100 + 0.0047036_rk * temp100 * temp100)

         k0 = exp(ln_k0)   ! mol kg-1 atm-1

         ! ---------------- Surface ocean pCO2 ----------------
         ! co2_kg is CO2* in mol kg-1, convert to uatm
         pco2_sea = 1.0e6_rk * co2_kg / (k0 + eps_h)


         temp1 = min(40.0_rk, max(-2.0_rk, temp))     ! Limit temperature to Schmidt polynomial validity range
         temp2 = temp1 * temp1
         temp3 = temp2 * temp1
         temp4 = temp2 * temp2

         ! --------------------------------------------------------------
         ! Schmidt number for uptake of CO2 in seawater
         ! Least squares fourth-order polynomial fit of Schmidt number versus temperature for seawater (35‰) at 
         ! temperatures from –2°C to 40°C.
         ! Wanninkhof (2014), Relationship between wind speed and gas exchange over the ocean revisited
         ! --------------------------------------------------------------
         sch_co2 = 2116.8_rk - 136.25_rk * temp1 + 4.7353_rk * temp2 - &
                   0.092307_rk * temp3 + 0.0007555_rk * temp4

         ! Gas transfer velocity at Sc = 660, Wanninkhof (2014)
         ! k660 [cm h-1] = 0.251 ws^2
         kgas = 0.251_rk * wind * wind       ! cm h-1
         kgas = kgas * cm_hr_to_m_s          ! m s-1
         ! CO2 exchange velocity (m/s)
         kgas = kgas * sqrt(660.0_rk / max(sch_co2, eps_h))         

         ! ---------------- Atmospheric pCO2 ----------------
         ! Atmospheric pressure from Pa to atm
         patm_atm = patm * pa_to_atm

         ! Water vapour pressure over seawater (atm)
         vapor_sw = exp(24.4543_rk - 67.4509_rk * (100.0_rk / tk) - &
                        4.8489_rk * log(tk / 100.0_rk) - 0.000544_rk * sal)

         ! co2_air_molmol is a dry-air mole fraction (mol mol-1)
         ! Convert to moist-air partial pressure in atm
         pco2_air_atm = co2_air_molmol * max(patm_atm - vapor_sw, 0.0_rk)     ! partial pressure of CO2 in moist air

         ! Fugacity correction for atmospheric CO2
         B     = -1636.75_rk + 12.0408_rk * tk - 0.0327957_rk * tk*tk + &
                  3.16528e-5_rk * tk*tk*tk
         delta = 57.7_rk - 0.118_rk * tk
         x2    = (1.0_rk - pco2_air_atm)**2

         fugcoeff = exp(patm_atm * (B + 2.0_rk * delta * x2) / (82.05736_rk * tk))

         ! Fugacity of atmospheric CO2 (atm)
         fco2_air_atm = pco2_air_atm * fugcoeff

         ! Air-sea CO2 driving difference in atm (air fugacity minus sea pCO2)
         ! pco2_sea is in micro-atm, so it is converted to atm
         delta_pco2_atm = fco2_air_atm - pco2_sea * 1.0e-6_rk

         ! ---------------- Air-sea CO2 flux ----------------
         ! Positive flux means transfer from air to sea
         flux_co2 = 1000.0_rk * kgas * rho * k0 * delta_pco2_atm   ! mmol C m-2 s-1

         ! Reduce flux under sea ice
         ice_frac = min(max(ice, 0.0_rk), 1.0_rk)
         flux_co2 = flux_co2 * (1.0_rk - ice_frac)

         _ADD_SURFACE_FLUX_(self%id_dic, flux_co2)

         _SET_SURFACE_DIAGNOSTIC_(self%id_co2_flux, flux_co2 * secs_per_day)
         _SET_SURFACE_DIAGNOSTIC_(self%id_pco2_sea, pco2_sea)
         _SET_SURFACE_DIAGNOSTIC_(self%id_delta_pco2, 1.0e6_rk * delta_pco2_atm)

      _SURFACE_LOOP_END_

   end subroutine do_surface

end module carbonate_chemistry