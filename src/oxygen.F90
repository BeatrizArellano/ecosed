#include "fabm_driver.h"

module oxygen

   use fabm_types
   use molecular_diff, only: DIFF_O2CO2_AB, A_O2, B_O2
   implicit none
   private

   ! ------------------------------------------------------------------
   ! Air-sea oxygen exchange
   !
   ! Inputs:
   !   - wind speed [m s-1]
   !   - in situ temperature [degC]
   !   - practical salinity [PSU]
   !   - dissolved oxygen [mmol m-3]
   !
   ! Computes:
   !   - Schmidt number for O2
   !   - k660 from wind speed using Wanninkhof (2014)
   !   - gas transfer velocity for O2
   !   - oxygen saturation concentration
   !   - air-sea oxygen flux
   !
   ! Flux sign convention: positive into the ocean
   !
   ! Notes:
   ! - O2 saturation follows García and Gordon (1992), as in MEDUSA2.0.
   ! - Gas transfer scaling follows:
   !       kO2 = k660 * (Sc_O2 / 660)^(-1/2)
   ! ------------------------------------------------------------------

   type, extends(type_base_model), public :: type_oxygen
      ! --- State variables
      type(type_state_variable_id) :: id_o2

      ! --- Environmental dependencies
      type(type_dependency_id) :: id_temp
      type(type_dependency_id) :: id_sal

      ! --- Surface dependencies
      type(type_horizontal_dependency_id) :: id_wind
      type(type_horizontal_dependency_id) :: id_ice_fraction

      ! --- Diagnostics
      type(type_surface_diagnostic_variable_id) :: id_o2sat
      type(type_surface_diagnostic_variable_id) :: id_o2flux

      ! --- Parameters
      logical  :: apply_ice_cover
   contains
      procedure :: initialize
      procedure :: do_surface
   end type type_oxygen

contains

   subroutine initialize(self, configunit)
      class(type_oxygen), intent(inout), target :: self
      integer,            intent(in)            :: configunit

      call self%get_parameter(self%apply_ice_cover, 'apply_ice_cover', '-', 'Reduce gas exchange by sea-ice fraction', default=.false.)

      ! Oxygen
      call self%register_state_variable(self%id_o2, 'o2', 'mmol m-3', 'Dissolved oxygen', initial_value=300.0_rk, minimum=0.0_rk)
      call self%set_variable_property(self%id_o2, 'is_solute', .true.)
      call self%set_variable_property(self%id_o2, 'diff_method', DIFF_O2CO2_AB)
      call self%set_variable_property(self%id_o2, 'A', A_O2)
      call self%set_variable_property(self%id_o2, 'B', B_O2)


      ! Environmental dependencies
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_sal,  standard_variables%practical_salinity)
      call self%register_dependency(self%id_wind, standard_variables%wind_speed)
      if (self%apply_ice_cover) then
        ! Optional sea-ice fraction
        call self%register_dependency(self%id_ice_fraction, standard_variables%ice_area_fraction, required=.true.)
      end if
      
      ! Diagnostics
      call self%register_diagnostic_variable(self%id_o2sat,  'o2_sat', 'mmol m-3', 'Oxygen saturation concentration')
      call self%register_diagnostic_variable(self%id_o2flux, 'o2_flux', 'mmol m-2 d-1', 'Air-sea oxygen flux, positive into ocean')

   end subroutine initialize


   !-------------------------------------------------
   ! Compute oxygen saturation concentration (mmol m-3)
   ! and air-sea oxygen flux (mmol m-2 s-1) at 1 atm
   ! from temperature, salinity, and gas transfer velocity.
   ! Valid range for O2 saturation formulation:
   !   -1.9 <= T <= 40 °C,  0 <= S <= 42 PSU
   !
   ! O2 saturation follows García and Gordon (1992),
   ! Air-sea exchange is scaled from k660 using the
   ! Schmidt number for O2, with positive flux defined
   ! into the ocean.
   !-------------------------------------------------
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class(type_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: o2, temp, sal
      real(rk) :: wind, k_o2, sc_o2, k660
      real(rk) :: o2_sat, o2_flux
      real(rk) :: ice_fraction, open_water

      ! Garcia & Gordon (1992) coefficients for O2 saturation
      real(rk), parameter :: a0 = 2.00907_rk
      real(rk), parameter :: a1 = 3.22014_rk
      real(rk), parameter :: a2 = 4.05010_rk
      real(rk), parameter :: a3 = 4.94457_rk
      real(rk), parameter :: a4 = -2.56847e-1_rk
      real(rk), parameter :: a5 = 3.88767_rk
      real(rk), parameter :: b0 = -6.24523e-3_rk
      real(rk), parameter :: b1 = -7.37614e-3_rk
      real(rk), parameter :: b2 = -1.03410e-2_rk
      real(rk), parameter :: b3 = -8.17083e-3_rk
      real(rk), parameter :: c0 = -4.88682e-7_rk

      real(rk) :: ts
      real(rk) :: ln_cstar
      real(rk) :: o2_ml_l

      real(rk), parameter :: ml_per_mol  = 22391.6_rk
      real(rk), parameter :: sec_per_day = 86400.0_rk
      real(rk), parameter :: cmh_to_ms   = 1.0_rk / (100.0_rk * 3600.0_rk)

      _SURFACE_LOOP_BEGIN_

         _GET_(self%id_o2,   o2)
         _GET_(self%id_temp, temp)
         _GET_(self%id_sal,  sal)

         _GET_SURFACE_(self%id_wind, wind)       ! m s-1
         if (self%apply_ice_cover) then
            _GET_SURFACE_(self%id_ice_fraction, ice_fraction)
         else
            ice_fraction = 0.0_rk
         end if

         ice_fraction = min(max(ice_fraction, 0.0_rk), 1.0_rk)
         open_water   = 1.0_rk - ice_fraction

         ! --------------------------------------------------------------
         ! Schmidt number for uptake of O2 in seawater
         ! Least squares fourth-order polynomial fit of Schmidt number versus temperature for seawater at 
         ! temperatures from –2°C to 40°C.
         ! Wanninkhof (2014), Relationship between wind speed and gas exchange over the ocean revisited
         ! --------------------------------------------------------------
         sc_o2 = 1920.4_rk - 135.6_rk*temp + 5.2122_rk*temp*temp &
                - 0.10939_rk*temp**3 + 9.3777e-4_rk*temp**4

         ! Gas transfer velocity at Sc = 660, Wanninkhof (2014)
         ! k660 [cm h-1] = 0.251 ws^2
         k660 = 0.251_rk * wind * wind
         k660 = k660 * cmh_to_ms    ! convert to m s-1

         ! Oxygen exchange velocity (m/s)
         k_o2 = k660 * (660.0_rk / sc_o2)**0.5_rk

         ! --------------------------------------------------------------
         ! O2 saturation concentration (Garcia & Gordon, 1992)
         ! temp in degC, sal in PSU
         ! formula returns ml l-1, then converted to mmol m-3
         ! --------------------------------------------------------------
         ts = log((298.15_rk - temp) / (273.15_rk + temp))

         ln_cstar = a0 + a1*ts + a2*ts**2 + a3*ts**3 + a4*ts**4 + a5*ts**5 &
                  + sal*(b0 + b1*ts + b2*ts**2 + b3*ts**3) + c0*sal*sal

         o2_ml_l = exp(ln_cstar)

         ! Saturation concentration
         ! ml l-1 -> mmol m-3
         o2_sat = (o2_ml_l / ml_per_mol) * 1.0e6_rk

         ! Positive into ocean when undersaturated
         o2_flux = open_water * k_o2 * (o2_sat - o2)

         _ADD_SURFACE_FLUX_(self%id_o2, o2_flux)

         _SET_SURFACE_DIAGNOSTIC_(self%id_o2sat,  o2_sat)
         _SET_SURFACE_DIAGNOSTIC_(self%id_o2flux, o2_flux * sec_per_day)

      _SURFACE_LOOP_END_

   end subroutine do_surface

end module oxygen