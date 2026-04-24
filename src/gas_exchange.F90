#include "fabm_driver.h"

module gas_exchange

   use fabm_types
   implicit none
   private

   ! ------------------------------------------------------------------
   ! Gas exchange velocity at Schmidt number 660
   !
   ! Supported parameterisations:
   !   1 = Nightingale et al. (2000) 
   !       k600 [cm h-1] = 0.222 U10^2 + 0.333 U10
   !       converted internally here to k660
   !
   !   2 = OCMIP-2
   !       k660 [cm h-1] = 0.337 U10^2
   !
   !   3 = Wanninkhof (2014) - Parameterised using a 6h resolution
   !       k660 [cm h-1] = 0.251 U10^2  
   !
   ! Notes:
   ! - U10 is wind speed [m s-1]
   ! - Output k660 is stored in [m s-1]
   ! - W14 is formally tied to the second moment of short-interval winds.
   !   If daily-mean winds are used, applying k660 = 0.251 * U10^2 is an
   !   approximation because unresolved subdaily wind variance is ignored.
   ! - This module was adapted from MEDUSA-2.0 (Yool et al., 2013)
   ! ------------------------------------------------------------------

   integer, parameter, public :: eqn_nightingale    = 1
   integer, parameter, public :: eqn_ocmip2         = 2
   integer, parameter, public :: eqn_wanninkhof2014 = 3

   type, extends(type_base_model), public :: type_gas_exchange
      ! --- Dependencies
      type(type_horizontal_dependency_id) :: id_wind

      ! --- Diagnostics
      type(type_surface_diagnostic_variable_id) :: id_k660

      ! --- Parameters
      integer :: parameterization
   contains
      procedure :: initialize
      procedure :: do_surface
   end type type_gas_exchange

contains

   subroutine initialize(self, configunit)
      class(type_gas_exchange), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      ! Parameterisation choice:
      !   1 = Nightingale et al. (2000)
      !   2 = OCMIP-2
      !   3 = Wanninkhof (2014)
      call self%get_parameter(self%parameterization, 'equation', '-', &
                              'Gas exchange parameterisation: 1=Nightingale (2000), 2=OCMIP-2, 3=Wanninkhof 2014', &
                              default=eqn_wanninkhof2014, minimum=1, maximum=3)

      ! Surface wind speed forcing [m s-1]
      call self%register_dependency(self%id_wind, standard_variables%wind_speed)

      ! Gas exchange velocity at Sc = 660 [m s-1]
      call self%register_diagnostic_variable(self%id_k660, 'k660', 'm s-1', 'Gas exchange velocity normalized to Schmidt number 660')
   end subroutine initialize


   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class(type_gas_exchange), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: wind
      real(rk) :: k660_cmh, k_ref_cmh
      real(rk) :: k660

      real(rk), parameter :: cmh_to_ms = 1.0_rk / (100.0_rk * 3600.0_rk)

      _SURFACE_LOOP_BEGIN_

         _GET_SURFACE_(self%id_wind, wind)
         wind = max(wind, 0.0_rk)

         select case (self%parameterization)

         case (eqn_nightingale)
            ! Nightingale et al. (2000): parameterisation defined for Sc = 600 (k600)
            ! k600 [cm h-1] = 0.222 ws^2 + 0.333 ws
            k_ref_cmh = (0.222_rk * wind*wind) + 0.333_rk * wind
            ! Convert k600 to k660 using Schmidt scaling (k ∝ Sc^{-1/2})
            ! k_gas = k_ref * (Sc_ref / Sc_gas)^{1/2}
            k660_cmh = k_ref_cmh * (600.0_rk / 660.0_rk)**0.5_rk

         case (eqn_ocmip2)
            ! OCMIP-2
            ! k660 [cm h-1] = 0.337 ws^2
            k660_cmh = 0.337_rk * wind * wind

         case (eqn_wanninkhof2014)
            ! Wanninkhof (2014)
            ! k660 [cm h-1] = 0.251 ws^2
            k660_cmh = 0.251_rk * wind * wind

         end select

         ! Convert from cm h-1 to m s-1
         k660 = k660_cmh * cmh_to_ms

         _SET_SURFACE_DIAGNOSTIC_(self%id_k660, k660)

      _SURFACE_LOOP_END_

   end subroutine do_surface

end module gas_exchange