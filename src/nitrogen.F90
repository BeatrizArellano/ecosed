#include "fabm_driver.h"

module nitrogen

   use fabm_types
   use molecular_diff, only: DIFF_ION_LINEAR, m0_NO3, m1_NO3, m0_NO2, m1_NO2, m0_NH4, m1_NH4
   implicit none
   private



   type, extends(type_base_model), public :: type_nitrogen
      ! --- State variables
      type(type_state_variable_id) :: id_no3

   contains
      procedure :: initialize
   end type type_nitrogen

contains

   subroutine initialize(self, configunit)
      class(type_nitrogen), intent(inout), target :: self
      integer,              intent(in)            :: configunit

      ! ---------------- State variables ----------------
      ! --- NO3 ---
      call self%register_state_variable(self%id_no3, 'no3', 'mmol m-3', 'Dissolved Nitrate', initial_value=5.0_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_no3, 'is_solute', .true.)
      call self%set_variable_property(self%id_no3, 'diff_method', DIFF_ION_LINEAR)
      call self%set_variable_property(self%id_no3, 'm0', m0_NO3)
      call self%set_variable_property(self%id_no3, 'm1', m1_NO3)
      
      ! Contribution to total Nitrogen
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_no3)

   end subroutine initialize
end module nitrogen