#include "fabm_driver.h"

! ---------------------------------------------------------------------
! Carbonate chemistry module
!
! Prognostic state variables:
!   - DIC : dissolved inorganic carbon      [mmol C  m-3]
!   - ALK : total alkalinity                [mmol eq m-3]
! ---------------------------------------------------------------------
module carbonate_chemistry

   use fabm_types
   use molecular_diff, only: DIFF_ION_LINEAR, m0_HCO3, m1_HCO3
   implicit none
   private

   type, extends(type_base_model), public :: type_carbonate_chemistry

      ! --- State variables
      type(type_state_variable_id) :: id_dic
      type(type_state_variable_id) :: id_alk

      ! --- Future optional couplings/state dependencies
      ! type(type_state_variable_id) :: id_po4
      ! type(type_state_variable_id) :: id_si
      ! type(type_state_variable_id) :: id_nh4
      ! type(type_state_variable_id) :: id_h2s

      ! --- Future dependencies
      ! type(type_dependency_id) :: id_temp
      ! type(type_dependency_id) :: id_sal
      ! type(type_dependency_id) :: id_density
      ! type(type_dependency_id) :: id_pressure

      ! --- Future diagnostics
      ! type(type_diagnostic_variable_id) :: id_ph
      ! type(type_diagnostic_variable_id) :: id_co2
      ! type(type_diagnostic_variable_id) :: id_hco3
      ! type(type_diagnostic_variable_id) :: id_co3

   contains
      procedure :: initialize
      !procedure :: do
   end type type_carbonate_chemistry

contains

   subroutine initialize(self, configunit)
      class(type_carbonate_chemistry), intent(inout), target :: self
      integer,                         intent(in)            :: configunit

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


   end subroutine initialize


   !subroutine do(self, _ARGUMENTS_DO_)
   !   class(type_carbonate_chemistry), intent(in) :: self
   !   _DECLARE_ARGUMENTS_DO_
!
   !   _LOOP_BEGIN_
!
   !      ! No carbonate source/sink terms yet.
   !      ! Chemistry, diagnostics and air-sea exchange will be added later.
!
   !   _LOOP_END_
!
   !end subroutine do

end module carbonate_chemistry