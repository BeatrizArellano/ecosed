#include "fabm_driver.h"

!-------------------------------------------------------------------------------------------------------
! Shared definitions for the pelagic ecosystem.
!
! This module does not represent a FABM model and is not instantiated in fabm.yaml.
! It defines aggregate variables that can receive contributions from multiple
! pelagic model instances and can subsequently be used as dependencies by
! other ecosystem modules.
!-------------------------------------------------------------------------------------------------------
module pelagic_common

   use fabm_types

   implicit none
   private

   public :: pom_production_n
   public :: total_chlorophyll

   !----------------------------------------------------------------------------------------------------
   ! Total production of particulate organic matter by pelagic biological processes.
   !
   ! Individual phytoplankton and zooplankton model instances can contribute
   ! their local POM-production diagnostics to this aggregate variable.
   !----------------------------------------------------------------------------------------------------
   type(type_interior_standard_variable), parameter :: pom_production_n = &
         type_interior_standard_variable(name='pom_production_n', &
                                         units='mmol N m-3 s-1', &
                                         aggregate_variable=.true.)


   !----------------------------------------------------------------------------------------------------
   ! Total chlorophyll concentration across all phytoplankton functional groups.
   !----------------------------------------------------------------------------------------------------
   type(type_interior_standard_variable), parameter :: total_chlorophyll = &
          type_interior_standard_variable(name='total_chlorophyll', &
                                          units='mg m-3', &
                                          aggregate_variable=.true.)

end module pelagic_common