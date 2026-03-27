module ecosed_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use gas_exchange
   use nitrogen
   use organic_matter
   use oxygen
   use pelagic_ecosystem

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: ecosed_model_factory

contains

   subroutine create(self, name, model)

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('gas_exchange');    allocate (type_gas_exchange::model)
         case ('nitrogen');        allocate (type_nitrogen::model)
         case ('organic_matter');  allocate (type_organic_matter::model)
         case ('oxygen');          allocate (type_oxygen::model)
         case ('pelagiceco');      allocate (type_pelagic_ecosystem::model)
         
         ! Add case statements for new models here
      end select

   end subroutine create

end module