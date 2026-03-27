add_library(fabm_models_ecosed OBJECT
            model_library.F90
            src/gas_exchange.F90
            src/molecular_diff.F90
            src/nitrogen.F90
            src/organic_matter.F90
            src/oxygen.F90
            src/pelagic_ecosystem.F90            
           )

target_link_libraries(fabm_models_ecosed PRIVATE fabm_base)
