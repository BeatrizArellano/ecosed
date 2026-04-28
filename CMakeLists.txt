add_library(fabm_models_ecosed OBJECT
            model_library.F90
            src/carbonate_chemistry.F90
            src/carbonate_system.F90
            src/iron.F90
            src/manganese.F90
            src/molecular_diff.F90
            src/nitrogen.F90
            src/om_degradation.F90
            src/oxygen.F90
            src/pelagic_ecosystem.F90 
            src/sulfur.F90           
           )

target_link_libraries(fabm_models_ecosed PRIVATE fabm_base)
