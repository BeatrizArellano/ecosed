module molecular_diff
    use fabm_types

    ! molecular_diff contains the parameters needed to compute in-situ molecular diffusivities (for solutes)
    ! depending on Temperature and salinity.
    ! These methods and parameters are based on Chapter 4 in Boudreau (1997) 
    ! and verified with the R package marelac.

    ! Methods to compute diffusivities
    integer, parameter, public :: DIFF_NONE=0, DIFF_O2CO2_AB=1, DIFF_ION_LINEAR=2, &
                                  DIFF_ARRHENIUS=3, DIFF_WILKE_CHANG=4, DIFF_STOKES_EINSTEIN=5

    ! Oxygen (method DIFF_O2CO2_AB)
    ! Boudreau (1997) 4.58–4.59: D° = (A + B*(TK/mu0))*1e-5 [cm^2/s] 
    real(rk), parameter, public :: A_O2 = 0.2604_rk
    real(rk), parameter, public :: B_O2 = 0.006383_rk

    ! NO3, NO2, NH4 and PO4 (method DIFF_ION_LINEAR)
    ! Coefficients m0, m1 are taken directly from Boudreau (1997) Tables 4.7 and 4.8. 
    real(rk), parameter, public :: m0_NO3 = 9.50_rk
    real(rk), parameter, public :: m1_NO3 = 0.388_rk
    ! NO2
    real(rk), parameter, public :: m0_NO2 = 10.3_rk
    real(rk), parameter, public :: m1_NO2 = 0.331_rk
    ! NH4
    real(rk), parameter, public :: m0_NH4 = 9.50_rk
    real(rk), parameter, public :: m1_NH4 = 0.413_rk
    ! PO4
    real(rk), parameter, public :: m0_PO4 = 2.62_rk
    real(rk), parameter, public :: m1_PO4 = 0.143_rk

    !--- Carbonate chemistry
    ! CO2   (method DIFF_O2CO2_AB)
    ! Eq. 4.59 (Boudreau 1997)
    real(rk), parameter, public :: A_CO2 = 0.1954_rk
    real(rk), parameter, public :: B_CO2 = 0.005089_rk
    ! HCO3  (method DIFF_ION_LINEAR)
    ! Table 4.8 (Boudreau 1997)
    real(rk), parameter, public :: m0_HCO3 = 5.06_rk
    real(rk), parameter, public :: m1_HCO3 = 0.275_rk
    ! CO3  (method DIFF_ION_LINEAR)
    ! Table 4.8 (Boudreau 1997)
    real(rk), parameter, public :: m0_CO3 = 4.33_rk
    real(rk), parameter, public :: m1_CO3 = 0.199_rk


end module molecular_diff