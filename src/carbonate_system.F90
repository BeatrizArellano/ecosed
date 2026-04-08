!---------------------------------------------------------------------------------------------------
! Carbonate system routines for seawater biogeochemistry.
!
! This module provides the core routines needed to solve the marine carbonate system from dissolved
! inorganic carbon (DIC) and total alkalinity. It computes salinity-derived background totals
! (boron, sulfate, fluoride), evaluates equilibrium constants at in-situ temperature, salinity and
! pressure, and solves for hydrogen ion concentration [H+] on the total pH scale.
!
! The implementation uses a safeguarded Newton–Raphson solver for pH with dynamic bracketing and 
! fallback to bisection when the Newton iteration stagnates. 
! Uses CO2SYS-based formulations for equilibrium constants and pressure corrections, and 
! includes extended alkalinity contributions from phosphate, silicate, ammonia and sulfide.
! Internally, the solver is formulated on the total pH scale, while KS is retained on the free scale.
!---------------------------------------------------------------------------------------------------
module carbonate_system

   use fabm_types
   implicit none
   private

   ! ---------------- Module parameters ----------------
   real(rk), parameter, public :: eps_h = 1.0e-20_rk
   real(rk), parameter, public :: calcon_ref = 1.03e-2_rk   ! Reference [Ca2+] at S=35, mol kg-1
   real(rk), parameter         :: rgas = 83.14472_rk   
   integer,  parameter         :: max_iter_h = 20
   real(rk), parameter         :: h_rel_tol = 1.0e-4_rk

   public :: compute_background_totals_from_salinity
   public :: compute_equilibrium_constants
   public :: solve_hplus

contains

   ! Compute salinity-derived total concentrations of boron (TB),
   ! sulfate (TS), and fluoride (TF) in mol kg-1.
   subroutine compute_background_totals_from_salinity(sal, tb, ts, tf)
      real(rk), intent(in)  :: sal             ! PSU
      real(rk), intent(out) :: tb, ts, tf      
      real(rk) :: chlorinity

      ! Chlorinity from practical salinity
      ! Wooster et Aal. (1969)
      chlorinity = sal / 1.80655_rk

      ! Total borate [mol kg-1]
      ! Uppstrom, L., Deep-Sea Research 21:161-162, 1974
      tb = 0.0004157_rk * sal / 35.0_rk

      ! Total sulfate [mol kg-1]
      ! Morris and Riley (1966), via chlorinity
      ts = (0.14_rk / 96.062_rk) * chlorinity

      ! Total fluoride [mol kg-1]
      ! Riley (1965), via chlorinity
      tf = (0.000067_rk / 18.998_rk) * chlorinity

   end subroutine compute_background_totals_from_salinity


   ! ------------------------------------------------------------------------------------------
   ! Compute equilibrium constants for the carbonate system at in-situ T, S, and P.
   ! Constants are first evaluated in their native (free or SWS) scales, pressure-corrected,
   ! and then converted to a consistent total pH scale using factors derived from KS and KF.
   ! All returned constants are solver-ready and consistent with H+ defined on the total scale.
   ! KS is retained on its free-scale basis because it defines the conversion between free and total H scales.
   ! Units follow CO2SYS conventions (mol kg-1 or derived).
   ! Adapted from CO2Sys and RADIv2 for the inclusion of ammonia and sulfide. 
   ! ------------------------------------------------------------------------------------------
   subroutine compute_equilibrium_constants(temp, sal, pres, ts, tf, &
                                            k1, k2, kb, kw, ks, kf, &
                                            kp1, kp2, kp3, ksi, knh4, kh2s, ksp_ca)
      real(rk), intent(in)  :: temp, sal, pres, ts, tf
      real(rk), intent(out) :: k1, k2, kb, kw, ks, kf
      real(rk), intent(out) :: kp1, kp2, kp3, ksi, knh4, kh2s
      real(rk), intent(out) :: ksp_ca

      real(rk) :: tk, logtk, sqrt_sal, sal15
      real(rk) :: ion_s, sqrt_ion_s
      real(rk) :: pbar, rt
      real(rk) :: delta_v, kappa, ln_fac
      real(rk) :: ks_free, kf_free
      real(rk) :: sws_to_tot_np, free_to_tot_np
      real(rk) :: sws_to_tot_p, free_to_tot_p
      real(rk) :: ln_ks, ln_kf, ln_kb, ln_kw
      real(rk) :: ln_kp1, ln_kp2, ln_kp3, ln_ksi
      real(rk) :: pknh4, ln_kh2s
      real(rk) :: p_k1, p_k2
      real(rk) :: log_ksp_ca      

      tk         = temp + 273.15_rk           ! Celsius to Kelvin
      logtk      = log(tk)
      sqrt_sal   = sqrt(max(sal, 0.0_rk))
      sal15      = sal * sqrt_sal
      pbar       = pres / 10.0_rk - 1.0_rk    ! absolute dbar -> bar relative to the surface, ~0 at surface
      rt         = rgas * tk          
           
      ! Ionic strength
      ! Adapted from CO2sys: From the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4
      ion_s      = 19.924_rk * sal / (1000.0_rk - 1.005_rk * sal)
      sqrt_ion_s = sqrt(ion_s)

      ! ------------------------------------------------------------
      ! KS: Bisulfate dissociation constant (Dickson, 1990), on the free pH scale. 
      ! Adapted from CO2Sys
      ! ------------------------------------------------------------
      ln_ks = -4276.1_rk / tk + 141.328_rk - 23.093_rk * logtk + &
              (-13856.0_rk / tk + 324.57_rk - 47.986_rk * logtk) * sqrt_ion_s + &
              (35474.0_rk / tk - 771.54_rk + 114.723_rk * logtk) * ion_s - &
              (2698.0_rk/tk) * ion_s**1.5_rk + (1776.0_rk/tk) * ion_s * ion_s + &
              log(1.0_rk - 0.001005_rk * sal)    ! convert to mol/kg-SW
      ks_free = exp(ln_ks)      

      ! ------------------------------------------------------------
      ! KF: Fluoride dissociation constant (Dickson & Riley, 1990), on the free pH scale.
      ! Adapted from CO2Sys
      ! ------------------------------------------------------------
      ln_kf = 1590.2_rk / tk - 12.641_rk + 1.525_rk * sqrt_ion_s + &
              log(1.0_rk - 0.001005_rk * sal)
      kf_free = exp(ln_kf)

      ! ------------------------------------------------------------
      ! pH scale conversion factors non corrected for pressure
      ! The internal system is kept in the total scale
      ! ------------------------------------------------------------
      free_to_tot_np = 1.0_rk + ts/ks_free 
      sws_to_tot_np  = (1.0_rk + ts/ks_free) / (1.0_rk + ts/ks_free + tf/kf_free)

      ! Pressure correction for ks on the free scale (Millero, 1995). Adapted from CO2sys
      delta_v = -18.03_rk + 0.0466_rk * temp + 0.000316_rk * temp * temp
      kappa   = (-4.53_rk + 0.090_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      ks_free = ks_free * exp(ln_fac)

      ! Pressure correction for kf on the free scale (Millero, 1995). Adapted from CO2sys
      delta_v = -9.78_rk - 0.0090_rk * temp - 0.000942_rk * temp * temp
      kappa   = (-3.91_rk + 0.054_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      kf_free = kf_free * exp(ln_fac)

      ! ------------------------------------------------------------
      ! K1, K2: carbonic acid dissociation (Lueker, Dickson, Keeling, 2000)
      ! This fit is on the total scale
      ! Adapted from CO2Sys
      ! ------------------------------------------------------------
      p_k1 = 3633.86_rk/tk - 61.2172_rk + 9.6777_rk*logtk - 0.011555_rk*sal + 0.0001152_rk*sal*sal
      p_k2 = 471.78_rk /tk + 25.9290_rk - 3.16967_rk*logtk - 0.01781_rk*sal + 0.0001122_rk*sal*sal
      
      k1 = 10.0_rk ** (-p_k1)   ! In the total pH scale
      k2 = 10.0_rk ** (-p_k2)

      ! Convert to SWS scale for pressure correction
      k1 = k1 / sws_to_tot_np
      k2 = k2 / sws_to_tot_np

      ! Pressure correction for k1 (Millero, 1995)
      delta_v = -25.5_rk + 0.1271_rk * temp
      kappa   = (-3.08_rk + 0.0877_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      k1      = k1 * exp(ln_fac)

      ! Pressure correction for k2 (Millero, 1995), fit on the SWS scale
      delta_v = -15.82_rk - 0.0219_rk * temp
      kappa   = (1.13_rk - 0.1475_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      k2      = k2 * exp(ln_fac)

      ! ------------------------------------------------------------
      ! KB: Boric acid dissociation (Dickson, 1990), fit on the total pH scale
      ! ------------------------------------------------------------
      ln_kb = (-8966.90_rk - 2890.53_rk * sqrt_sal - 77.942_rk * sal + 1.728_rk * sal15 - 0.0996_rk*sal*sal) / tk + &
              (148.0248_rk + 137.1942_rk * sqrt_sal + 1.62142_rk * sal) + &
              (-24.4344_rk - 25.085_rk * sqrt_sal - 0.2474_rk * sal) * logtk + &
              0.053105_rk * sqrt_sal * tk
      kb = exp(ln_kb)       ! This is in the total pH scale

      kb = kb / sws_to_tot_np            ! Convert to SWS scale for pressure correction

      ! Pressure correction for kb (Millero, 1979), fit on the SWS pH scale
      delta_v = -29.48_rk + 0.1622_rk * temp - 0.002608_rk * temp * temp
      kappa   = -2.84_rk / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      kb      = kb * exp(ln_fac)

      ! ------------------------------------------------------------
      ! KW: water dissociation (Millero, 1995) fit on the SWS scale
      ! ------------------------------------------------------------
      ln_kw = -13847.26_rk/tk + 148.9802_rk - 23.6521_rk*logtk + &
              (118.67_rk/tk - 5.977_rk + 1.0495_rk*logtk) * sqrt_sal - 0.01615_rk*sal
      kw = exp(ln_kw)                  ! this is on the SWS pH scale in (mol/kg-SW)^2

      ! Pressure correction for kw (Millero, 1983), adapted from CO2sys
      delta_v = -20.02_rk + 0.1119_rk * temp - 0.001409_rk * temp * temp
      kappa   = (-5.13_rk + 0.0794_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      kw      = kw * exp(ln_fac)

      ! ------------------------------------------------------------
      ! Phosphate constants (Yao and Millero, 1995), fit on the SWS scale
      ! ------------------------------------------------------------
      ln_kp1 = -4576.752_rk/tk + 115.540_rk - 18.453_rk*logtk + &
               (-106.736_rk/tk + 0.69171_rk) * sqrt_sal + &
               (-0.65643_rk/tk - 0.01844_rk) * sal
      kp1 = exp(ln_kp1)                    ! this is on the SWS pH scale in mol/kg-H2O

      ln_kp2 = -8814.715_rk/tk + 172.1033_rk - 27.927_rk*logtk + &
               (-160.340_rk/tk + 1.3566_rk) * sqrt_sal + &
               (0.37335_rk/tk - 0.05778_rk) * sal
      kp2 = exp(ln_kp2)                   ! this is on the SWS pH scale in mol/kg-H2O

      ln_kp3 = -3070.75_rk/tk - 18.126_rk + &
               (17.27039_rk/tk + 2.81197_rk) * sqrt_sal + &
               (-44.99486_rk/tk - 0.09984_rk) * sal
      kp3 = exp(ln_kp3)                   ! this is on the SWS pH scale in mol/kg-H2O

      ! Pressure correction for Kp1, kp2 and kp3 on the SWS? scale (Millero, 1995). Adapted from CO2sys
      delta_v = -14.51_rk + 0.1211_rk * temp - 0.000321_rk * temp * temp
      kappa   = (-2.67_rk + 0.0427_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      kp1     = kp1 * exp(ln_fac)
      ! kp2
      delta_v = -23.12_rk + 0.1758_rk * temp - 0.002647_rk * temp * temp
      kappa   = (-5.15_rk + 0.090_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      kp2     = kp2 * exp(ln_fac)
      ! kp3
      delta_v = -26.57_rk + 0.2020_rk * temp - 0.003042_rk * temp * temp
      kappa   = (-4.08_rk + 0.0714_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      kp3     = kp3 * exp(ln_fac)

      ! ------------------------------------------------------------
      ! Silicate constant (Yao and Millero, 1995)
      ! ------------------------------------------------------------
      ln_ksi = -8904.2_rk / tk + 117.400_rk - 19.334_rk * logtk + &
               (-458.79_rk / tk + 3.5913_rk) * sqrt_ion_s + &
               (188.74_rk / tk - 1.5998_rk) * ion_s + &
               (-12.1652_rk / tk + 0.07871_rk) * ion_s * ion_s + &
               log(1.0_rk - 0.001005_rk * sal)
      ksi = exp(ln_ksi)                   ! this is on the SWS pH scale in mol/kg-H2O

      ! Pressure correction for ksi (Millero, 1995). Adapted from CO2sys that uses values from boric acid.
      delta_v = -29.48_rk + 0.1622_rk * temp - 0.002608_rk * temp * temp
      kappa   = -2.84_rk / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      ksi     = ksi * exp(ln_fac)

      ! ------------------------------------------------------------
      ! Ammonia dissociation constant (Yao and Millero, 1995), fit on the total scale
      ! Clegg and Whitfield style fit used in CO2SYS/RADI extension
      ! ------------------------------------------------------------
      pknh4 = 9.244605_rk - 2729.33_rk * (1.0_rk / 298.15_rk - 1.0_rk / tk) + &
              (0.04203362_rk - 11.24742_rk / tk) * sal**0.25_rk + &
              (-13.6416_rk + 1.176949_rk * sqrt(tk) - 0.02860785_rk * tk + 545.4834_rk / tk) * sqrt_sal + &
              (-0.1462507_rk + 0.0090226468_rk * sqrt(tk) - 0.0001471361_rk * tk + 10.5425_rk / tk) * sal**1.5_rk + &
              (0.004669309_rk - 0.0001691742_rk * sqrt(tk) - 0.5677934_rk / tk) * sal*sal + &
              (-2.354039e-05_rk + 0.009698623_rk / tk) * sal**2.5_rk
      knh4 = 10.0_rk ** (-pknh4)                         ! total scale, mol/kg-H2O
      knh4 = knh4 * (1.0_rk - 0.001005_rk * sal)         !  mol/kg-SW
      knh4 = knh4 / sws_to_tot_np                        ! converts to SWS pH scale

      ! Pressure correction adapted from RADI/CO2sys
      delta_v = -26.43_rk + 0.0889_rk * temp - 0.000905_rk * temp * temp
      kappa   = (-5.03_rk + 0.0814_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      knh4    = knh4 * exp(ln_fac)

      ! ------------------------------------------------------------
      ! Sulfide dissociation constant (Millero et al., 1988), fit on the total pH scale
      ! Adapted from RADI/CO2sys 
      ! ------------------------------------------------------------
      ln_kh2s = 225.838_rk - 13275.3_rk / tk - 34.6435_rk * logtk + 0.3449_rk * sqrt_sal - 0.0274_rk * sal
      kh2s = exp(ln_kh2s)                                ! total scale
      kh2s = kh2s / sws_to_tot_np                        ! convert to SWS pH scale

      ! Pressure correction adapted from RADI/CO2sys
      delta_v = -11.07_rk - 0.0090_rk * temp - 0.000942_rk * temp * temp
      kappa   = (-2.89_rk + 0.054_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      kh2s    = kh2s * exp(ln_fac)

      ! ------------------------------------------------------------
      ! Calcite solubility product (Mucci et al., 1983)
      ! ------------------------------------------------------------
      log_ksp_ca = -171.9065_rk - 0.077993_rk * tk + 2839.319_rk / tk + &
                   71.595_rk * logtk / log(10.0_rk) + &
                   (-0.77712_rk + 0.0028426_rk * tk + 178.34_rk / tk) * sqrt_sal - &
                   0.07711_rk * sal + 0.0041249_rk * sal15
      ksp_ca = 10.0_rk ** log_ksp_ca                    ! this is in (mol/kg-SW)^2

      ! Pressure correction for calcite (Millero, 1979), Adapted from CO2sys
      delta_v = -48.76_rk + 0.5304_rk * temp
      kappa   = (-11.76_rk + 0.3692_rk * temp) / 1000.0_rk
      ln_fac  = (-delta_v + 0.5_rk * kappa * pbar) * pbar / rt
      ksp_ca  = ksp_ca * exp(ln_fac)

      ! ------------------------------------------------------------
      ! Convert to the total pH scale
      ! ------------------------------------------------------------      
      ! The internal system is kept in the total scale
      
      free_to_tot_p = 1.0_rk + ts/ks_free 
      sws_to_tot_p  = (1.0_rk + ts/ks_free) / (1.0_rk + ts/ks_free + tf/kf_free)

      ks = ks_free                     ! ks is kept in the free scale because it is the bridge between Hydrogen scales  
      kf = kf_free * free_to_tot_p     ! Convert from the free-scale to total-scale (pressure corrected) 
      k1 = k1 * sws_to_tot_p           ! Convert from SWS to the total pH scale 
      k2 = k2 * sws_to_tot_p   
      kb = kb * sws_to_tot_p    
      kw = kw * sws_to_tot_p           ! Convert from the SWS to the Total scale  
      kp1 = kp1 * sws_to_tot_p
      kp2 = kp2 * sws_to_tot_p
      kp3 = kp3 * sws_to_tot_p       
      ksi = ksi * sws_to_tot_p
      knh4 = knh4 * sws_to_tot_p                    
      kh2s = kh2s * sws_to_tot_p

   end subroutine compute_equilibrium_constants


   !----------------------------------------------------------------------------------
   ! Calculates an initial guess for [H+] on the total pH scale.
   !
   ! Inputs:
   !   dic : dissolved inorganic carbon [mol kg-1]
   !   alk : total alkalinity            [mol eq kg-1]
   !   tb  : total boron                [mol kg-1]
   !   k1, k2, kb : equilibrium constants on the solver pH scale
   !
   ! This is a cold-start estimate based on the carbonate-borate alkalinity approximation.
   ! It is intended to provide a robust initial value for the full safeguarded solver,
   ! not to represent the complete alkalinity system.
   !--------------------------------------------------------------------------------
   real(rk) function initial_hplus_guess(dic, alk, tb, k1, k2, kb)
      real(rk), intent(in) :: dic, alk, tb, k1, k2, kb
      
      real(rk) :: ca1, ba1
      real(rk) :: disc, sqrtd, hmin
      real(rk) :: a2, a1, a0
      real(rk) :: alk_cb

      ! For the analytical initial guess, approximate total alkalinity as being
      ! dominated by the carbonate-borate system.
      alk_cb = alk

      ! Set extreme alkalinity cases using bounds from the carbonate-borate system
      if (alk_cb <= 0.0_rk) then
         initial_hplus_guess = 1.0e-3_rk    ! pH~3
         return
      else if (alk_cb >= (2.0_rk * dic + tb)) then
         initial_hplus_guess = 1.0e-10_rk   ! pH~10
         return
      end if

      ca1 = dic / (alk_cb + eps_h)            ! Normalised ratios of DIC and boron to alkalinity.
      ba1 = tb  / (alk_cb + eps_h)

      ! Coefficients of the cubic polynomial approximation of the
      ! carbonate–borate alkalinity equation as a function of H+.
      a2 = kb * (1.0_rk - ba1) + k1 * (1.0_rk - ca1)
      a1 = k1 * kb * (1.0_rk - ba1 - ca1) + k1 * k2 * (1.0_rk - 2.0_rk * ca1)
      a0 = k1 * k2 * kb * (1.0_rk - ba1 - 2.0_rk * ca1)

      ! Discriminant of the quadratic equation for finding the local minimum of the cubic. 
      disc = a2*a2 - 3.0_rk*a1

      if (disc > 0.0_rk) then
         sqrtd = sqrt(disc)
         ! Estimate the location of the local minimum of the polynomial.
         if (a2 < 0.0_rk) then
            hmin = (-a2 + sqrtd) / 3.0_rk
         else
            hmin = -a1 / (a2 + sqrtd + eps_h)
         end if

         ! Compute the root using a second-order Taylor expansion around the local minimum of the cubic approximation.
         initial_hplus_guess = hmin + sqrt(max(0.0_rk, -(a0 + hmin*(a1 + hmin*(a2 + hmin))) / (sqrtd + eps_h)))

      else
         ! If no real minimum exists, fall back to a neutral pH guess.
         initial_hplus_guess = 10.0_rk**(-8.1_rk)
      end if

      initial_hplus_guess = max(initial_hplus_guess, 1.0e-10_rk)

   end function initial_hplus_guess


   !----------------------------------------------------------------------------------
   !> Solve for hydrogen ion concentration [H+] (mol kg-1) from DIC and alkalinity.
   !!
   !! This routine computes [H+] by solving the nonlinear alkalinity balance equation:
   !!
   !!    ALK_tracer = ALK_speciation(H)
   !!
   !! where ALK_speciation(H) is the total alkalinity implied by chemical speciation
   !! of dissolved inorganic carbon and other acid-base systems at a given H+.
   !!
   !! The solver uses a safeguarded Newton method in relative/log(H+) form, adapted
   !! from the PISCES carbonate chemistry solver. The approach combines:
   !!
   !!   - A Newton-like update expressed as a relative change in H+ (log-space step)
   !!   - Dynamic bracketing to ensure the root remains within a physically valid interval
   !!   - A fallback to a geometric-mean (pH-midpoint) step when convergence slows
   !!
   !! INPUT:
   !!   dic      : Dissolved inorganic carbon              [mol kg-1]
   !!   alk      : Total alkalinity (tracer value)         [mol eq kg-1]
   !!   po4      : Total phosphate                         [mol kg-1]
   !!   si       : Total silicate                          [mol kg-1]
   !!   nh4      : Total ammonia                           [mol kg-1]
   !!   h2s      : Total sulfide                           [mol kg-1]
   !!   tb       : Total boron                             [mol kg-1]
   !!   ts       : Total sulfate                           [mol kg-1]
   !!   tf       : Total fluoride                          [mol kg-1]
   !!
   !!   k1, k2   : Carbonic acid dissociation constants
   !!   kb       : Boric acid dissociation constant
   !!   kw       : Water dissociation constant
   !!   ks       : Bisulfate dissociation constant (free scale)
   !!   kf       : Fluoride dissociation constant
   !!   kp1-3    : Phosphate dissociation constants
   !!   ksi      : Silicate dissociation constant
   !!   knh4     : Ammonia dissociation constant
   !!   kh2s     : Sulfide dissociation constant
   !!
   !!   h_init   : Initial guess for H+ (mol kg-1), typically from previous timestep
   !!
   !! OUTPUT:
   !!   solve_hplus : Hydrogen ion concentration [mol kg-1]
   !!                 Returns -1.0 if convergence is not achieved within max_iter_h
   !----------------------------------------------------------------------------------
   real(rk) function solve_hplus(dic, alk, po4, si, nh4, h2s, tb, ts, tf, &
                                 k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s, h_init)
      real(rk), intent(in) :: dic, alk, po4, si, nh4, h2s
      real(rk), intent(in) :: tb, ts, tf
      real(rk), intent(in) :: k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s
      real(rk), intent(in) :: h_init

      integer  :: iter
      real(rk) :: h, h_prev, h_lo, h_hi
      real(rk) :: alk_res, dalkdh, h_lnfactor, h_delta
      real(rk) :: alk_inf, alk_sup, free_to_tot, delta
      real(rk) :: res_absmin

      ! ---------------- Initial guess ----------------
      if (h_init > 0.0_rk) then
         h = h_init                  ! Use previous value as an initial guess, if available
      else
         h = initial_hplus_guess(dic, alk, tb, k1, k2, kb)
      end if

      ! Lower and upper alkalinity bounds from extreme protonation states.
      alk_inf = -po4 - ts - tf
      alk_sup = 2.0_rk * dic + tb + 2.0_rk * po4 + si + nh4 + h2s

      ! Conversion factor between free and total proton scales.
      free_to_tot = 1.0_rk + ts / ks

      ! Lower H+ bound associated with the minimum admissible alkalinity.
      delta = (alk - alk_inf)**2 + 4.0_rk * kw / free_to_tot
      if (alk >= alk_inf) then
         h_lo = 2.0_rk * kw / ((alk - alk_inf) + sqrt(delta))
      else
         h_lo = free_to_tot * (-(alk - alk_inf) + sqrt(delta)) / 2.0_rk
      end if

       ! Upper H+ bound associated with the maximum admissible alkalinity.
      delta = (alk - alk_sup)**2 + 4.0_rk * kw / free_to_tot
      if (alk <= alk_sup) then
         h_hi = free_to_tot * (-(alk - alk_sup) + sqrt(delta)) / 2.0_rk
      else
         h_hi = 2.0_rk * kw / ((alk - alk_sup) + sqrt(delta))
      end if

      ! Keeping the interval within broad physically reasonable limits.
      h_lo = max(h_lo, 1.0e-11_rk)
      h_hi = min(h_hi, 1.0e-3_rk)

      ! Guard against impossible cases.
      if (h_lo >= h_hi) then
         h_lo = 1.0e-11_rk
         h_hi = 1.0e-3_rk
      end if

      ! Ensure the initial guess lies inside the bracket.
      h = max(min(h_hi, h), h_lo)

      res_absmin = huge(1.0_rk)            ! To store the smallest absolute residual

      ! ------ Newton-Raphson iterations with a bisection-method fallback ----------------
      ! Approach adapted from PISCES
      do iter = 1, max_iter_h

         h_prev = h

         ! Alkalinity residual: difference between actual alkalinity and the 
         ! alkalinity that would result if water had the given concentration of H+
         ! Alk_res(H) = Alk_speciation(H) - Alk_tracer
         ! This is used to find the value for H+ that yields a residual close to zero. 
         alk_res = alk_total(h, dic, po4, si, nh4, h2s, tb, ts, tf, &
                             k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s) - alk

         ! Derivative of the alkalinity with respect to H+
         ! Slope of the tangent at the current guess
         dalkdh = dalk_total_dh(h, dic, po4, si, nh4, h2s, tb, ts, tf, &
                                k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s)

         ! Adapt bracketing interval using the sign of the residual, so the root stays inside
         if (alk_res > 0.0_rk) then
            ! If alkalinity from speciation is high, then the root should be at a larger H+
            h_lo = h_prev
         else if (alk_res < 0.0_rk) then
            ! If alkalinity from speciation is low, then the root should be at a smaller H+
            h_hi = h_prev
         end if

         if (abs(alk_res) >= 0.5_rk * res_absmin) then
            ! If the current residual is still at least half as large as the best one found so far, 
            ! perform a more conservative bisection step in pH space = geometric mean in H+
            h = sqrt(h_lo * h_hi)                     ! Geometric mean in H+, equivalent to midpoint in pH
            h_lnfactor = (h - h_prev) / h_prev        ! Relative change in H+
         else
            ! Newton method in relative/log(H+) space
            h_lnfactor = -alk_res / (dalkdh * h_prev)      ! relative change needed to reduce the residual to zero

            if (abs(h_lnfactor) > 1.0_rk) then
               h = h_prev * exp(h_lnfactor)                ! If the change is large, h is estimated in the log-space
            else
               h_delta = h_lnfactor * h_prev               ! If relative change is small, the correction is applied linearly
               h = h_prev + h_delta
            end if
            ! If the solution found lies outside the safe interval, find a safer intermediate root inside the interval. 
            if (h < h_lo) then
               h = sqrt(h_prev * h_lo)                      ! Geometric mean
               h_lnfactor = (h - h_prev) / h_prev
            end if

            if (h > h_hi) then
               h = sqrt(h_prev * h_hi)
               h_lnfactor = (h - h_prev) / h_prev
            end if
         end if

         res_absmin = min(abs(alk_res), res_absmin)

         if (abs(h_lnfactor) < h_rel_tol) then
            solve_hplus = max(h, eps_h)
            return
         end if

      end do

      ! If the maximum number of iterations has been reached without finding a root
      solve_hplus = -1.0_rk

   end function solve_hplus

   real(rk) function alk_total(h, dic, po4, si, nh4, h2s, tb, ts, tf, &
                               k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s)
      real(rk), intent(in) :: h, dic, po4, si, nh4, h2s
      real(rk), intent(in) :: tb, ts, tf
      real(rk), intent(in) :: k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s

      alk_total = alk_carbonate(h, dic, k1, k2) &
                + alk_borate(h, tb, kb)         &
                + alk_water(h, kw, ts, ks)      &
                + alk_sulfate(h, ts, ks)        &
                + alk_fluoride(h, tf, kf)       &
                + alk_phosphate(h, po4, kp1, kp2, kp3) &
                + alk_silicate(h, si, ksi)      &
                + alk_ammonia(h, nh4, knh4)     &
                + alk_sulfide(h, h2s, kh2s)

   end function alk_total

   real(rk) function dalk_total_dh(h, dic, po4, si, nh4, h2s, tb, ts, tf, &
                                   k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s)
      real(rk), intent(in) :: h, dic, po4, si, nh4, h2s
      real(rk), intent(in) :: tb, ts, tf
      real(rk), intent(in) :: k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s

      dalk_total_dh = dalk_carbonate_dh(h, dic, k1, k2) &
                    + dalk_borate_dh(h, tb, kb)         &
                    + dalk_water_dh(h, kw, ts, ks)      &
                    + dalk_sulfate_dh(h, ts, ks)        &
                    + dalk_fluoride_dh(h, tf, kf)       &
                    + dalk_phosphate_dh(h, po4, kp1, kp2, kp3) &
                    + dalk_silicate_dh(h, si, ksi)      &
                    + dalk_ammonia_dh(h, nh4, knh4)     &
                    + dalk_sulfide_dh(h, h2s, kh2s)

   end function dalk_total_dh   

   ! ------------------------------------------------------------------------------------------
   ! Alkalinity contributions expressed as functions of H+ on the total pH scale.
   ! except Ks, which is retained on the free scale and used to convert
   ! between free and total proton conventions.
   ! ------------------------------------------------------------------------------------------

   ! Carbonate alkalinity from DIC speciation: 
   ! Derived from equilibrium partitioning of CO2, HCO3-, CO3-- using K1, K2.
   ! CO2 ​+ H2​O ↔ H+ + HCO3−​ and HCO3− ​↔ H+ + CO32−​
   ! Denominator is total DIC normalization; numerator is charge-weighted sum (HCO3- + 2*CO3--).
   pure real(rk) function alk_carbonate(h, dic, k1, k2)
      real(rk), intent(in) :: h, dic, k1, k2
      real(rk) :: denom

      denom = h*h + k1*h + k1*k2
      alk_carbonate = dic * k1 * (h + 2.0_rk*k2) / denom
   end function alk_carbonate

   ! d(Alk_carbonate)/dh from analytical differentiation of DIC speciation using quotient rule.
   pure real(rk) function dalk_carbonate_dh(h, dic, k1, k2)
      real(rk), intent(in) :: h, dic, k1, k2
      real(rk) :: denom, numer

      denom = h*h + k1*h + k1*k2
      numer = k1*k2*k1 + h*(4.0_rk*k1*k2 + h*k1)
      dalk_carbonate_dh = -dic * numer / (denom*denom)
   end function dalk_carbonate_dh

   ! Borate alkalinity from B(OH)3 ↔ B(OH)4- equilibrium using KB.
   ! Only B(OH)4- contributes (charged species), so Alk = [B(OH)4-].
   pure real(rk) function alk_borate(h, tb, kb)
      real(rk), intent(in) :: h, tb, kb

      alk_borate = tb * kb / (kb + h)
   end function alk_borate

   ! d(Alk_borate)/dh from analytical differentiation of borate equilibrium expression.
   pure real(rk) function dalk_borate_dh(h, tb, kb)
      real(rk), intent(in) :: h, tb, kb

      dalk_borate_dh = -tb * kb / (kb + h)**2
   end function dalk_borate_dh

   ! Water alkalinity from OH- minus free [H+], with total-to-free proton conversion via sulfate.
   ! Uses KW/h for OH- and h/(1 + TS/KS) for free proton contribution.
   pure real(rk) function alk_water(h, kw, ts, ks)
      real(rk), intent(in) :: h, kw, ts, ks
      real(rk) :: free_to_tot

      free_to_tot = 1.0_rk + ts / ks
      alk_water = kw / h - h / free_to_tot
   end function alk_water

   ! d(Alk_water)/dh from analytical differentiation of water equilibrium expression.
   pure real(rk) function dalk_water_dh(h, kw, ts, ks)
      real(rk), intent(in) :: h, kw, ts, ks
      real(rk) :: free_to_tot

      free_to_tot = 1.0_rk + ts / ks
      dalk_water_dh = -kw / (h*h) - 1.0_rk / free_to_tot
   end function dalk_water_dh


   ! Sulfate alkalinity from HSO4- ↔ SO4-- equilibrium, referenced to zero level HSO4-.
   ! Total-to-free proton conversion enters through (1 + TS/KS) because h is on the total scale.
   pure real(rk) function alk_sulfate(h, ts, ks)
      real(rk), intent(in) :: h, ts, ks
      real(rk) :: free_to_tot, denom

      free_to_tot = 1.0_rk + ts / ks
      denom = ks * free_to_tot + h
      alk_sulfate = ts * (ks * free_to_tot / denom - 1.0_rk)
   end function alk_sulfate

   ! d(Alk_sulfate)/dh from analytical differentiation of sulfate equilibrium expression.
   pure real(rk) function dalk_sulfate_dh(h, ts, ks)
      real(rk), intent(in) :: h, ts, ks
      real(rk) :: free_to_tot, denom

      free_to_tot = 1.0_rk + ts / ks
      denom = ks * free_to_tot + h
      ! Double check this derivative and if ks shoudl be converted to the total scale
      dalk_sulfate_dh = -ts * ks * free_to_tot / (denom*denom)
   end function dalk_sulfate_dh

   ! Fluoride alkalinity from HF ↔ F- equilibrium, referenced to fully protonated HF.
   ! Alk = [F-] - TF = -[HF], so contribution is always ≤ 0.   
   pure real(rk) function alk_fluoride(h, tf, kf)
      real(rk), intent(in) :: h, tf, kf

      alk_fluoride = tf * (kf / (kf + h) - 1.0_rk)
   end function alk_fluoride

   ! d(Alk_fluoride)/dh from analytical differentiation of fluoride equilibrium expression.
   pure real(rk) function dalk_fluoride_dh(h, tf, kf)
      real(rk), intent(in) :: h, tf, kf

      dalk_fluoride_dh = -tf * kf / (kf + h)**2
   end function dalk_fluoride_dh

   ! Phosphate alkalinity from H3PO4/H2PO4-/HPO4--/PO4--- speciation using KP1, KP2, KP3.
   ! Zero level referenced to fully protonated H3PO4; numerator is the charge-weighted species sum.
   pure real(rk) function alk_phosphate(h, po4, kp1, kp2, kp3)
      real(rk), intent(in) :: h, po4, kp1, kp2, kp3
      real(rk) :: numer, denom

      if (po4 <= 0.0_rk) then      ! If PO4 is not present or is 0
         alk_phosphate = 0.0_rk
         return
      end if

      numer = 3.0_rk*kp1*kp2*kp3 + h*(2.0_rk*kp1*kp2 + h*kp1)
      denom = kp1*kp2*kp3 + h*(kp1*kp2 + h*(kp1 + h))

      alk_phosphate = po4 * (numer / denom - 1.0_rk)
   end function alk_phosphate

   ! d(Alk_phosphate)/dh from analytical differentiation of phosphate speciation using quotient rule.
   pure real(rk) function dalk_phosphate_dh(h, po4, kp1, kp2, kp3)
      real(rk), intent(in) :: h, po4, kp1, kp2, kp3
      real(rk) :: denom, numer

      if (po4 <= 0.0_rk) then
         dalk_phosphate_dh = 0.0_rk
         return
      end if

      denom = kp1*kp2*kp3 + h*(kp1*kp2 + h*(kp1 + h))
      numer = kp1*kp2*kp1*kp2*kp3 + h*(4.0_rk*kp1*kp1*kp2*kp3 + &
              h*(9.0_rk*kp1*kp2*kp3 + kp1*kp1*kp2 + h*(4.0_rk*kp1*kp2 + h*kp1)))
      dalk_phosphate_dh = -po4 * numer / (denom*denom)
   end function dalk_phosphate_dh

   ! Silicate alkalinity from H4SiO4 ↔ H3SiO4- equilibrium using KSI.
   ! Only the deprotonated species contributes, so Alk = [H3SiO4-].
   pure real(rk) function alk_silicate(h, si, ksi)
      real(rk), intent(in) :: h, si, ksi

      if (si <= 0.0_rk) then
         alk_silicate = 0.0_rk
         return
      end if

      alk_silicate = si * ksi / (ksi + h)
   end function alk_silicate

   ! d(Alk_silicate)/dh from analytical differentiation of silicate equilibrium expression.
   pure real(rk) function dalk_silicate_dh(h, si, ksi)
      real(rk), intent(in) :: h, si, ksi

      if (si <= 0.0_rk) then
         dalk_silicate_dh = 0.0_rk
         return
      end if

      dalk_silicate_dh = -si * ksi / (ksi + h)**2
   end function dalk_silicate_dh

   ! Ammonia alkalinity from NH4+ ↔ NH3 equilibrium using KNH4.
   ! Assumes nh4 is total ammonia; only the base form NH3 contributes to alkalinity.
   pure real(rk) function alk_ammonia(h, nh4, knh4)
      real(rk), intent(in) :: h, nh4, knh4

      if (nh4 <= 0.0_rk) then
         alk_ammonia = 0.0_rk
         return
      end if

      alk_ammonia = nh4 * knh4 / (knh4 + h)
   end function alk_ammonia

   ! d(Alk_ammonia)/dh from analytical differentiation of ammonia equilibrium expression.
   pure real(rk) function dalk_ammonia_dh(h, nh4, knh4)
      real(rk), intent(in) :: h, nh4, knh4

      if (nh4 <= 0.0_rk) then
         dalk_ammonia_dh = 0.0_rk
         return
      end if

      dalk_ammonia_dh = -nh4 * knh4 / (knh4 + h)**2
   end function dalk_ammonia_dh

   ! Sulfide alkalinity from H2S ↔ HS- equilibrium using KH2S.
   ! Assumes h2s is total sulfide; only the deprotonated species HS- is included here.
   pure real(rk) function alk_sulfide(h, h2s, kh2s)
      real(rk), intent(in) :: h, h2s, kh2s

      if (h2s <= 0.0_rk) then
         alk_sulfide = 0.0_rk
         return
      end if

      alk_sulfide = h2s * kh2s / (kh2s + h)
   end function alk_sulfide

   ! d(Alk_sulfide)/dh from analytical differentiation of sulfide equilibrium expression.
   pure real(rk) function dalk_sulfide_dh(h, h2s, kh2s)
      real(rk), intent(in) :: h, h2s, kh2s

      if (h2s <= 0.0_rk) then
         dalk_sulfide_dh = 0.0_rk
         return
      end if

      dalk_sulfide_dh = -h2s * kh2s / (kh2s + h)**2
   end function dalk_sulfide_dh   

   !--------------------------------------------------------------
   ! Testing subroutine to compare against the results of pyCO2sys
   !--------------------------------------------------------------
   subroutine test_carbonate_solver()
      integer, parameter :: ncases = 5
      integer :: i

      character(len=24) :: names(ncases)

      ! ---------------- Input test cases ----------------
      ! State/environment values in model-like units
      real(rk) :: temp(ncases), sal(ncases), pres(ncases), rho(ncases)
      real(rk) :: dic(ncases), alk(ncases), po4(ncases), si(ncases), nh4(ncases), h2s(ncases)

      ! ---------------- Converted concentrations ----------------
      real(rk) :: dic_kg, alk_kg, po4_kg, si_kg, nh4_kg, h2s_kg

      ! ---------------- Background seawater totals ----------------
      real(rk) :: tb, ts, tf

      ! ---------------- Equilibrium constants ----------------
      real(rk) :: k1, k2, kb, kw, ks, kf
      real(rk) :: kp1, kp2, kp3, ksi, knh4, kh2s
      real(rk) :: ksp_ca, ca

      ! ---------------- Solver/speciation outputs ----------------
      real(rk) :: h, ph
      real(rk) :: co2_kg, hco3_kg, co3_kg
      real(rk) :: co2_out, hco3_out, co3_out
      real(rk) :: omega_ca
      real(rk) :: denom

      ! ---------------- Initial guess ----------------
      real(rk) :: h_init

      ! ------------------------------------------------------------------
      ! Five test cases
      ! Units:
      !   temp : degC
      !   sal  : PSU
      !   pres : dbar
      !   rho  : kg m-3
      !   tracers : mmol m-3
      ! ------------------------------------------------------------------
      names = [character(len=20) :: 'surface_ref', 'deep_ref', 'reduced_case', &
                                    'ammonia_case', 'sulfide_case']

      temp = [15.0_rk,   4.0_rk,   10.0_rk, 10.0_rk, 10.0_rk]
      sal  = [35.0_rk,  35.0_rk,   35.0_rk, 35.0_rk, 35.0_rk]
      pres = [10.0_rk, 2000.0_rk, 100.0_rk, 100.0_rk, 100.0_rk]
      rho  = [1025.0_rk, 1040.0_rk, 1028.0_rk, 1028.0_rk, 1028.0_rk]

      dic  = [2100.0_rk, 2200.0_rk, 2400.0_rk, 2300.0_rk, 2300.0_rk]
      alk  = [2300.0_rk, 2350.0_rk, 2500.0_rk, 2400.0_rk, 2400.0_rk]
      po4  = [0.0_rk,      1.5_rk,    2.0_rk,    0.0_rk,    0.0_rk]
      si   = [0.0_rk,     20.0_rk,   15.0_rk,    0.0_rk,    0.0_rk]
      nh4  = [0.0_rk,      0.0_rk,   10.0_rk,   20.0_rk,    0.0_rk]
      h2s  = [0.0_rk,      0.0_rk,    5.0_rk,    0.0_rk,   20.0_rk]

      write(*,*)
      write(*,'(A)') '==========================='
      write(*,'(A)') 'Testing carbonate_system'
      write(*,'(A)') '==========================='

      do i = 1, ncases

         ! Convert from mmol m-3 to mol kg-1
         dic_kg = dic(i) / (1000.0_rk * rho(i))
         alk_kg = alk(i) / (1000.0_rk * rho(i))
         po4_kg = po4(i) / (1000.0_rk * rho(i))
         si_kg  = si(i)  / (1000.0_rk * rho(i))
         nh4_kg = nh4(i) / (1000.0_rk * rho(i))
         h2s_kg = h2s(i) / (1000.0_rk * rho(i))

         ! Background totals from salinity
         call compute_background_totals_from_salinity(sal(i), tb, ts, tf)

         ! Equilibrium constants
         call compute_equilibrium_constants(temp(i), sal(i), pres(i), ts, tf, &
                                            k1, k2, kb, kw, ks, kf, &
                                            kp1, kp2, kp3, ksi, knh4, kh2s, ksp_ca)

         ! Cold-start initial guess
         h_init = -1.0_rk

         ! Solve carbonate system
         h = solve_hplus(dic_kg, alk_kg, po4_kg, si_kg, nh4_kg, h2s_kg, &
                         tb, ts, tf, &
                         k1, k2, kb, kw, ks, kf, kp1, kp2, kp3, ksi, knh4, kh2s, h_init)

         ph = -log10(max(h, eps_h))

         ! Carbonate speciation
         denom   = h*h + k1*h + k1*k2
         co2_kg  = dic_kg * h*h   / denom
         hco3_kg = dic_kg * k1*h  / denom
         co3_kg  = dic_kg * k1*k2 / denom

         ! Convert back to mmol m-3 for easier comparison with model diagnostics
         co2_out  = 1000.0_rk * rho(i) * co2_kg
         hco3_out = 1000.0_rk * rho(i) * hco3_kg
         co3_out  = 1000.0_rk * rho(i) * co3_kg

         ! Calcite saturation state
         ca = calcon_ref * sal(i) / 35.0_rk
         omega_ca = ca * co3_kg / (ksp_ca + eps_h)

         write(*,*)
         write(*,'(A,A)') 'CASE: ', trim(names(i))
         write(*,'(A,F10.4,2X,A,F10.4,2X,A,F10.4,2X,A,F10.4)') &
            'T=', temp(i), 'S=', sal(i), 'P=', pres(i), 'rho=', rho(i)
         write(*,'(A,F12.5,2X,A,F12.5,2X,A,F12.5,2X,A,F12.5,2X,A,F12.5,2X,A,F12.5)') &
            'DIC=', dic(i), 'ALK=', alk(i), 'PO4=', po4(i), 'Si=', si(i), 'NH4=', nh4(i), 'H2S=', h2s(i)

         write(*,'(A,ES16.8,2X,A,ES16.8,2X,A,ES16.8)') &
            'TB=', tb, 'TS=', ts, 'TF=', tf

         write(*,'(A,ES16.8,2X,A,ES16.8,2X,A,ES16.8,2X,A,ES16.8)') &
            'K1=', k1, 'K2=', k2, 'KB=', kb, 'KW=', kw
         write(*,'(A,ES16.8,2X,A,ES16.8,2X,A,ES16.8,2X,A,ES16.8)') &
            'KS=', ks, 'KF=', kf, 'KP1=', kp1, 'KP2=', kp2
         write(*,'(A,ES16.8,2X,A,ES16.8,2X,A,ES16.8,2X,A,ES16.8)') &
            'KP3=', kp3, 'KSi=', ksi, 'KNH4=', knh4, 'KH2S=', kh2s
         write(*,'(A,ES16.8)') 'Ksp_calcite=', ksp_ca

         write(*,'(A,ES16.8,2X,A,F12.6)') 'H=', h, 'pH=', ph
         write(*,'(A,F12.6,2X,A,F12.6,2X,A,F12.6,2X,A,F12.6)') &
            'CO2*=', co2_out, 'HCO3-=', hco3_out, 'CO3--=', co3_out, 'OmegaCa=', omega_ca

      end do

      write(*,*)
      write(*,'(A)') 'End carbonate_system test'
      write(*,'(A)') '=============================================================='
      write(*,*)

   end subroutine test_carbonate_solver

end module carbonate_system