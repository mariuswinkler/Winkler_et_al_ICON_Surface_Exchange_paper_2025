#%%
# Functions

import numpy as np
import src.constants as c

def stability_function_mom(RIB, hz0, tc):
    if RIB >= 0:
        stab_fun = 1 / (1 + 10 * RIB * (1 + 8 * RIB))
    else:
        hz0_fac = (max(hz0, 1)**(1/3) - 1)**1.5 
        stab_fun = 1 + 10 * abs(RIB) / (1 + 75 * tc * hz0_fac * np.sqrt(abs(RIB)))
    return stab_fun

def stability_function_heat(RIB, hzh, tc):
    if RIB >= 0:
        stab_fun = 1 / (1 + 10 * RIB * (1 + 8 * RIB))
    else:
        hzh_fac = (max(hzh, 1)**(1/3) - 1)**1.5 
        stab_fun = 1 + 15 * abs(RIB) / (1 + 75 * tc * hzh_fac * np.sqrt(abs(RIB)))
    return stab_fun

def businger_heat(z0, z1, L):
    if L > 0:  # Stable
        zeta = z1 / L
        zeta0 = z0 / L
        if zeta > 1: 
            psi = -c.bsh * np.log(zeta) - zeta + 1
            psi0 = -c.bsh * np.log(zeta0) - zeta0 + 1
            factor = (np.log(L / z0) + c.bsh - psi + psi0) / c.ckap
        else:
            psi = -c.bsh * zeta
            psi0 = -c.bsh * zeta0
            factor = (np.log(z1 / z0) - psi + psi0) / c.ckap
    elif L < 0:  # Unstable
        zeta = z1 / L
        zeta0 = z0 / L
        lamda = np.sqrt(1 - c.buh * zeta)
        lamda0 = np.sqrt(1 - c.buh * zeta0)
        psi = 2 * (np.log(1 + lamda) - c.ln2)
        psi0 = 2 * (np.log(1 + lamda0) - c.ln2)
        factor = (np.log(z1 / z0) - psi + psi0) / c.ckap
    else:  # Neutral
        factor = np.log(z1 / z0) / c.ckap
    return factor

def businger_mom(z0, z1, L):
    if L > 0:  # Stable
        zeta = z1 / L
        zeta0 = z0 / L
        if zeta > 1: 
            psi = -c.bsm * np.log(zeta) - zeta + 1
            psi0 = -c.bsm * np.log(zeta0) - zeta0 + 1
            factor = (np.log(L / z0) + c.bsh - psi + psi0) / c.ckap
        else:
            psi = -c.bsm * zeta
            psi0 = -c.bsm * zeta0
            factor = (np.log(z1 / z0) - psi + psi0) / c.ckap
    elif L < 0:  # Unstable
        zeta = z1 / L
        zeta0 = z0 / L
        lamda = np.sqrt(np.sqrt(1 - c.bum * zeta))
        lamda0 = np.sqrt(np.sqrt(1 - c.bum * zeta0))
        psi = 2 * np.log(1 + lamda) + np.log(1 + lamda**2) - \
              2 * np.arctan(lamda) + c.pi_2 - 3 * c.ln2
        psi0 = 2 * np.log(1 + lamda0) + np.log(1 + lamda0**2) - \
               2 * np.arctan(lamda0) + c.pi_2 - 3 * c.ln2
        factor = (np.log(z1 / z0) - psi + psi0) / c.ckap
    else:  # Neutral
        factor = np.log(z1 / z0) / c.ckap
    return factor

def sfc_exchange_coefficients(dz, pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, min_wind_threshold=1.0):
    zepsec = 0.028
    zcons17 = 1.0 / c.ckap**2

    mwind = max(mwind, min_wind_threshold)

    z_mc = dz
    # First guess for tch and tcm using bulk approach
    RIB = c.grav * (thetam1 - theta_sfc) * (z_mc - rough_m) / (theta_sfc * mwind**2)
    tcn_mom = (c.ckap / np.log(z_mc / rough_m))**2
    tcm = tcn_mom * stability_function_mom(RIB, z_mc / rough_m, tcn_mom)
    stab_func_mom_out = stability_function_mom(RIB, z_mc / rough_m, tcn_mom)
    
    tcn_heat = c.ckap**2 / (np.log(z_mc / rough_m) * np.log(z_mc / rough_m))
    tch = tcn_heat * stability_function_heat(RIB, z_mc / rough_m, tcn_heat)
    stab_funkheat_out = stability_function_heat(RIB, z_mc / rough_m, tcn_heat)
    
    # Now iterate
    for itr in range(5):
        shfl_local = tch * mwind * (theta_sfc - thetam1)
        lhfl_local = tch * mwind * (qsat_sfc - pqm1)
        bflx1 = shfl_local + (c.vtmpc1 * theta_sfc * lhfl_local)
        ustar = np.sqrt(tcm) * mwind

        obukhov_length = -ustar**3.0 * theta_sfc * c.rgrav / (c.ckap * bflx1)

        inv_bus_mom = 1.0 / businger_mom(rough_m, z_mc, obukhov_length)
        tch = inv_bus_mom / businger_heat(rough_m, z_mc, obukhov_length)
        tcm = inv_bus_mom**2.0

    cH = tch
    cD = tcm
    cH_neutral = c.ckap / max(zepsec, np.sqrt(tcn_heat))
    cD_neutral = c.ckap / max(zepsec, np.sqrt(tcn_mom))

    return cD, cH, cD_neutral, cH_neutral, RIB, mwind, stab_func_mom_out, stab_funkheat_out