# Constants

# c.f. https://gitlab.dkrz.de/icon/icon-mpim/-/blob/master/src/shared/mo_physical_constants.f90
grav   = 9.80665
ckap   = 0.4
rv     = 461.51
rd     = 287.04
vtmpc1 = rv/(rd-1.0)
rgrav  = 1.0/grav

# c.f. https://gitlab.dkrz.de/icon/icon-mpim/-/blob/master/src/atm_phy_aes/tmx/mo_vdf_diag_smag.f90?ref_type=heads#L57
bsm = 5.0  #Businger Stable Momentum
bum = 16.0   #Businger Untable Momentum
bsh = 5.0  #Businger Stable Heat
buh = 16.0   #Businger Untable Heat

# c.f. https://gitlab.dkrz.de/icon/icon-mpim/-/blob/master/src/shared/mo_math_constants.f90
pi_2 = 1.57079632679489661923132169163975144
ln2 = 0.693147180559945309417232121458176568