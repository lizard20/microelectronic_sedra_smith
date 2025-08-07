'''
    constants.py

    Constants for use of the semiconductor library

    Author: Aldo Nunez Tovar
    Date: May-2025
'''

from astropy import constants as const
from astropy import units as u

#################### Main constants ########################

atoms_si = 5e22 * u.dimensionless_unscaled  # 5x10²²Number of atoms in Si
T = 300 * u.Kelvin                          # [°K]. Room temperature in °K. Or  27 °C
B  = 7.3e15 * 1/(u.cm**3 * u.Kelvin**1.5)   # 7.3x10¹⁵ [1/(cm^3 K^3/2)]. Material dependent parameter. B for intrinsic Si
B_gaas = 3.56e14 * 1/(u.cm**3 * u.Kelvin**1.5)  # 3.5x10¹⁴ [1/(cm^3 K^3/2)]. Material dependent parameter. B_gaas for GaAs
E_g = 1.12 * u.eV        # [eV]. Band energy separation for Si
E_g_gaas = 1.42 * u.eV   # [eV]. Band energy separation for GaAs
k_J_K = const.k_B       # 1.380649×10 −23 [J/K]
k = k_J_K.to('eV/K')    # 8.62x10⁻⁵ [ev/°K] - Boltzmann's constant
q = const.e.to(u.C)     # 1.60217663x10⁻¹⁹ [Coulomb] Electron charge
eps_0 = const.eps0.to(u.F/u.cm)     #  [F/cm] - Vacuum permitivity
eps_si = 11.7 * eps_0               # [F/cm] - Silicon permitivity
mu_p = 480 * u.cm**2/(u.Volt * u.second)        # [cm²/Vs]. Holes mobility in intrinsic Si
mu_n = 1_350 * u.cm**2/(u.Volt * u.second)      # [cm²/Vs]. Electrons mobility intrinsic Si
mu_p_p = 400 * u.cm**2/(u.Volt * u.second)      # [cm²/Vs].  Holes mobility in doped Si
mu_n_p = 1_110 * u.cm**2/(u.Volt * u.second)    # [cm²/Vs]. Electrons mobility doped Si
n = 1 * u.dimensionless_unscaled    # Unitless parameter
#n = 2 * u.dimensionless_unscaled
N_A = 1e16 * 1/u.cm**3          # 1x10¹⁶ [carriers/cm³]. Acceptor atoms density
N_D = 1e17 *1/u.cm**3           # 1x10¹⁷ [carriers/cm³]. Donors atoms density
ni_si = 1.5e10 * 1/u.cm**3      #  1.5x10¹⁰ [carriers/cm³].  Density of charge carriers for intrinsic Si at 300 °K
D_n = 35 * u.cm**2 / u.second   # [cm²/s]. Diffusivity of electtrons in intrinsic  Si
D_p = 12 * u.cm**2 / u.second   # [cm²/s]. Diffusivity of holes in intrinsic  Si




