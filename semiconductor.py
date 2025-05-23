# -*- coding: utf-8 -*-
'''
    semiconductor.py

    This library translates the main equations 
    from chapter 3 of the Sedra's book into the Python 
    programming language as an aid to solving problems.

    References:
    The equations were taken from  chapter 3 
    of the book:  Microelectronic Circuits
    Sedra and Smith
    7th. edition

    Version: 0.1

    Author: Aldo Nunez Tovar
    Date: May-2025
'''

import math as m
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

################################################

################# Methods #########################
def n_i(T: u.Quantity = T) -> u.Quantity:
    '''
        Parameters:
            T: Temperature. Expected units: Kelvin (K).
               Default value: 300 °K
        Return:
            n_i: Intrinsic carrier density. Expected units: 1/cm**3.

        Calculate the intrinsic carriers density in intrinsic Si
        n_i = B * (T**1.5) * exp(-Eg/2*k*T)    (Eq. 3.2)
    '''
    # Validate input parameters
    if not T.unit.is_equivalent(u.Kelvin):
        raise u.UnitsError("Input temperature T must be in units equivalent to Kelvin.")

    exponent = -E_g / (2.0 * k * T)
    return  B * (T**1.5) * m.exp(exponent)

def n_i_gaas(T: u.Quantity = T) -> u.Quantity:
    '''
        Parameters:
            T: Temperature. Expected units: Kelvin (K).
               Default value: 300 °K
        Return:
            n_i: Intrinsic carrier density in Galium Arsenide. Expected units: 1/cm**3.

        Calculate the intrinsic carriers density in GaAs
        ni_gaas = B_gaas*(T**1.5)*e(-Eg_gaas/2*k*T)
    '''
    # Validate input parameters
    if not T.unit.is_equivalent(u.Kelvin):
         raise u.UnitsError("Input temperature T must be in units equivalent to Kelvin.")

    exponent = -E_g_gaas / (2.0 * k * T)
    return  B_gaas * (T**1.5) * m.exp(exponent)

def p_n(T: u.Quantity = T,
        N_D: u.Quantity = N_D) -> u.Quantity:
    '''
        Parameters:
            T: Temperature. Expected units: Kelvin (K).
               Default value: 300 °K
            N_D: Donor atoms concentration. Expected units: atoms/cm³.
                 Default value: 1x10¹⁷ atoms/cm³
        Return:
            p_n: Hole concentration. Expected units: 1/cm**3.

        Concentration of holes, minority charge carriers in  N-type material:
        p_n = ni² / N_D    (Eq. 3.5)
    '''
    # Validate input parameters
    if not T.unit.is_equivalent(u.Kelvin):
         raise u.UnitsError("Input temperature T must be in units equivalent to Kelvin.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
         raise u.UnitsError("Input atoms density N_D must be in units equivalent to 1/cm**3.")

    return  (n_i(T))**2 / N_D

def n_p(T: u.Quantity = T,
        N_A: u.Quantity = N_A) -> u.Quantity:
    '''
        Parameters:
            T: Temperature. Expected units: Kelvin (K).
               Default value: 300 °K
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
                 Default value: 1x10¹⁷ atoms/cm³
        Return:
            n_p: electron concentration. Expected units: 1/cm**3.

        Concentration of electrons, minority charge carriers in  P-type material:
        n_p = ni² / N_A    (Eq. 3.7)
    '''
    # Validate input parameters
    if not T.unit.is_equivalent(u.Kelvin):
         raise u.UnitsError("Input temperature T must be in units equivalent to Kelvin.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
         raise u.UnitsError("Input atoms density N_A must be in units equivalent to 1/cm**3.")

    return  (n_i(T))**2 / N_A


def sigma(n: u.Quantity = ni_si,
        p: u.Quantity = ni_si, 
        mu_n: u.Quantity = mu_n, 
        mu_p: u.Quantity = mu_p ) -> u.Quantity:
    '''
    Parameters:
      n: Electron concentration. Expected units: 1/cm**3.
         Default value:  1.5x10¹⁰ atoms/cm³
      p: Hole concentration. Expected units: 1/cm**3.
         Default value:  1.5x10¹⁰ atoms/cm³
      mu_n: Electron mobility. Expected units: cm²/Vs.
            Default value: 1350 cm²/Vs
      mu_p: Hole mobility. Expected units: cm²/Vs.
            Default value: 480 cm²/Vs
    Return:
      sigma: conductivity: siemens/cm

    Calculate the conductivity in semiconductor material:
        sigma =  q * (n * mu_n + p * mu_p) (Eq. 3.16)
    '''
    # Validate input parameters
    if not n.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input electron concentration n must be in units equivalent to 1/cm**3.")
    if not p.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input hole concentration p must be in units equivalent to 1/cm**3.")
    if not mu_n.unit.is_equivalent(u.cm**2/(u.Volt * u.second)):
        raise u.UnitsError("Input electron mobility mu_n must be in units equivalent to cm²/Vs.")
    if not mu_p.unit.is_equivalent(u.cm**2/(u.Volt * u.second)):
        raise u.UnitsError("Input hole mobility mu_p must be in units equivalent to cm²/Vs.")

    return  q * (n * mu_n + p * mu_p)


def rho(n: u.Quantity = ni_si,
        p: u.Quantity = ni_si, 
        mu_n: u.Quantity = mu_n, 
        mu_p: u.Quantity = mu_p ) -> u.Quantity:
    '''
    Parameters:
      n: Electron concentration. Expected units: 1/cm**3.
         Default value:  1.5x10¹⁰ atoms/cm³
      p: Hole concentration. Expected units: 1/cm**3.
         Default value:  1.5x10¹⁰ atoms/cm³
      mu_n: Electron mobility. Expected units: cm²/Vs.
            Default value: 1350 cm²/Vs
      mu_p: Hole mobility. Expected units: cm²/Vs.
            Default value: 480 cm²/Vs
    Return:
      rho: Resistivity: ohms-cm

    Calculate the resistivity in semiconductor material:
        rho = 1 / (q * (n * mu_n + p * mu_p)) (Eq. 3.17)
    '''
    # Validate input parameters
    if not n.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input electron concentration n must be in units equivalent to 1/cm**3.")
    if not p.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input hole concentration p must be in units equivalent to 1/cm**3.")
    if not mu_n.unit.is_equivalent(u.cm**2/(u.Volt * u.second)):
        raise u.UnitsError("Input electron mobility mu_n must be in units equivalent to cm²/Vs.")
    if not mu_p.unit.is_equivalent(u.cm**2/(u.Volt * u.second)):
        raise u.UnitsError("Input hole mobility mu_p must be in units equivalent to cm²/Vs.")

    return  (1/(q * (n * mu_n + p * mu_p)))

def vp_drift(E: u.Quantity, 
             mu_p: u.Quantity =  mu_n) -> u.Quantity:
    '''
        Parameters:
            E: Electric field. Expected units: V/cm.
            mu_p: Holes mobility. Expected units: cm²/Vs.
                  Default value: 480 cm²/Vs
        Return:
            vp_drift: Drift velocity of holes. Expected units: cm/s.
        
        Calculate the holes drift velocity:
        vp_drift = -mu_p * E   (3.9)
  '''
    
    # Validate input parameters
    if not E.unit.is_equivalent(u.Volt/u.cm):
        raise u.UnitsError("Input electric field E must be in units equivalent to V/cm.")
    if not mu_p.unit.is_equivalent(u.cm**2/(u.Volt * u.second)):
        raise u.UnitsError("Input electron mobility mu_p must be in units equivalent to cm²/Vs.")
    
    return mu_p * E

def vn_drift(E: u.Quantity, 
             mu_n: u.Quantity =  mu_n) -> u.Quantity:
    '''
        Parameters:
            E: Electric field. Expected units: V/cm.
            mu_n: Electron mobility. Expected units: cm²/Vs.
                  Default value: 1350 cm²/Vs
        Return:
            vn_drift: Drift velocity of electrons. Expected units: cm/s.
        
        Calculate the electron drift velocity:
        vn_drift = -mu_n * E   (3.9)
  '''
    
    # Validate input parameters
    if not E.unit.is_equivalent(u.Volt/u.cm):
        raise u.UnitsError("Input electric field E must be in units equivalent to V/cm.")
    if not mu_n.unit.is_equivalent(u.cm**2/(u.Volt * u.second)):
        raise u.UnitsError("Input electron mobility mu_n must be in units equivalent to cm²/Vs.")
    
    return -mu_n * E

def J_sigma(E: u.Quantity, 
      sigma: u.Quantity ) -> u.Quantity:
    '''
        Parameters:
            E: Electric field. Expected units: V/cm.
            sigma: conductivity. Expected units: siemens/cm.
        Return:
            J: Drift current density. Expected units: A/cm²
    
        Calculate the drift current density from conductivity, sigma
        J = sigma * E   (3.14)
    ''' 
    # Validate input parameters
    if not E.unit.is_equivalent(u.Volt/u.cm):
        raise u.UnitsError("Input electric field E must be in units equivalent to V/cm.")
    if not sigma.unit.is_equivalent(u.siemens / u.cm):
        raise u.UnitsError("Input resistivity rho must be in units equivalent to mho/cm.")
    
    return (sigma * E).to(u.A * 1/u.cm**2)



def J_rho(E: u.Quantity, 
      rho: u.Quantity ) -> u.Quantity:
    '''
        Parameters:
            E: Electric field. Expected units: V/cm.
            rho: Resistivity. Expected units: ohms-cm.
        Return:
            J: Drift current density. Expected units: A/cm²
    
        Calculate the drift current density from ressitivity, rho
        J = E / rho   (3.15)
    ''' 
    # Validate input parameters
    if not E.unit.is_equivalent(u.Volt/u.cm):
        raise u.UnitsError("Input electric field E must be in units equivalent to V/cm.")
    if not rho.unit.is_equivalent(u.Ohm * u.cm):
        raise u.UnitsError("Input resistivity rho must be in units equivalent to ohms-cm.")
    
    return (E / rho).to(u.A * 1/u.cm**2)


def J_p(D_p: u.Quantity,
        dp_dx: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            D_p: diffusion constant of holes. Expected unit: cm²/s
            dp_dx: concentration gradient of holes. Expected units:  1/cm⁴
        Return:
            J_p, holes current density. Expected units A/cm²

        Calculate the holes current density
        J_p = -q*D_p*dp(x)/dx    (3.19)
    '''
    
    # Validate input parameters
    if not D_p.unit.is_equivalent(u.cm**2/u.second):
        raise u.UnitsError("Input diffusion constant D_p must be in units equivalent to cm²/s.")
    if not dp_dx.unit.is_equivalent(1/u.cm**4):
        raise u.UnitsError("Input concentration gradient dp_dx must be in units equivalent to 1/cm⁴.")

    return q * D_p * dp_dx

def J_n(D_n: u.Quantity,
        dn_dx: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            D_n: diffusion constant of electrons. Expected unit: cm²/s
            dn_dx: concentration gradient of electrons. Expected units:  1/cm⁴
        Return:
            J_n, electrons current density. Expected units A/cm²

        Calculate the electons current density
        J_n = -q * D_n* dn(x)/dx    (3.19)
    '''
    
    # Validate input parameters
    if not D_n.unit.is_equivalent(u.cm**2/u.second):
        raise u.UnitsError("Input diffusion constant D_n must be in units equivalent to cm²/s.")
    if not dn_dx.unit.is_equivalent(1/u.cm**4):
        raise u.UnitsError("Input concentration gradient dn_dx must be in units equivalent to 1/cm⁴.")

    return q * D_n * dn_dx


def V_O(N_A: u.Quantity = N_A,
       N_D: u.Quantity  = N_D) -> u.Quantity:
    '''
        Parameters:
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³. 
                Default value: N_A = 1e16 atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³.
                Default value: N_D = 1e17 atoms/cm³
        Returns:
            V_O, built-in voltage. Expected units: Volts (V)

        Calculate the built-in Voltage:
        V_O = VT * ln(N_A * N_D/ni²)   (3.22)
  '''
    # Validate input parameters
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")

    return  V_T() * m.log(N_A * N_D / n_i()**2)


def Q_n(N_D: u.Quantity,
        x_n: u.Quantity,
        A: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            N_D: Donor atoms concentration. Expected units: atoms/cm³
            x_n: Width of the depletion layer in the N region. Expected units: cm
            A: Area of the depletion region. Expected units: cm²
        Return:
            Q_n: Charge  stored. Expected units: Coulomb (C)

        Calculate the charge stored on the N  side
        |Q+| = q * A * x_n * N_D                  (3.23)
    '''
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")
    if not x_n.unit.is_equivalent(u.cm):
        raise u.UnitsError("Input width of the depletion region x_n must be in units equivalent to cm.")
    if not A.unit.is_equivalent(u.cm**2):
        raise u.UnitsError("Input area A must be in units equivalent to cm².")

    return q * A * x_n * N_D

def Q_p(N_A: u.Quantity,
        x_p: u.Quantity,
        A: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            N_A: Acceptors atoms concentration. Expected units: atoms/cm³
            x_p: Width of the depletion layer in the P region. Expected units: cm
            A: Area of the depletion region. Expected units: cm²
        Return:
            Q_p: Charge  stored. Expected units: Coulomb (C)

        Calculate the charge stored on the P  side
        |Q-| = q * A * x_p * N_A    (3.24)
    '''
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")
    if not x_p.unit.is_equivalent(u.cm):
        raise u.UnitsError("Input width of the depletion region x_n must be in units equivalent to cm.")
    if not A.unit.is_equivalent(u.cm**2):
        raise u.UnitsError("Input area A must be in units equivalent to cm².")

    return q * A * x_p * N_A


def W(N_A: u.Quantity = N_A, 
      N_D: u.Quantity = N_D) -> u.Quantity:
    '''
        Parameters:
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
                Default value: N_A = 1e16 atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³
                Default value: N_A = 1e17 atoms/cm³
        Return:
            W: Width of the depletion region. Expected units: cm

        Calculate the width of the depletion region
        W = x_n + x_p = sqrt((2 * eps_si/q) * (1/N_A + 1/ N_D) * V_O) (3.26)
    '''
    # Validate input parameters
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")

    return ((2 * eps_si / q) * (1/N_A + 1/N_D) * V_O(N_A, N_D))**0.5




def x_n(W: u.Quantity,
        N_A: u.Quantity,
        N_D: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            W: Width of the depletion region. Expected units: cm
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³
        Return:
            x_n, width of the depletion layer in the N region [m]

        Calculate the width of the depletion layer in the N region:
        x_n  = W * N_A/(N_A + N_D) (3.27)

    '''
    if not W.unit.is_equivalent(u.cm):
        raise u.UnitsError("Input width of the depletion region W must be in units equivalent to cm.")
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")

    return W * N_A/(N_A + N_D)

def x_p(W: u.Quantity,
        N_A: u.Quantity,
        N_D: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            W: Width of the depletion region. Expected units: cm
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³
        Return:
            x_p, width of the depletion layer in the P region [m]

        Calculate the width of the depletion layer in the P region:
        x_p  = W * N_A/(N_A + N_D) (3.28)
    '''
    if not W.unit.is_equivalent(u.cm):
        raise u.UnitsError("Input width of the depletion region W must be in units equivalent to cm.")
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")

    return W * N_A/(N_A + N_D)

def Q_J(N_A: u.Quantity, 
         N_D: u.Quantity, 
         A:   u.Quantity) -> u.Quantity:
    '''
        Parameters:
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³
            A: Area of the depletion region. Expected units: cm²  
        Return 
             Q_J: Charge  stored. Expected units: Coulomb (C)

        Calculate the charge on either side of the depletion region
        when a reverse voltage is applied:
        Q_J = A* sqrt(2 * eps_si * q *((N_A * N_D)/(N_A + N_D))*(V_O))  (3.30)
    '''
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.") 
    if not A.unit.is_equivalent(u.cm**2):
        raise u.UnitsError("Input area A must be in units equivalent to cm².")
  
    return (A * (2 * eps_si * q * ((N_A * N_D)/(N_A + N_D)) * (V_O(N_A, N_D)))**0.5).to(u.C)



def W_R(V_R: u.Quantity,
        N_A: u.Quantity = N_A,
        N_D: u.Quantity = N_D) -> u.Quantity:
    '''
        Parameters:
            V_R: Reverse potential. Expected units: Volts (V)
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³
        Return:
            W_R: Width of the depletion region. Expected units: cm

        Calculate the width of the depletion region, W, when a reverse 
        voltage, V_R, is applied
        W = x_n + x_p = sqrt((2 * eps_si / q) * (1/N_A + 1/ N_D) * (V_O + V_R)) (3.31)
    '''
    if not V_R.unit.is_equivalent(u.Volt):
        raise u.UnitsError("Input reverse potential V_R must be in units equivalent to Volt.")
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_D must be in units equivalent to 1/cm**3.")
  
    return  (((2 * eps_si / q) * (1/N_A + 1/N_D) * (V_O(N_A, N_D) + V_R))**0.5).to(u.cm)

def V_T(T: u.Quantity = T) -> u.Quantity:
    '''
        Parameters:
            T: Temperature. Expected units: Kelvin (K)
        Return:
            V_T: Thermal voltaje. Expected units: Volts (V)

        Calculate the thermal voltage:
        V_T = k * T / q
    '''
    # Validate input parameters
    if not T.unit.is_equivalent(u.Kelvin):
        raise u.UnitsError("Input temperature T must be in units equivalent to Kelvin.")

    return (k_J_K * T / q).to(u.V)

def I(V: u.Quantity,
      I_S: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            V: Voltage. Expected units: Volts (V)
            I_S: Saturation current. Expected units: Amperes (A) 
        Return:
            I: Diode curren. Expected units: (A)
        
        Diode current  current when a forward bias voltage is applied
        I = I_S * (exp(V / V_T) - 1)    (3.40)
  '''

    # Validate input parameters
    if not V.unit.is_equivalent(u.Volt):
        raise u.UnitsError("Input voltage V must be in units equivalent to Volt.")
    if not I_S.unit.is_equivalent(u.ampere):
        raise u.UnitsError("Input voltage V must be in units equivalent to Ampere.")

    return I_S * (m.exp(V / V_T()) - 1)


def QJ_R(N_A: u.Quantity, 
         N_D: u.Quantity, 
         A:   u.Quantity,
         V_R: u.Quantity ) -> u.Quantity: 
    '''
        Parameters:
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³
            A: Area of the depletion region. Expected units: cm²
            V_R: Reverse voltage. Expected units: Volts (V)
        Return:
            Qj: Charge  stored. Expected units: Coulomb (C)

        Calculate the charge on either side of the depletion region
        when a reverse voltage is applied:
        Qj = A* sqrt(2 * eps_si * q *((N_A * N_D)/(N_A + N_D))*(V_O + V_R) )  (3.32)
  '''
    # Validate input parameters
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")
    if not A.unit.is_equivalent(u.cm**2):
        raise u.UnitsError("Input area A must be in units equivalent to cm².")
    if not V_R.unit.is_equivalent(u.Volt):
        raise u.UnitsError("Input reverse voltage V_R must be in units equivalent to Volt.")
    
    sqrt = (2 * eps_si * q * ((N_A * N_D)/(N_A + N_D)) * (V_O(N_A, N_D) + V_R))**0.5
    return (A * sqrt).to(u.C)


def C_JO(A: u.Quantity,
         N_A: u.Quantity,
         N_D: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            A  : Area of the depletion region. Expected units: cm²
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³
        Return:
            C_JO: Capacitance at the junction. Expected units: Farad (F)

        Calculate the capacitance 
        C_JO = A * sqrt((eps_si * q)/ 2) * (N_A * N_D / (N_A + N_D) * (1/V_O) ) (3.48)
    '''
    # Validate input parameters
    if not A.unit.is_equivalent(u.cm**2):
        raise u.UnitsError("Input area A must be in units equivalent to cm².")
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")

    sqrt = (((eps_si * q) / 2) * ((N_A * N_D) / (N_A + N_D)) * (1 / V_O(N_A, N_D)))**0.5
    return (A * sqrt).to(u.Farad)

def I_S(A: u.Quantity, 
        N_A: u.Quantity, 
        N_D: u.Quantity,
        L_n: u.Quantity,
        L_p: u.Quantity,
        D_n:u.Quantity,
        D_p:u.Quantity) -> u.Quantity:
    '''
        Parameters:
            A: Area of the depletion region. Expected units: cm²
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³
            L_n: Diffusion Length of electrons in the P material. Expected units: microns
            L_p: Diffusion Length of holes in the N material microns. Expected units: microns
            D_n: Diffusion constant of electrons in the P material. Expected units: cm²/s
            D_p: Diffusion constant of holes in the N material. Expected units: cm²/s
        Return:
            I_S: Saturation current. Expected units: A

        Calculate the saturation current:
        I_S = A * q * ni²(D_p/(L_p * N_D) + Dn/(Ln * N_A)) (3.41)
  '''
    # Validate input parameters
    if not A.unit.is_equivalent(u.cm**2):
        raise u.UnitsError("Input area A must be in units equivalent to cm².")
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")
    if not L_n.unit.is_equivalent(u.um):
        raise u.UnitsError("Input diffusion length L_n must be in units equivalent to microns.")
    if not L_p.unit.is_equivalent(u.um):
        raise u.UnitsError("Input diffusion length L_p must be in units equivalent to microns.")
    if not D_n.unit.is_equivalent(u.cm**2/u.second):
        raise u.UnitsError("Input diffusion constant D_n must be in units equivalent to cm²/s.")
    if not D_p.unit.is_equivalent(u.cm**2/u.second):
        raise u.UnitsError("Input diffusion constant D_p must be in units equivalent to cm²/s.")

    return (A * q * n_i()**2 * (D_p/(L_p.to(u.cm) * N_D) + D_n/(L_n.to(u.cm) * N_A))).to(u.ampere)



def C_J(V_R: u.Quantity,
        A: u.Quantity,
        N_A: u.Quantity,
        N_D: u.Quantity ) -> u.Quantity:
    '''
        Parameters:
            V_R: Reverse voltage. Expected units: Volts (V)
            A  : Area of the depletion region. Expected units: cm²
            N_A: Acceptor atoms concentration. Expected units: atoms/cm³
            N_D: Donor atoms concentration. Expected units: atoms/cm³  
        Return:
            C_J: Capacitance on the junction. Expected units: Farad (F)
    
        Calculate the capacitance at the junction 
        C_J = C_JO/(1 + V_R/V_O]])**m  (3.49)
    '''
    # Validate input parameters
    if not V_R.unit.is_equivalent(u.Volt):
        raise u.UnitsError("Input voltage V_R must be in units equivalent to volts.")
    if not A.unit.is_equivalent(u.cm**2):
        raise u.UnitsError("Input area A must be in units equivalent to cm².")
    if not N_A.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input acceptor atoms density N_A must be in units equivalent to 1/cm**3.")
    if not N_D.unit.is_equivalent(1/u.cm**3):
        raise u.UnitsError("Input donor atoms density N_D must be in units equivalent to 1/cm**3.")
  
    C_J = C_JO(A, N_A,N_D)/(1 + V_R / V_O(N_A, N_D))**0.5
    return C_J.to(u.Farad)


def tau_p(L_p: u.Quantity,
        D_p: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            L_p: Diffusion Length of holes in the P material. Expected units: microns
            D_p: Diffusion constant of holes in the P material. Expected units: cm²/s
        Return
            tau_p: Lifetime of holes. Expected units: s

        Calculate the minority carrier lifetime of holes
        tau_p = (L_p)² / D_p  (3.51)
  '''
    # Validate input parameters
    if not L_p.unit.is_equivalent(u.um):
        raise u.UnitsError("Input diffusion length L_p must be in units equivalent to microns.")
    if not D_p.unit.is_equivalent(u.cm**2/u.second):
        raise u.UnitsError("Input diffusion constant D_p must be in units equivalent to cm²/s.")
    
    return ((L_p.to(u.cm))**2) / D_p


def C_d(tau_T: u.Quantity,
        V_T: u.Quantity,
        I: u.Quantity) -> u.Quantity:
    '''
        Parameters:
            tau_T: Mean transit time of the junction. Expected units: s
            V_T: Thermal voltage. Expected units: Volts (V)
            I: Forward bias current. Expected units: Ampere (A)
        Return:
            C_d: Incremental diffusion capacitance. Expected units: F

        Calculate the incremental Diffusion capacitance:
        C_d = (tau_T / V_T) * I  (3.57)
    '''
    # Validate input parameters
    if not tau_T.unit.is_equivalent(u.second):
        raise u.UnitsError("Input mean transit time tau_T must be in units equivalent to s.")
    if not V_T.unit.is_equivalent(u.Volt):
        raise u.UnitsError("Input thermal voltage V_T must be in units equivalent to Volts.")
    if not I.unit.is_equivalent(u.ampere):
        raise u.UnitsError("Input forward bias current I must be in units equivalent to Ampere.")

    return ((tau_T / V_T ) * I).to(u.Farad)
