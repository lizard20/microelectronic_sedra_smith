'''
    Solved problems from chapter three

    Aldo Nunez
'''

import semiconductor as sc
from astropy import units as u 
import math as m


def problem_3_1():
    t__55 = -55 * u.deg_C
    t_0   = 0 * u.deg_C
    t_20  = 20 * u.deg_C
    t_75  = 75 * u.deg_C

    ni__55 = sc.n_i(t__55.to(u.Kelvin, equivalencies=u.temperature()))
    ni_0   = sc.n_i(t_0.to(u.Kelvin, equivalencies=u.temperature()))
    ni_20  = sc.n_i(t_20.to(u.Kelvin, equivalencies=u.temperature()))
    ni_75  = sc.n_i(t_75.to(u.Kelvin, equivalencies=u.temperature()))

    fraction__55 = sc.atoms_si/ni__55
    fraction_0 = sc.atoms_si/ni_0
    fraction_20 = sc.atoms_si/ni_20
    fraction_75 = sc.atoms_si/ni_75

    print(f"Problem 3.1:")
    print(f"ni(-55°C) = {ni__55:1.2e}")
    print(f"One out of every {fraction__55:1.2e} atoms is ionized at -55 °C")
    print(f"ni(0°C) = {ni_0:1.2e}")
    print(f"One out of every {fraction_0:1.2e} atoms is ionized at 0 °C")
    print(f"ni(20°C) = {ni_20:1.2e}")
    print(f"One out of every {fraction_20:1.2e} atoms is ionized at 20 °C")
    print(f"ni(75°C) = {ni_75:1.2e}")
    print(f"One out of every {fraction_75:1.2e} atoms is ionized at 75 °C")

def problem_3_3():
    N_A = 5e18 * 1/u.cm**3
    T = 300 * u.Kelvin

    n_p = sc.n_p(T, N_A)

    print(f"Problem 3.3:")
    print(f"p_p({T:.0f}) = {N_A:1.0e}")
    print(f"n_p({T:.0f}) = {n_p:1.2f}")

def problem_3_5():
    N_D = 1e17 * 1/u.cm**3 
    T_1 = 27 * u.deg_C  
    T_2 = 125 * u.deg_C 

    n_n_1 = N_D
    n_n_2 = N_D
    p_n_1 = sc.p_n(T_1.to(u.Kelvin, equivalencies=u.temperature()), N_D)
    p_n_2 = sc.p_n(T_2.to(u.Kelvin, equivalencies=u.temperature()), N_D)

    print(f"Problem 3.5:")
    print(f"n_n({T_1.value :.0f}°C) = {n_n_1}")
    print(f"n_n({T_2.value:.0f}°C) = {n_n_2}")
    print(f"p_n({T_1.value:.0f}°C) = {p_n_1:1.2e}")
    print(f"p_n({T_2.value:.0f}°C) = {p_n_2:1.2e}")

def problem_3_7():
    V = 3 * u.Volt
    x = 10 * u.um
    mu_n = 1350 * u.cm**2/(u.Volt * u.second)
    mu_p = 480 * u.cm**2/(u.Volt * u.second)

    E = -V/x
    v_n_drift = sc.vn_drift(E, mu_n)
    v_p_drift = sc.vn_drift(E, mu_p)

    print(f"Problem 3.7:")
    print(f"vn_drift =  {v_n_drift.to(u.cm/u.second):1.2e}")
    print(f"vp_drift =  {v_p_drift.to(u.cm/u.second):1.2e}")

def problem_3_12():
    N_A = 5e16 * 1/u.cm**3
    N_D = N_A
    A = 20 * u.um**2    

    V_O = sc.V_O(N_A, N_D)

    W = sc.W(N_A, N_D)
    x_n = sc.x_n(W.to(u.cm), N_A, N_D)
    x_p = sc.x_p(W.to(u.cm), N_A, N_D)
    Q_J = sc.Q_J(N_A, N_D, A.to(u.cm**2))

    print(f"Problem 3.12:")
    print(f"V_O =  {V_O.to(u.mV):3.2f}")
    print(f"W =  {W.to(u.um): 1.2f}")
    print(f"x_n =  {x_n.to(u.um):1.1f}")
    print(f"x_p =  {x_p.to(u.um):1.1f}")
    print(f"Q_J =  {Q_J:1.2e}")

def problem_3_14():
    x = 0.1 * u.um
    A = 10 * 10 * u.um**2
    N = 1e18 * 1/u.cm**3

    # we can also use Q_n
    Q = sc.Q_p(N, x.to(u.cm), A.to(u.cm**2))

    print(f"Problem 3.14:")
    print(f"Q = {Q.to(u.pC) :0.2f}")


def problem_3_16():
    factor = 10 * u.dimensionless_unscaled
    
    V_O = sc.V_T() * m.log(factor)

    print(f"Problem 3.16:")
    print(f"V_O = {V_O.to(u.mV):0.2f}")


def problem_3_20():
    V = 750e-3 * u.Volt
    N_A = 1e17 * 1/u.cm**3
    N_D = 1e16 * 1/u.cm**3
    A = 100 * u.um**2
    L_p = 5 * u.um
    L_n  = 10 * u.um
    D_p = 10 * u.cm**2/u.second
    D_n = 18 * u.cm**2/u.second

    I_S = sc.I_S(A, N_A, N_D, L_n, L_p, D_n, D_p)
    I = sc.I(V, I_S)

    print(f"Problem 3.20:")
    print(f"I_S = {I_S:1.2e}")
    print(f"I = {I.to(u.mA):1.2f}")

def problem_3_24():
    N_A = 1e17 * 1/u.cm**3
    N_D = 1e16 * 1/u.cm**3
    A  = 100 * u.um**2
    V_R = 3 * u.Volt

    C_JO = sc.C_JO(A.to(u.cm**2), N_A, N_D)
    C_J = sc.C_J(V_R, A.to(u.cm**2), N_A, N_D)
    
    print(f"Problem 3.24:")
    print(f"C_JO = {C_JO.to(u.fF):1.2f}")
    print(f"C_J= {C_J.to(u.fF):1.2f}")

def problem_3_27():
    I_1 = 1e-3 * u.ampere
    C_d_1 = 5e-12 * u.Farad
    I_2 = 0.1e-3 * u.ampere

    tau_T_1 = sc.V_T() * C_d_1 / I_1 
    C_d_2= sc.C_d(tau_T_1, sc.V_T(), I_2)
    tau_T_2 = sc.V_T() * C_d_2 / I_2 

    print(f"Problem 3.27:")
    print(f"C_d = {C_d_2.to(u.pF): 1.1f}")
    print(f"tau_T = {tau_T_2.to(u.ps): 1.2f}")

def main():
    problem_3_1()
    input("\nPress Enter")
    problem_3_3()
    input("\nPress Enter")
    problem_3_5()
    input("\nPress Enter")
    problem_3_7()
    input("\nPress Enter")
    problem_3_12()
    input("\nPress Enter")
    problem_3_14()
    input("\nPress Enter")
    problem_3_16()
    input("\nPress Enter")
    problem_3_20()
    input("\nPress Enter")
    problem_3_24()
    input("\nPress Enter")
    problem_3_27()

if __name__ == "__main__":
    main()


