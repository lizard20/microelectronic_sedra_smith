'''
    Solved exercises from chapter three

    Author: Aldo Nunez
'''

import semiconductor as sc
from astropy import units as u 

def exercise_3_1():
  T_1 = 50 * u.Kelvin
  T_2 = 350 * u.Kelvin

  n_i_1 = sc.n_i(T_1)
  n_i_2 = sc.n_i(T_2)

  print(f"Exercise 3.1:")
  print(f"n_i({T_1:.0f}) = {n_i_1:0.2e}")
  print(f"n_i({T_2:.0f}) = {n_i_2:0.2e}")

def exercise_3_2():
    T = 350 * u.Kelvin   
    N_D = 1e17 * 1/u.cm**3  

    n_n = sc.N_D
    p_n = sc.p_n(T, N_D)

    print(f"Exercise 3.2:")
    print(f"n_n({T:.0f}) = {n_n:.0e}")
    print(f"p_n({T:.0f}) = {p_n:1.2e}")

def exercise_3_3():
    T = 300 * u.Kelvin
    factor = 1e6 * u.dimensionless_unscaled

    N_A = sc.n_i(T) * factor

    print(f"Exercise 3.3:")
    print(f"N_A({T:.0f}) = {N_A:0.1e}")

def exercise_3_4():
    x = 2 * u.um
    V = 1 * u.Volt
    N_D = 1e16 * 1/u.cm**3
    mu_n = 1_350 * u.cm**2/(u.Volt * u.second)
    A = 0.25 * u.um**2

    # Part a), calculate the Electric field, E
    E = -V/x.to(u.cm)
    vn_drift  = sc.vn_drift(E, mu_n)
    
    # Part b), calculate the time, t
    t = x.to(u.cm) / vn_drift

    # Part c), calculate the drift current density, J
    p_n = sc.p_n(N_D = N_D)
    sigma = sc.sigma(n = N_D, p = p_n)
    J = sc.J_sigma(E, sigma)

    # Part d), calculate the current, I
    I = J * A.to(u.cm**2)

    print(f"Exercise 3.4:")
    print("Intermediate parameters:")
    print(f"E = {E:0.2e}")
    print(f"p_n = {p_n:0.2e}")
    print(f"sigma = {sigma.to(u.siemens / u.cm):0.2e}")
    print(f"----------------------------\n")
    print(f"vn_drift = {vn_drift: 0.2e}")
    print(f"t = {t.to(u.ps):0.1f}")
    print(f"J = {J: 0.2e}")
    print(f"I = {I.to(u.uA): .0f}")

def exercise_3_5():
    n_0 = 1e17 * 1/u.cm**3
    W = 1 * u.um
    D_n = 35 * u.cm**2 /u.second    
    I_n = 1e-3 * u.ampere

    # gradient
    dn_dx = -n_0 / W.to(u.cm)
    J_n = sc.J_n(D_n, dn_dx).to(u.A / u.um**2)
    A = I_n / abs(J_n)

    print(f"Exercise 3.5:")
    print(f"J_n = {J_n.to(u.uA/u.um**2):0.2f}")
    print(f"A = {A:0.2f}")

def exercise_3_6():
    D_n = sc.V_T() * sc.mu_n
    D_p = sc.V_T() * sc.mu_p

    print(f"Exercise 3.6:")
    print(f"D_n = {D_n.to(u.cm**2/(u.second)):0.1f}")
    print(f"D_p = {D_p:0.2f}")

def exercise_3_11():
    A = 1e-4 * u.cm**2
    N_A = 1e18 * 1/u.cm**3
    N_D = 1e16 * 1/u.cm**3
    L_p = 5 * u.um
    L_n = 10 * u.um
    D_p = 10 * u.cm**2 / u.second
    D_n = 18 * u.cm**2 / u.second
    V = 0.605 * u.Volt

    I_S = sc.I_S(A, N_A, N_D/2, L_n, L_p, D_n, D_p)
    I = sc.I(V, I_S)

    print(f"Exercise 3.11:")
    print(f"I_S = {I_S: 0.2e}")
    print(f"I = {I.to(u.mA): 0.2f}")

def exercise_3_12(): 
    N_A = 1e18 * 1/u.cm**3
    N_D = 1e16 * 1/u.cm**3
    V_F = -0.605 * u.Volt    

    W = sc.W_R(V_F, N_A, N_D).to(u.um)

    print(f"Exercise 3.12:")
    print(f"W = {W: 0.3f}")

def exercise_3_13():
    N_A = 1e18 * 1/u.cm**3
    N_D = 1e16 * 1/u.cm**3
    V_R = 2 * u.Volt
    A = 1e-4 * u.cm**2   
    L_p = 5  * u.um    
    L_n = 10 * u.um   
    D_p = 10 * u.cm**2 / u.second
    D_n = 18 * u.cm**2 / u.second

    W = sc.W_R(V_R, N_A, N_D).to(u.um)
    Q_J = sc.QJ_R(N_A, N_D, A, V_R)
    I_S = sc.I_S(A, N_A, N_D, L_n, L_p, D_n, D_p)

    print(f"Exercise 3.13:")
    print(f"W = {W.to(u.um): 0.3f}")
    print(f"Q_J = {Q_J.to(u.pC): 0.2f}")
    print(f"I_S = {I_S: 0.2e}")

def exercise_3_14():
    V_R = 2 * u.Volt
    N_A = 1e18 * 1 / u.cm**3
    N_D = 1e16 * 1 / u.cm**3
    A   = 1e-4 * u.cm**2

    C_JO = sc.C_JO(A, N_A, N_D)
    C_J = sc.C_J(V_R, A, N_A, N_D)

    print(f"Exercise 3.14:")
    print(f"C_JO = {C_JO.to(u.pF): 0.1f}")
    print(f"C_J =  {C_J.to(u.pF): 0.2f}")

def exercise_3_16():
    D_p = 10 * u.cm**2/u.s
    L_p = 5 * u.um
    I_p = 0.0991e-3 * u.ampere

    tau_p = sc.tau_p(L_p, D_p)  
    C_d = sc.C_d(tau_p, sc.V_T(), I_p)

    print(f"Exercise 3.16:")
    print(f"tau_p = {tau_p.to(u.ns):2.1f}")
    print(f"C_d = {C_d.to(u.pF):2.1f}")

def main():
    exercise_3_1()
    input("\nPress Enter")
    exercise_3_2()
    input("\nPress Enter")
    exercise_3_3()
    input("\nPress Enter")
    exercise_3_4()
    input("\nPress Enter")
    exercise_3_5()
    input("\nPress Enter")
    exercise_3_6()
    input("\nPress Enter")
    exercise_3_11()
    input("\nPress Enter")
    exercise_3_12()
    input("\nPress Enter")
    exercise_3_13()
    input("\nPress Enter")
    exercise_3_14()
    input("\nPress Enter")
    exercise_3_16()

if __name__ == "__main__":
    main()


