# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 11:39:27 2021

@author: chris
"""
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import fsolve

def choose_comps(comp_names: list[str], file: str, sheet: str, columns: list[str], rows: int) -> pd.DataFrame:
    """reads specified columns from excel sheet and returns rows with data for
    relevant substances as a pandas dataframe
    """

    choose = pd.read_excel(file, sheet_name=sheet, usecols=columns,
                            skiprows=rows, index_col=0)
    choose = choose.loc[comp_names,:]
    return choose

def compile_data(file: str, comp_names: list[str]) -> tuple[pd.DataFrame]:
    """compiles data from excel file for given components
    """

    #general data
    sheet: str = "allgemein"
    columns = ["Substance","Formula","Tc","Pc","Dh0f"]
    df_gen = choose_comps(comp_names, file, sheet, columns, 3)

    # vapour pressure data
    sheet = "Dampfdruck"
    columns = ["Substance","Formula","A","B","C","D"]
    df_vap = choose_comps(comp_names, file, sheet, columns, 8)

    # enthalpy data
    sheet = "dhv"
    columns = ["Substance","Formula","A","B","C","D","E"]
    df_enth = choose_comps(comp_names, file, sheet, columns, 8)

    # cp data
    sheet = "cp ideal"
    columns = ["Substance","Formula","A","B","C","D","E","F","G"]
    df_cp = choose_comps(comp_names, file, sheet, columns, 13)
    df_cp.fillna(0, inplace = True)

    return df_gen, df_vap, df_enth, df_cp

def estimate_flow(n: int, F: float, D: float, nu: float, f: float, n_feed: int, G_j: float, L_j: float) -> None:
    """calculates molar flow rates G_j and L_j by assuming constant molar 
    overflow as a first estimate, used at the start of the calculation
    """
    
    # general mass balance
    G_j[0] = D          # putting head flow in G vector
    L_j[-1] = F - D    # general balance: putting bottom flow in L vector
    
    # balance around head
    L_j[0] = D * nu     # putting head reflux in L vector
    G_j[1] = D + L_j[0] 
    
    # loop over all plates
    for j in range(1, n-1):
        
        # rectifying section
        if 0 <= j <= n_feed-2:
            if n_feed >= 2:     # else no rectifying section in column
                L_j[j] = L_j[0]
                G_j[j+1] = G_j[1]
                
        # vapour feed
        if j == n_feed-1:
            G_j[n_feed] = G_j[1]
            L_j[n_feed-1] = L_j[0]
        
        # stripping section
        if j >= n_feed:
            L_j[j] = L_j[0] + (1 - f) * F
            G_j[j+1] = L_j[j] - L_j[-1]

def calc_vap(df_gen: pd.DataFrame, df_vap: pd.DataFrame, T: np.ndarray, p_vap: pd.DataFrame) -> np.ndarray:
    """calculates vapour pressures p_vap for all components on all plates
    using 2.5,5 - Wagner equation and estimated temperature
    """
    
    A = df_vap.loc[:,"A"].to_numpy()
    B = df_vap.loc[:,"B"].to_numpy()
    C = df_vap.loc[:,"C"].to_numpy()
    D = df_vap.loc[:,"D"].to_numpy()
    Tc = df_gen.loc[:,"Tc"].to_numpy()
    Pc = df_gen.loc[:,"Pc"].to_numpy()
    
    p = np.copy(p_vap)
    T_r = np.ones(p.shape)
    T_r = (T_r.T * T).T / Tc
    t = 1-T_r
    
    p = (A * t + B * t**1.5 + C * t**2.5 + D * t**5) / T_r
    p = np.exp(p.astype(float)) * Pc.astype(float)
    
    return p

def solve_system(A_ji: np.ndarray, v_ji: np.ndarray, L_i: np.ndarray, n: int, m: int) -> None:
    """solves system of linear equations by inverting the matrix of
    coefficients, saves solution in v_ji
    """
    
    for i in range(1, m+1):             # loop over components
        
        # initialising vectors to use as diagonals in matrix
        lower_diag = np.delete(A_ji[:, i-1], n-1)
        upper_diag = np.ones(n-1)
        main_diag = -(1 + A_ji[:, i-1])
        
        # initialising matrix for component i
        M_i = (np.zeros((n,n)) + np.diag(main_diag, 0) + np.diag(upper_diag, 1)
               + np.diag(lower_diag, -1))
              
        #inverting matrix and solving equations
        v_ji[:, i-1] = np.dot(np.linalg.inv(M_i), (-1)*L_i[:, i-1])
        
def print_current(G_j: np.ndarray, L_j: np.ndarray, v_ji: np.ndarray, l_ji: np.ndarray, y_ji: np.ndarray, x_ji: np.ndarray, T: np.ndarray, theta: float) -> None:
    """prints current flow rates, temperatures and theta
    """
    print("Temperature on plates [K]: \t", T)
    print("Total vapour molar flow [mol/h]: \t", G_j)
    print("Component vapour molar flow [mol/h]: \n", v_ji)
    print("Vapour molar flow composition: \n", y_ji)
    print("Total liquid molar flow [mol/h]: \t", L_j)
    print("Component liquid molar flow [mol/h]: \n", l_ji)
    print("Liquid molar flow composition: \n", x_ji)
    print("Current theta weight factor: \t", theta)

def theta_fsolve(theta: float, d_icalc: np.ndarray, b_icalc: np.ndarray, F: float, z_F: np.ndarray, D: float) -> float:
    """generates function for fsolve to calculate theta, function value is 0
    for the correct theta
    """
    
    f = np.sum(F * z_F / (1 + theta * b_icalc / d_icalc)) - D
    
    return f

def T_fsolve(T, x_jicor, df_gen, df_vap, p_vap, p_tot) -> float:
    """generates function for fsolve to calculate T from VLE equilibrium;
    function value is 0 for correct T
    """
    
    if T.size == 1:     # for scalar temperature
        g = np.sum(calc_vap(df_gen, df_vap, T, p_vap) / p_tot * x_jicor) - 1

    else:               # for temperature vector
        g = np.sum(calc_vap(df_gen, df_vap, T, p_vap) 
                   / p_tot * x_jicor, axis=1) - 1
    
    return g

def calc_dhv(T, Tc, R, A, B, C, D, E) -> float:
    """calculates vapourisation enthalpy of one component at given temperature
    using data from Excel sheet, returns scalar
    """
    
    tau = 1 - T/Tc
    dhv = R*Tc * (A*tau**(1/3) + B*tau**(2/3) + C*tau + D*tau**2 + E*tau**6)
    return dhv

def calc_cp(T, R, A, B, C, D, E, F, G) -> float:
    """calculates cp of one component at given temperature, returns scalar
    """
    
    y = T / (A + T)
    cp = R * (B + (C-B)*y**2 * (1 + (y-1) * (D + E*y + F*y**2 + G*y**3)))
    return cp

def pure_enthalpy(df_gen, df_cp, df_enth, T: np.ndarray, h_gas: np.ndarray, h_liq: np.ndarray, R: float) -> None:
    """calculates liquid and vapour enthalpy of all components using data from
    excel sheet, fills matrices h_gas and h_liq
    """
    
    # converting dataframes to arrays
    dh0f = df_gen.loc[:,"Dh0f"].to_numpy()
    Tc = df_gen.loc[:,"Tc"].to_numpy()
    coeff_cp = df_cp[["A", "B", "C", "D", "E", "F", "G"]].to_numpy()
    coeff_dhv = df_enth[["A", "B", "C", "D", "E"]].to_numpy()
    
    dhv = np.zeros(h_gas.shape)
    T_ST = 298.15 # K

    if isinstance(T, (int, float)):
        T = [T]
    
    # loop over h_gas
    for j in range(np.size(h_gas,0)):       # loop over plates
        for i in range(np.size(h_gas,1)):   # loop over components
            
            # integrating cp over T
            quad_tmp = quad(calc_cp, T_ST, T[j],
                            args=(R, coeff_cp[i,0], coeff_cp[i,1], coeff_cp[i,2],
                                  coeff_cp[i,3], coeff_cp[i,4], coeff_cp[i,5],
                                  coeff_cp[i,6]))
            # vapour enthalpy
            h_gas[j,i] = dh0f[i] + quad_tmp[0]
            # delta enthalpy
            dhv[j,i] = calc_dhv(T[j], Tc[i], R, coeff_dhv[i,0],
                                coeff_dhv[i,1], coeff_dhv[i,2],
                                coeff_dhv[i,3], coeff_dhv[i,4])
            # liquid enthalpy
            h_liq[j,i] = h_gas[j,i] - dhv[j,i]

def rachford_rice(f: float, z: np.ndarray, K: np.ndarray) -> float:
    """generates function for fsolve to calculate vapour amount in feed using
    feed temperature and K constants, returns 0 for correct vapour ratio
    """
    c = 1/(1-K)
    d = (c[0]-c)/(c[-1]-c[0])
    def F(a: float) -> float:
        return np.sum(z*a/(d+a*(1+d)))
    def G(a: float) -> float:
        return (1+a)/a*F(a)
    def H(a: float) -> float:
        return -a*F(a)
    aL = z[0]/(1-z[0])
    aR = (1-z[-1])/z[-1]
    fal = F(aL)
    far = F(aR)
    # fal, far
    a0 = aL - fal / (far-fal) * (aR-aL)
    a = f
    if G(a0)>0:
        a, = fsolve(G, a0)
    else:
        a, = fsolve(H, a0)
    V = (c[0] + a*c[-1]) / (1+a)
    return V
    
def define_feed(f0: float, z_F: np.ndarray, K_feed: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """calculates amount and composition of feed phases
    """
    f = rachford_rice(f0, z_F, K_feed[0])
    if f<=0: 
        f = 0
        x_F = z_F
        y_F = np.zeros(z_F.shape)
    elif f>=1:
        f = 1
        y_F = z_F
        x_F = np.zeros(z_F.shape)
    else:
        x_F = z_F / (1 + f*(K_feed -1))
        y_F = x_F * K_feed
    return f, x_F, y_F