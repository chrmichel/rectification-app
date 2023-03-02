import numpy as np
import matplotlib.pyplot as plt
from configparser import ConfigParser
from rect_modules import *

FILE = "Eingabedatei.xlsx"
R = 8.3141 # J/molK

class Rectification:

    def __init__(self, input_dict):
        self.comps: list[str] = input_dict["comps"]
        self.m = len(self.comps)
        self.z_F = np.array(input_dict["z_f"])
        self.n = int(input_dict["n"]) + 2
        self.T_F = float(input_dict["t_f"])
        self.n_F = int(input_dict["n_f"])
        self.total_boil = bool(input_dict["boil"])
        self.total_cond = bool(input_dict["cond"])
        self.F = float(input_dict["f"])
        self.D = float(input_dict["d"])
        self.nu = float(input_dict["nu"])
        self.f = float(input_dict["eta"])
        self.p_tot = float(input_dict["p_tot"])

        self.data = compile_data(FILE, self.comps)
        self.setup()
        self.results: dict = {}
    

    def setup(self):
        self.G_j = np.zeros(self.n)
        self.L_j = np.zeros(self.n)
        self.T_j = self.T_F * np.ones(self.n)
        # vapour pressure of all comps on all plates [bar]
        self.p_vap = np.zeros((self.n, self.m))
        # Equilibrium constants [-]
        self.K_ji = np.zeros((self.n, self.m))  
        # Absorption coefficients [-]
        self.A_ji = np.zeros((self.n, self.m))  
        # vapuor comp flow rates [mol/h]
        self.v_ji = np.zeros((self.n, self.m))  
        # liqiud comp flow rates [mol/h]
        self.l_ji = np.zeros((self.n, self.m)) 
        # vapour composition [-] 
        self.y_ji = np.zeros((self.n, self.m))  
        # liquid composition [-]
        self.x_ji = np.zeros((self.n, self.m))  

        self.theta = 0.0                # weight factor
        # corrected vapuor comp flow rates [mol/h]
        self.v_jicor = np.zeros((self.n, self.m))  
        # corrected liqiud comp flow rates [mol/h]
        self.l_jicor = np.zeros((self.n, self.m))  
        # corrected vapour composition [-]
        self.y_jicor = np.zeros((self.n, self.m))
        # corrected liquid composition [-]  
        self.x_jicor = np.zeros((self.n, self.m))  

        # enthalpy of vapour components [J/mol]
        self.h_gas = np.zeros((self.n, self.m))    
        # enthalpy of liquid components [J/mol] 
        self.h_liq = np.zeros((self.n, self.m))
        # vapourization enthalpy [J/mol]
        self.delta_hv = np.zeros((self.n, self.m))     

        # vapour pressure of feed components [bar]
        self.p_feed = np.zeros((1, self.m))
        # feed liquid composition [-]
        self.x_F = np.zeros_like(self.z_F)
        # gas liquid composition [-]
        self.y_F = np.zeros_like(self.z_F)
        # enthalpy of vapour feed components [bar]
        self.h_feed_gas = np.zeros((1, self.m))
        # enthalpy of liquid feed components [bar]
        self.h_feed_liq = np.zeros((1, self.m))
        # feed enthalpy [J/mol]
        self.hF = 0.
        self.hFG = 0.
        self.hFL = 0.
        # vapour composition [-]
        self.y_old = np.zeros((self.n, self.m))     
        # liquid composition [-]
        self.x_old = np.zeros((self.n, self.m))     


    def calc_feed(self):
        self.p_feed = calc_vap(self.data[0], self.data[1],
                                self.T_F, self.p_feed)
        self.K_F = self.p_feed / self.p_tot   # feed equilibrium constants [-]
        self.f, self.x_F, self.y_F = define_feed(self.f, self.z_F, self.K_F)
        pure_enthalpy(self.data[0], self.data[3], self.data[2], self.T_F,
                        self.h_feed_gas, self.h_feed_liq, R)
        self.hFG = np.sum(self.y_F * self.h_feed_gas)
        self.hFL = np.sum(self.x_F * self.h_feed_liq)
        self.hF = self.f * self.hFG + (1 - self.f) * self.hFL


    def run(self):
        # estimate molar flow on all plates
        estimate_flow(self.n, self.F, self.D, self.nu, self.f, self.n_F,
                        self.G_j, self.L_j)
        # print(self.G_j, self.L_j)
        # set up iteration loop
        dx = 1; dy = 1; counter = 0; eps = 1e-4
        while ((dx>eps or dy>eps) and counter<20):
            
            # calculate feed flash and enthalpy
            self.calc_feed()
            # calculate vapor pressure
            self.p_vap = calc_vap(self.data[0], self.data[1],
                                    self.T_j, self.p_vap)

            # prepare system of linear equations
            L_i = np.zeros((self.n, self.m))
            L_i[self.n_F, :] = self.F * (1 - self.f) * self.x_F
            L_i[self.n_F - 1, :] = self.F * self.f * self.y_F
            # print(L_i)

            self.K_ji = self.p_vap / self.p_tot
            self.A_ji = (self.L_j / self.G_j / self.K_ji.T).T
            if self.total_cond:
                self.A_ji[0] = self.L_j[0] / self.G_j[0]
            if self.total_boil:
                self.A_ji[-1] = self.L_j[-1] / self.G_j[-1]

            solve_system(self.A_ji, self.v_ji, L_i, self.n, self.m)
            self.l_ji = self.A_ji * self.v_ji

            self.y_ji = (self.v_ji.T / self.G_j).T
            self.x_ji = (self.l_ji.T / self.L_j).T

            ######################################
            # theta covergence                   #
            ######################################

            d_iber = self.v_ji[0]
            b_iber = self.l_ji[-1]

            self.theta = fsolve(theta_fsolve, self.theta,
                                  (d_iber, b_iber, self.F, self.z_F, self.D))
            d_ikor = self.F * self.z_F / (1 + self.theta * (b_iber/d_iber))
            self.l_jicor = d_ikor/d_iber * self.l_ji
            self.v_jicor = d_ikor/d_iber * self.v_ji

            for k in range(self.n):
                self.y_jicor[k] = self.v_jicor[k]/np.sum(self.v_jicor[k])
                self.x_jicor[k] = self.l_jicor[k]/np.sum(self.l_jicor[k])

            ######################################
            # boiling point method               #
            ######################################
            
            self.T_j = fsolve(T_fsolve, self.T_j, (self.x_jicor, self.data[0],
                                self.data[1], self.p_vap, self.p_tot))
            pure_enthalpy(self.data[0], self.data[3], self.data[2], self.T_j,
                            self.h_gas, self.h_liq, R)
            
            h_head = np.sum(self.x_jicor[0] * self.h_liq[0])
            h_bottom = np.sum(self.x_jicor[-1] * self.h_liq[-1])
            hfGxf = np.sum(self.x_jicor[self.n_F-1] * self.h_gas[self.n_F-1])

            hjGxj = np.sum(self.h_gas * self.x_jicor, axis=1)
            hjGxh = np.sum(self.h_gas * self.x_jicor[0], axis=1)
            hjLxb = np.sum(self.h_liq * self.x_jicor[-1], axis=1)
            hjLxj = np.sum(self.h_liq * self.x_jicor, axis=1)
            hjGyj = np.zeros(self.n)
            hjLyj = np.zeros(self.n)

            for j in range(self.n - 1):
                hjGyj[j] = np.sum(self.h_gas[j] * self.y_jicor[j+1])
                hjLyj[j] = np.sum(self.h_liq[j] * self.y_jicor[j+1])

            ######################################
            # energy balances                    #
            ######################################
            
            Q_H = self.L_j[0] * (hjGxj[0] - hjLxj[0]) \
                 + self.G_j[0] * (hjGxh[0] - h_head)
            Q_B = self.G_j[0] * h_head + self.L_j[-1] * h_bottom \
                 + Q_H - self.F * self.hF

            for j in range(self.n - 1):

                if 1 <= j <= self.n_F-2:
                    if self.n_F >= 2:
                        self.L_j[j] = (self.G_j[0] * (h_head - hjGxh[j])
                            + Q_H) / (hjGxj[j] - hjLxj[j])
                        self.G_j[j+1] = self.L_j[j] + self.G_j[0]
                
                elif j == self.n_F - 1:
                    self.L_j[j] = (self.G_j[0] * (h_head - hjGxh[j]) + Q_H
                            + self.F * self.f * (hfGxf - self.hFG)) \
                        / (hjGxj[j] - hjLxj[j])
                    self.G_j[j+1] = self.L_j[j] + self.G_j[0]
                
                elif j >= self.n_F:
                    self.G_j[j+1] = (self.L_j[-1] * (hjLxb[j] - h_bottom)\
                         + Q_B) / (hjGyj[j] - hjLyj[j])
                    self.L_j[j] = self.G_j[j+1] + self.L_j[-1]

            for k in range(self.n):
                self.v_ji[k] = self.y_jicor[k] * self.G_j[k]
                self.l_ji[k] = self.x_jicor[k] * self.L_j[k]

            # check convergence
            dx = np.amax(np.abs(self.x_old - self.x_jicor))
            dy = np.amax(np.abs(self.y_old - self.y_jicor))

            self.x_old = np.copy(self.x_jicor)
            self.y_old = np.copy(self.y_jicor)

            counter += 1
        
        self.results["T"] = self.T_j
        self.results["x"] = self.x_jicor
        self.results["y"] = self.y_jicor


if __name__ == '__main__':
    cfg = ConfigParser(converters={'list': lambda x:
                            [i.strip() for i in x.split()]})
    cfg.read("config.ini")
    input_dict = dict(cfg.items("DATA"))
    input_dict["comps"] = [i.strip() for i in input_dict["comps"].split()]
    input_dict["z_f"] = [float(i.strip()) for i in input_dict["z_f"].split()]
    r = Rectification(input_dict)
    r.run()