# Electron tunneling by Hishinuma's formula
# Gap: fixed value
# Transmission: Non-equilibrium Green's Function (NEGF)

# gold work function:
# P. A. Anderson, Phys. Rev. 115 553 (1959)
# gold fermi energy:
# N. W. Ashcroft and N. D. Mermin, Solid State Physics, 1976

# Coded by Takuro TOKUNAGA
# Last modified: November 10 2020

import numpy as np
import time
import pandas as pd
import sys

start = time.time()

# Unit conversion:
ucev_to_j = 1.602176620898*np.power(10.,-19) # electron volt to joule
ucnano = 1.0*np.power(10.,-9)                # nm to m
ucangs = 1.0*np.power(10.,-10)               # angstrom to m
uccm = 1.0*np.power(10.,-2)                  # cm to m

# physical constants
kb = 1.38064852*np.power(10.,-23)           # Boltzmann constant
ch = 6.62607004*np.power(10.,-34)           # Planck constant
rh = ch/(2*np.pi)                           # reduced Planck constant
c = 299792458                               # light speed, [m]
eps_0 = 8.854187817*np.power(10.,-12)       # [F/m]
sq = 1.602176620898*np.power(10.,-19)       # elementary cahrge [C], small q
me = 9.10938356*np.power(10.,-31)           # electron mass, [kg]
criteria = 708                              # python's expnential arguments limit

# Parameters for bias voltage
volt = -0.6     # [V], load (vias) voltage
#volt = 0.0     # [V], load (vias) voltage
se = 1         # [eV/V], single electron
volt = se*volt # [eV]

# temperature
#temperature_L = 1575 # [K], emitter
#temperature_R = 1000 # [K], collector
#temperature_L = 305 # [K], emitter
#temperature_R = 300 # [K], collector
temperature_L = 295 # [K], emitter
temperature_R = 195 # [K], collector
#temperature_L = 280 # [K], emitter
#temperature_R = 120 # [K], collector
delta_temperature = temperature_L-temperature_R # [K]

# fermi energy
#Ef_L = 0        # [eV], emitter (Devon's Ph.D. thesis)
#Ef_R = volt     # [eV], collector:
Ef_gold = 5.51 # [eV], relative to the emitter fermi energy (Devon's Ph.D. thesis, P.36)
Ef_L = Ef_gold   # [eV]
Ef_R = Ef_gold   # [eV]

# area
#radius = 450*ucnano # [m], 170~250 [nm]
radius = 30*ucnano # [m], 170~250 [nm]
Atip = np.pi*np.power(radius,2.0) # [m2]
#print(str(Atip))

## Parameters (end) ##

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# define Ez left integrand
def NL(arg_Ez): # [eV] emitter
    #term1 = me*kb*temperature_L/(2*np.power(np.pi,2.0)*np.power(rh,3.0))
    term1 = 4*np.pi*me*kb*temperature_L/np.power(ch,3.0)

    #term2 = -(arg_Ez-volt-Ef_L)*ucev_to_j/(kb*temperature_L) # [-] eV to J
    term2 = -(arg_Ez-Ef_L)*ucev_to_j/(kb*temperature_L) # [-] eV to J

    # overflow avoidance
    if term2 > criteria:
        integrand = 0
    else:
        integrand = term1*np.log(1+np.exp(term2)) # P.35, Devon's Ph.D. thesis

    return integrand

# define Ez right integrand
def NR(arg_Ez): # [eV] collector
    #term1 = me*kb*temperature_R/(2*np.power(np.pi,2.0)*np.power(rh,3.0))
    term1 = 4*np.pi*me*kb*temperature_R/np.power(ch,3.0)
    term2 = -(arg_Ez-Ef_R)*ucev_to_j/(kb*temperature_R) # [-] eV to J

    # overflow avoidance
    if term2 > criteria:
        integrand = 0
    else:
        integrand = term1*np.log(1+np.exp(term2)) # P.35, Devon's Ph.D. thesis

    return integrand

# define total integrand, electron tunneling flux
def total(arg_Ez, arg_trans): # [eV]
    term1 = (arg_Ez*ucev_to_j+kb*temperature_L)*NL(arg_Ez)
    term2 = (arg_Ez*ucev_to_j+kb*temperature_R)*NR(arg_Ez-volt)
    term3 = term1 - term2
    integrand = arg_trans*term3 # Eq (3.10)

    return integrand


def fluxNEGF(arg_filenum):

    # file open
    ET = pd.read_csv("../green_d/transmission/trans_"+str(arg_filenum)+".txt", sep=' ', header=None)
    ET.columns = ["E", "T"] # Energy & Transmission
    row, col = ET.shape # row & column of matorix

    # table for integral calculation
    Ez_table = np.zeros(row, dtype='float64')
    trans_table = np.zeros(row, dtype='float64')

    for i in range(0, row):
        Ez_table[i] = ET.iat[i,0] # x line, energy [eV]
        trans_table[i] = ET.iat[i,1] # y line, transmission [-]

    # z integral calculation
    # Ez discretization (don't change the order)
    Ezmin = Ez_table[0]           # [eV], integral lower: change here, depends on condition
    Ez = Ezmin                    # [eV]
    dEz = Ez_table[1]-Ez_table[0] # [eV]

    # flux
    flux_qe = 0 # [W/m2], initialization

    for i in range(0, row):

        # integration for Ez, quantum tunneling flux, Eq (3.10)
        if i==0 or i==row-1:
            flux_qe = flux_qe + total(Ez_table[i], trans_table[i])*(dEz*ucev_to_j)*0.5 # [W/m2]
        else:
            flux_qe = flux_qe + total(Ez_table[i], trans_table[i])*(dEz*ucev_to_j) # [W/m2]

    print("HTC:{:.2f}".format(flux_qe/delta_temperature) + "[W/m2K]") # [W/m2K]

    return flux_qe/delta_temperature

# time display
#elapsed_time = time.time()-start
#print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
