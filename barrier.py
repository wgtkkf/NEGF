# Electron tunneling transmission by Non-equilibrium Green's function
# T. L. Westover and T. S. Fisher Phys. Rev. B. 77 115426 (2008)
# Coded by Takuro TOKUNAGA
# Last modified: January 25 2020

import numpy as np
import time
from scipy import special, optimize # used in the wz_total

start = time.time()

# Unit conversion:
ucev_to_j = 1.602176620898*np.power(10.,-19) # electron volt to joule
ucjto_ev = 1/ucev_to_j                       # joule to electron volt
ucnano = 1.0*np.power(10.,-9)                # nm to m
ucangs = 1.0*np.power(10.,-10)               # angstrom to m

# physical constants
kb = 1.38064852*np.power(10.,-23)           # Boltzmann constant
rh = 6.62607004*np.power(10.,-34)/(2*np.pi) # reduced Planck constant
c = 299792458                               # light speed
eps_0 = 8.854187817*np.power(10.,-12)       # [F/m]
sq = 1.602176620898*np.power(10.,-19)       # [C] elementary cahrge, small q
me = 9.10938356*np.power(10.,-31)           # electron mass, [kg]
ac = 6.02214076*np.power(10.,23)            # Avogadro constant, [mol-1]

# work function & voltage
phi_e = 5.10    # [eV], emitter work function 5.8
phi_c = 5.10    # [eV], collector work function 4.7
volt = -0.600    # [V], load (vias) voltage
#volt = 0.000    # [V], load (vias) voltage
se = 1         # [eV/V], single electron
volt = se*volt # [eV]
Ef = 5.51      # [eV], gold fermi energy

# gap
small = 0.01*ucnano           # zero [m]
gapmin = small                # [m]
number = 1000                 # [-]

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# total profile
def wz_total(arg_x, arg_d):
    # ideal barrier profile
    ibp = phi_e-(phi_e-phi_c-volt)*(arg_x/arg_d) # [eV]

    # space charge
    scbp = 0 # in case of nanosacle, Devon's Ph.D. thesis

    # image charge barrier profile
    term1 = np.power(sq,2.0)/(16*np.pi*eps_0*arg_d) # [J]
    term1 = term1/ucev_to_j # [eV]
    term2 = -2*special.digamma(1) # [-]
    term3 = special.digamma(arg_x/arg_d) # [-]
    term4 = special.digamma(1-arg_x/arg_d) # [-]
    icbp = term1*(term2+term3+term4) # [eV]
    #icbp = 0

    # Fermi energy is the reference
    total = ibp+scbp+icbp

    # Fermi energy
    #total = total + Ef

    return total # [eV]

def wz_total_d(arg_gap, arg_dmax, arg_counter):
    # tables for output
    gap_table = np.zeros(number+1, dtype='float64')
    Ez_gap_table = np.zeros(number+1, dtype='float64')
    counter = 0
    dgap = (arg_dmax-arg_gap)/number # [m]

    # output
    f1 = open("../green_d/barrier/bprofile_"+str(arg_counter)+".txt", 'w')  # write mode

    while arg_gap < arg_dmax: # [m]
        wz = wz_total(arg_gap, arg_dmax) # [eV]

        # for graph & green's function
        gap_table[counter] = arg_gap/ucnano # [nm]
        Ez_gap_table[counter] = wz # [eV]

        f1.write(str(gap_table[counter])) # [nm]
        f1.write(str(' '))
        f1.write(str(Ez_gap_table[counter])) # [eV]
        f1.write('\n')

        # gap update
        arg_gap = arg_gap + dgap
        counter = counter + 1

    # file close
    f1.close

    return 0

# time display
#elapsed_time = time.time()-start
#print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
