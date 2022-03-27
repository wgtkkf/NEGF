# Electron tunneling transmission by Non-equilibrium Green's function
# T. L. Westover and T. S. Fisher Phys. Rev. B. 77 115426 (2008)
# Coded by Takuro TOKUNAGA
# Last modified: January 25 2020

import numpy as np
import pandas as pd
from numpy import linalg as LA
import time

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
U1 = 0         # [eV], left, fermi energy is referenced
Un = 0         # [eV], right, fermi energy is referenced
number = 1000  # [-]

# parameters: Green's function matrix related
mass = me                         # [kg], electron mass
lc = 4.065*ucangs*0.5             # [m]
st = np.power(rh/lc,2.0)/(2*mass) # small t, hopping
st = st*ucjto_ev                  # [eV]

def transNEGF(arg_filenum):
    # output file by barrier.py open
    barrier_table = pd.read_csv("../green_d/barrier_r/bprofile_r_"+str(arg_filenum)+".txt", sep=' ', header=None)
    row,col = barrier_table.shape

    # initialization
    gap_table = np.zeros(row, dtype='float64')
    Ez_gap_table = np.zeros(row, dtype='float64')
    trans_table = np.zeros(number, dtype='float64')
    Ez_table = np.zeros(number, dtype='float64')

    # Green's function matrix related
    Gamma1=np.zeros((row,row), dtype=np.complex) # eq. 10, size:number*number
    Gamma2=np.zeros((row,row), dtype=np.complex) # eq. 10, size:number*number
    Gd=np.zeros((row,row), dtype=np.complex)     # eq. 11, size:number*number
    imatrix=np.identity(row, dtype=np.complex)  # size:number*number
    # Hl, Sigma1, 2: Defined in the energy loop
    delta = 1.0*np.power(10.,-6.0)

    for i in range(0, row):
        gap_table[i] = barrier_table.iat[i,0]
        Ez_gap_table[i] = barrier_table.iat[i,1]

    Ezmin = 0  # Fermi level
    Ez = Ezmin
    Ezmax = max(Ez_gap_table)
    dEz = (Ezmax-Ezmin)/number
    #print(Ezmin)
    #print(Ezmax)
    #print(dEz)

    # file open
    f1 = open("../green_d/transmission/trans_"+str(arg_filenum)+".txt", 'w')

    for i in range(0, number):
        # matrix initialization
        Hlongitudinal=np.zeros((row,row), dtype=np.complex)  # eq. 11, size:number*number
        Sigma1=np.zeros((row,row), dtype=np.complex) # eq. 11, size:number*number
        Sigma2=np.zeros((row,row), dtype=np.complex) # eq. 11, size:number*number

        # Hlongitudinal construction
        for j in range(0,row): # 0 ~ row-1
            # diagonal components
            Hlongitudinal[j][j] = (2*st + Ez_gap_table[j])# - Ez # [eV]
            if j < row-1:
                # tridiagonal components
                Hlongitudinal[j][j+1] = -st # [eV]
                Hlongitudinal[j+1][j] = -st # [eV]

        #print(Hlongitudinal)

        # Sigma1, only non zero: (0,0)
        argument1 = 1-(Ez-U1)/(2*st) # [eV/eV]
        #print(argument1)
        if argument1 > 1 or argument1 < -1:
            component1 = 0
        else:
            term1 = np.arccos(argument1) # k1*a
            component1 = -st*np.exp(1j*term1) # k1*a, [eV]

        Sigma1[0][0] = component1 # [eV]
        #print(component1)

        # Sigma2, only non zero: (number,number)
        argument2 = 1-(Ez-Un)/(2*st) # [eV/eV]
        if argument2 > 1 or argument2 < -1:
            component2 = 0
        else:
            term2 = np.arccos(argument2)
            component2 = -st*np.exp(1j*term2) # kn*a

        Sigma2[row-1][row-1] = component2 # kn*a, [eV]

        #print(Sigma1[0][0])
        #print(str(Sigma2[row-1][row-1]))

        # Gamma1
        Gamma1 = 1j*(Sigma1-np.matrix.getH(Sigma1)) # eq.12
        # Gamma2
        Gamma2 = 1j*(Sigma2-np.matrix.getH(Sigma2)) # eq.12

        # Green's function
        Gd_temp = ((Ez+1j*delta)*imatrix - Hlongitudinal - Sigma1 - Sigma2)
        Gd = np.linalg.inv(Gd_temp) # eq.11

        # transmission, eq.10
        trans_temp = LA.multi_dot([Gamma1,Gd,Gamma2,np.matrix.getH(Gd)])
        trans = np.matrix.trace(trans_temp)
        trans_table[i] = trans.real

        Ez_table[i] = Ez

        # output
        f1.write(str(Ez_table[i])) # [-]
        f1.write(str(' '))
        f1.write(str(trans_table[i])) # [-]
        f1.write('\n')

        #print(Ez_table[i])
        Ez = Ez + dEz

    # close
    f1.close

    return 0

# time display
#elapsed_time = time.time()-start
#print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
