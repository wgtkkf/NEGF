# Electron tunneling transmission & flux by Non-equilibrium Green's function
# Coded by Takuro TOKUNAGA
# Last modified: November 16 2020

import numpy as np
import time
import sys
sys.path.append('../green_d/')

from barrier  import wz_total_d       # import function
from barrier_convert  import convert  # import function
from trans  import transNEGF          # import function
from flux  import fluxNEGF          # import function

start = time.time()

# Unit conversion:
ucev_to_j = 1.602176620898*np.power(10.,-19) # electron volt to joule
ucjto_ev = 1/ucev_to_j                       # joule to electron volt
ucnano = 1.0*np.power(10.,-9)                # nm to m
ucangs = 1.0*np.power(10.,-10)               # angstrom to m

# physical constants
kb = 1.38064852*np.power(10.,-23)            # Boltzmann constant
rh = 6.62607004*np.power(10.,-34)/(2*np.pi)  # reduced Planck constant
c = 299792458                                # light speed
eps_0 = 8.854187817*np.power(10.,-12)        # [F/m]
sq = 1.602176620898*np.power(10.,-19)        # [C] elementary cahrge, small q
me = 9.10938356*np.power(10.,-31)            # electron mass, [kg]
ac = 6.02214076*np.power(10.,23)             # Avogadro constant, [mol-1]
counter = 0                                  # [-]

### parameters from here
## memo: change only gapmin, gapmax, and volt if necessary
# gap
small = 0.01*ucnano           # zero [m]
gapmin = 0.2*ucnano           # [m]
gap = gapmin                  # [m]
gapmax = 1.0*ucnano-small     # [m]

# bias
volt = -0.600    # [V], load (vias) voltage
#volt = 0.000    # [V], load (vias) voltage
se = 1         # [eV/V], single electron
volt = se*volt # [eV]
### parameters until here

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# main start
begin()

# file open
f1 = open('flux.txt', 'w')
f2 = open('../green_d/barrier/info.txt', 'w')  # write mode

f2.write('filenumber[-] gapmax[nm]')
f2.write('\n')

# barrier profile generation
while gap < gapmax+(2*small): # [m]

    wz_total_d(small, gap-small, counter)

    # dgap update
    if gapmax<2.0*ucnano:
        dgap = 0.1*ucnano
    elif gapmax>=2.0*ucnano:
        dgap = 1.0*ucnano

    # file output
    f2.write(str(counter)) # [-]
    f2.write(str(' '))
    f2.write(str("{:.2f}".format(gap/ucnano))) # [nm]
    f2.write('\n')

    # gap update
    gap = gap + dgap
    counter += 1

for i in range(0, counter):

    # barrier conversion
    convert(volt, i)

    # transmission
    transNEGF(i)

    # flux
    HTC = fluxNEGF(i)

    # output
    f1.write(str(HTC)) # [W/m2K]
    f1.write('\n')

f1.close
f2.close

# end
end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
