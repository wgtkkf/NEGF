# Coded by Takuro TOKUNAGA
# Last modified: November 03 2020

import numpy as np                             # not used here
import time                                    # used
import pandas as pd                            # used

start = time.time()

def convert(arg_bias, arg_filenum):

    EB = pd.read_csv("../green_d/barrier/bprofile_"+str(arg_filenum)+".txt", sep=' ', header=None)
    EB.columns = ["E", "B"] # Energy & Barrier
    row, col = EB.shape # row & column of matorix

    x = np.zeros(row, dtype='float64')
    y = np.zeros(row, dtype='float64')

    f1 = open("../green_d/barrier_r/bprofile_r_"+str(arg_filenum)+".txt", 'w')

    for i in range(0, row):
        x[i] = EB.iat[i,0] # x line
        y[i] = EB.iat[i,1] # y line

        if arg_bias == 0: # no bias

            # 0 [eV] both left and right side
            if y[i] < 0:
                y[i] = 0
            else:
                y[i] = y[i]

        elif arg_bias != 0: # with bias
            if i <=int(row*0.5) and y[i] < 0:
                y[i] = 0
            elif i > int(row*0.5) and y[i] <=arg_bias:
                y[i] = arg_bias
            else:
                y[i] = y[i]

        f1.write(str(x[i])) # change of variable
        f1.write(str(' '))
        f1.write(str(y[i]))
        f1.write('\n')

    # file close
    f1.close()

# time display
#elapsed_time = time.time()-start
#print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
