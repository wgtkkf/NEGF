# Coded by Takuro TOKUNAGA
# Last modified: March 11, 2020
# Updated: November 17, 2020

import time
import os
import glob

start = time.time()

def comments():
    print ("Removing old files..")

def begin():
    print ("begin")

def end():
    print ("end")

path = "/Users/Takuro/codes/enfht/electron/green_d/"

# main

comments()
begin()

for file in glob.glob(path + 'barrier/*.txt', recursive=True):
    os.remove(file)

for file in glob.glob(path + 'barrier_r/*.txt', recursive=True):
    os.remove(file)

for file in glob.glob(path + 'transmission/*.txt', recursive=True):
    os.remove(file)

for file in glob.glob(path + '/flux.txt', recursive=True):
    os.remove(file)

# end
end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
