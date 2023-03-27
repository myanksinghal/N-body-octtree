import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import numpy as np
from math import isclose
from astropy.table import Table
from tqdm import tqdm


df=pd.read_csv('test_data_file.csv')

N_particles=8200
time_array=[]
table_array=[]
number_timesteps=int(len(df)/N_particles)
for inx in tqdm(range(0,number_timesteps)):
    t=df[inx*N_particles:(inx+1)*N_particles]
    #print(t)
    time_array.append(t.iloc[0]['time'])
    #temp_pos_array=[]
    #mass_array=[]
    X=np.array([np.array(t['x']),np.array(t['y']),np.array(t['z'])]).transpose()
    #print(X.shape) 
    Mass=np.array(t['mass'])
    #print(Mass.shape)
    tf=Table(data=[Mass,X],names=['Mass','X'])
    table_array.append(tf)

print("Data Loaded\n")

def mass_fraction_radii(test,target_mass_frac):
    rads=np.zeros(len(test))
    for inx in tqdm(range(0,len(test))):
        #print(inx)
        t=test[inx]
        total_mass=np.sum(t['Mass'])
        com=np.zeros((3))
        for i in range(0,3):
            com[i]=np.sum(t['Mass']*t['X'][:,i])/total_mass
        com_X=t['X']-com
        com_len=np.linalg.norm(com_X,axis=1)
        rad=1
        while(True):
            temp_red=t[com_len<rad]
            mass_frac=np.sum(temp_red['Mass'])/total_mass
            #print(mass_frac/target_mass_frac)
            if (mass_frac/target_mass_frac<1.15 and mass_frac/target_mass_frac>0.985):
                break
            if mass_frac<target_mass_frac:
                rad=rad+0.01
            elif mass_frac>target_mass_frac:
                rad=rad-0.001
        rads[inx]=rad
    return rads



r_05=mass_fraction_radii(table_array,0.5)
r_025=mass_fraction_radii(table_array,0.25)
r_075=mass_fraction_radii(table_array,0.75)

plt.plot(time_array,r_025)
plt.plot(time_array,r_05)
plt.plot(time_array,r_075)
plt.yscale('log')

plt.savefig("mass_radii.png")
