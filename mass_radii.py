import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import numpy as np
from math import isclose
from astropy.table import Table
from tqdm import tqdm
from scipy.stats import gaussian_kde


df=pd.read_csv('test_data_file.csv')

N_particles=1000
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

def kde_density_center(t,recursions=5):
    X=np.array(t['X']).T
    #print(X.shape)
    Mass=np.array(t['Mass'])
    #print(Mass.shape)
    kde=gaussian_kde(X,weights=Mass)
    minima=X.T.min(axis=0)
    maxima=X.T.max(axis=0)
    i=0
    while(i<=recursions):
     #   print(i)
      #  print(maxima)
       # print(minima)
        space = [np.linspace(mini,maxi,50) for mini, maxi in zip(minima,maxima)]
        grid = np.meshgrid(*space)
        coords = np.vstack(map(np.ravel, grid))
        density = kde(coords)
        inx=np.argmax(density)
        coord=coords.T[inx]
        factor=100/10**(i)
        #print(factor)
        maxima=coord+factor
        minima=coord-factor
        i=i+1
        #print(coords.T[inx])
    return coord


def mass_fraction_radii(test,cod_array,target_mass_frac):
    rads=np.zeros(len(test))
    for inx in tqdm(range(0,len(test))):
        #print(inx)
        t=test[inx]
        total_mass=np.sum(t['Mass'])
        com=cod_array[inx]
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

cod_array=[]
for inx in tqdm(range(0,len(table_array))):
        t=table_array[inx]
        total_mass=np.sum(t['Mass'])
        cod=kde_density_center(t,recursions=5)
        cod_array.append(cod)

#Parallel(n_jobs=6)(delayed(process)(num) for num in tqdm(range(0,frames)))

r_01=mass_fraction_radii(table_array,cod_array,0.1)
r_05=mass_fraction_radii(table_array,cod_array,0.5)
r_025=mass_fraction_radii(table_array,cod_array,0.25)
r_075=mass_fraction_radii(table_array,cod_array,0.75)

plt.plot(time_array,r_01,label="10%")
plt.plot(time_array,r_025,label="25%")
plt.plot(time_array,r_05,label="50%")
plt.plot(time_array,r_075,label="75%")
plt.yscale('log')
plt.legend()

plt.savefig("mass_radii.png")
