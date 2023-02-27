import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
from tqdm import tqdm
from joblib import Parallel, delayed

df=pd.read_csv('test_data_file.csv')

save_path='anim/'
N_particles=10000
num_updates=100
frames=num_updates;
def process(num):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    del_t=0.1
    red_df=df.iloc[num*N_particles:(num+1)*N_particles]
    ax.clear()
    for i in range(0,len(red_df)):
        ax.scatter(red_df.iloc[i].x, red_df.iloc[i].y, red_df.iloc[i].z, s=red_df.iloc[i].mass)
    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    ax.set_zlim(-100,100)
    fig.savefig(f"{save_path}{num}.png")
    plt.close(fig)
    
    
Parallel(n_jobs=10)(delayed(process)(num) for num in tqdm(range(0,frames)))

