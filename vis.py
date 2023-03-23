import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
from tqdm import tqdm
from joblib import Parallel, delayed

df=pd.read_csv('test_data_file.csv')

save_path='anim/'
N_particles=1000
num_updates=168
scale=1
mass_scale=1000
frames=num_updates;
def process(num):
    plt.style.use(['dark_background'])
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    del_t=0.1
    red_df=df.iloc[num*N_particles:(num+1)*N_particles]
    ax.clear()
    ax.scatter(red_df.x, red_df.y, red_df.z, s=red_df.mass*mass_scale,c='white',marker='o',alpha=0.5)
    ax.set_xlim(-scale,scale)
    ax.set_ylim(-scale,scale)
    ax.set_zlim(-scale,scale)
    plt.axis('off')

    fig.savefig(f"{save_path}{num}.png",dpi=600)
    plt.close(fig)
    
    
Parallel(n_jobs=6)(delayed(process)(num) for num in tqdm(range(0,frames)))

