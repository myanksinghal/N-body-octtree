import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
from tqdm import tqdm
from joblib import Parallel, delayed
import seaborn as sns
df=pd.read_csv('test_data_file.csv')

save_path='anim/'
N_particles=8200*2
num_updates=int(len(df)/N_particles)
scale=400
mass_scale=1000
KDE_FLAG=False
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
    
    if(KDE_FLAG):
        plt.figure()
        sns.kdeplot(data=red_df, x='x', y='y', hue='mass', fill=False,palette='pastel',levels=7)
        plt.xlim([-scale,scale])
        plt.ylim([-scale,scale])
        plt.savefig(f"{save_path}kde_{num}.png")
        plt.close(fig)
    
    
Parallel(n_jobs=6)(delayed(process)(num) for num in tqdm(range(0,frames)))

