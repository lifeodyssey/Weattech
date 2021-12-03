import h5py
import numpy as np
import matplotlib.pyplot as plt


try:
    exc = None
    from mpl_toolkits.basemap import Basemap
except Exception as exc:
    print(exc)
    pass

import os
import glob
os.chdir('H:\weattech')
path = 'H:\weattech'

datalist=glob.glob('**.nc')
for f1 in datalist:
    f = h5py.File(f1, 'r') # open for read only
    SST=f['SST_ALL']
    x=f['x']
    y=f['y']
    x=np.array(x)
    y=np.array(y)

    SST=np.array(SST)
    SST=np.ma.masked_where(SST<=-5,SST)
    SST=np.ma.masked_where(SST>=45,SST)

    SST=np.ma.masked_where(x<124,SST)
    SST=np.ma.masked_where(SST>=45,SST)
    plt.figure(figsize=(12,12))
    plt.axis('off')
    plt.imshow(SST,cmap='rainbow')
    plt.savefig(os.path.splitext(f1)[0],dpi=300, facecolor='w', edgecolor='w', orientation='portrait')