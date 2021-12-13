import h5py
import numpy as np
import matplotlib.pyplot as plt

import os
import glob
os.chdir('H:\weattech')
path = 'H:\weattech'

datalist=glob.glob('**.nc')
# for f1 in datalist:
#     f = h5py.File(f1, 'r') # open for read only
#     SST=f['SST_ALL']
#     x=f['x']
#     y=f['y']
#     x=np.array(x)
#     y=np.array(y)
#
#     SST=np.array(SST)
#     SST=np.ma.masked_where(SST<=-5,SST)
#     SST=np.ma.masked_where(SST>=45,SST)
#
#     plt.figure(figsize=(12,12))
#     plt.axis('off')
#     #plt.imshow(SST,cmap='rainbow')
#     plt.imshow(SST[300:1218, 1000:2000], cmap='rainbow')
#     plt.savefig(os.path.splitext(f1)[0],dpi=300, facecolor='w', edgecolor='w', orientation='portrait')
def fy4disk(rawfile,dim):
    sz = np.fromfile(rawfile, dtype=float, count=dim*dim*2)
    latlon = np.reshape(sz,(dim,dim,2))

    lat = latlon[:,:,0]
    lon = latlon[:,:,1]

    lat[lat > 100] = np.nan
    lon[lon < 0  ] = lon[lon < 0  ] + 360.
    lon[lon > 361] = np.nan

    return lon, lat

rawfile = r'H:\gOOLE\FullMask_Grid_4000\FullMask_Grid_4000.raw'


dim = 2748 # 4km

    # Two ways of reading an binary file
# lon, lat = fy4raw(rawfile,dim)
lon, lat = fy4disk(rawfile,dim)
import cartopy.crs as ccrs
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
    fig = plt.figure(figsize=(16, 16))  # 设置画布
    plt.axis('off')
    proj = ccrs.PlateCarree()  # 创建一个投影
    ax = plt.axes(projection=proj)  # 创建一个画纸， 并指明投影类型
    extent = [70, 140, 0, 60]  # 显示范围
    ax.set_extent(extent, proj)
    cs = ax.contourf(lon, lat, SST, transform=proj, cmap='rainbow')
    ax.background_patch.set_visible(False)  # Background
    ax.outline_patch.set_visible(False)  # Borders
    plt.savefig(os.path.splitext(f1)[0],dpi=300, facecolor='w', edgecolor='w', orientation='portrait',transparent=True)