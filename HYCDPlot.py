import h5py
import os
import time
import platform
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from matplotlib import colors
import pyresample
import netCDF4 as nc
import os
import glob
# try:
#     exc = None
#     from mpl_toolkits.basemap import Basemap
# except Exception as exc:
#     print(exc)
#     pass
#
# # if exc is not None:
# proj_path = 'C:/Users/zhenjia/AppData/Local/Programs/Python/Python37/' \
#             'Lib/site-packages/osgeo/data/proj/epsg'.replace('/', os.sep)
#
# os.environ['PROJ_LIB'] = os.path.dirname(proj_path)
from mpl_toolkits.basemap import Basemap

## Generate list
os.chdir('H:\gOOLE\server')
path = 'H:\gOOLE\server'
datalist=glob.glob('*L2B*.h5')


def HY1CD_Rrs_Reader(file: str):
    f = h5py.File(file, 'r')
    Rrs = ['Rrs412',
           'Rrs443', 'Rrs490', 'Rrs520', 'Rrs565', 'Rrs670', 'Rrs750', 'Rrs865']
    Rrs_output = {'Rrs412': [],
                  'Rrs4443': [],
                  'Rrs490': [],
                  'Rrs520': [],
                  'Rrs560': [],
                  'Rrs670': [],
                  'Rrs750': [],
                  'Rrs865': [], }
    key = 'Geophysical Data'
    for i in range(0, 8):
        var = Rrs[i]
        sds_obj = f[f'{key}/{var}']
        dn_attrs = sds_obj.attrs
        dn = sds_obj[:]
        slope = dn_attrs['Slope'][0]
        offset = dn_attrs['Offset'][0]
        min_val = dn_attrs['ValidMin'][0]
        max_val = dn_attrs['ValidMax'][0]
        fill_val = dn_attrs['FillValue'][0]
        varvalue = dn * slope + offset
        # Here do not mask minus and zero value of Rrs, you can active it if you want
        varvalue = np.ma.masked_where(varvalue < 0., varvalue)
        key2, var2 = 'Sensor Band Parameters', 'Wavelength'

        wave = np.array(f[f'{key2}/{var2}'])
        Rkey = 'Rrs' + str(wave[i])
        Rrs_output[Rkey] = varvalue
    f.close()
    return Rrs_output

def HY1CD_Var_Reader(file: str):
    f = h5py.File(file, 'r')
    Var = ['CDOM',
           'chl_a', 'Kd490', 'nLw565', 'SSC', 'SST', 'taua865','TSM']
    Var_output = {'CDOM': [],
                  'chl_a': [],
                  'Kd490': [],
                  'nLw565': [],
                  'SSC': [],
                  'SST': [],
                  'taua865': [],
                  'TSM': [], }
    key = 'Geophysical Data'
    for i in range(0, 8):
        var = Var[i]
        sds_obj = f[f'{key}/{var}']
        dn_attrs = sds_obj.attrs
        dn = sds_obj[:]
        slope = dn_attrs['Slope'][0]
        offset = dn_attrs['Offset'][0]
        min_val = dn_attrs['ValidMin'][0]
        max_val = dn_attrs['ValidMax'][0]
        fill_val = dn_attrs['FillValue'][0]
        varvalue = dn * slope + offset
        # Here do not mask minus and zero value of Rrs, you can active it if you want
        varvalue = np.ma.masked_where(varvalue <=0., varvalue)
        # key2, var2 = 'Sensor Band Parameters', 'Wavelength'
        #
        # wave = np.array(f[f'{key2}/{var2}'])
        # Varkey = 'Rrs' + str(wave[i])
        Var_output[var] = varvalue
    f.close()
    return Var_output


def swath_resampling(src_data: np.ma.array, src_lon: np.array, src_lat: np.array,
                     trg_lon: np.array, trg_lat: np.array, search_radius: float):
    if len(trg_lon.shape) == 1:
        grid_def = pyresample.geometry.SwathDefinition(*np.meshgrid(trg_lon, trg_lat))
    else:
        grid_def = pyresample.geometry.SwathDefinition(lons=trg_lon, lats=trg_lat)

    # source grid with original swath data
    swath_def = pyresample.geometry.SwathDefinition(lons=src_lon, lats=src_lat)

    # resample (here we use nearest. Bilinear, gaussian and custom defined methods are available)
    # for more, visit https://pyresample.readthedocs.io/en/latest/
    result = pyresample.kd_tree.resample_nearest(swath_def, src_data, grid_def, epsilon=0.5,
                                                 fill_value=None, radius_of_influence=search_radius)
    return result, grid_def


def plot_geo_image(sds: np.ma.array, lon: np.ndarray, lat: np.ndarray, log10: bool = True, title: str = None,
                   label: str = None,
                   caxis: list = None, lon_range: list = None, lat_range: list = None, save_image: str = None,
                   dpi: int = 100):
    if len(lon.shape) == 1:
        print('MeshGridding...')
        lon, lat = np.meshgrid(lon, lat)

    lon_0 = (lon.min() + lon.max()) / 2
    lat_0 = (lat.min() + lat.max()) / 2



    if (lon_range is not None) and (lat_range is not None):
        m = Basemap(llcrnrlon=min(lon_range), llcrnrlat=min(lat_range),
                    urcrnrlon=max(lon_range), urcrnrlat=max(lat_range),
                    resolution='i', lon_0=lon_0, lat_0=lat_0, projection='tmerc')
    else:
        m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
                    urcrnrlon=lon.max(), urcrnrlat=lat.max(),
                    resolution='i', lon_0=lon_0, lat_0=lat_0, projection='tmerc')
    x2d, y2d = m(lon, lat)

    fig = plt.figure(figsize=(8, 8 * m.aspect))
    ax = fig.add_axes([0.08, 0.1, 0.7, 0.7], facecolor='white')
    # changed to facecolor 8 October 2019

    if (lon_range is not None) and (lat_range is not None):
        parallels = np.linspace(min(lat_range), max(lat_range), 4)
        meridians = np.linspace(min(lon_range), max(lon_range), 3)
    else:
        parallels = meridians = None

    if caxis is not None:
        cmn, cmx = min(caxis), max(caxis)
    else:
        cmn, cmx = sds.min(), sds.max()

    ncl = 150
    if log10 is True:
        norm = colors.LogNorm(vmin=cmn, vmax=cmx)
    else:
        bounds = np.linspace(cmn, cmx, ncl)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncl)

    p = m.pcolor(x2d, y2d, sds, norm=norm, cmap=plt.cm.jet)

    if title is not None:
        plt.title(title)

    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('vertical', size="3%", pad=0.05)
    cax = plt.axes([0.85, 0.1, 0.05, 0.7])  # setup colorbar axes

    cb = plt.colorbar(format='%5.2f', cax=cax)  # draw colorbar
    if label is not None:
        cb.set_label("%s" % label)
    plt.sca(ax)  # make the original axes current again
    plt.clim(cmn, cmx)

    m.drawcoastlines()
    m.drawcountries()

    if save_image is not None:
        plt.savefig(save_image, dpi=dpi, facecolor='w', edgecolor='w', orientation='portrait')
        plt.close()
for file in datalist:
    print(file)
    Var=HY1CD_Var_Reader(file)
    f = h5py.File(file, 'r')
    key = 'Navigation Data'
    var = 'Longitude'
    Lons = np.array(f[f'{key}/{var}'])

    key = 'Navigation Data'
    var = 'Latitude'
    Lats = np.array(f[f'{key}/{var}'])
    x = np.arange(105.9021, 117.6794, .01)  # 1 km grid, full resolution did not work but 500 m was OK!
    y = np.arange(22.1228, 14.9353, -.01)

    var_re, grid = swath_resampling(Var['chl_a'], Lons, Lats, x, y, 5000)

    lon_grid, lat_grid = grid.get_lonlats()
    plot_geo_image(var_re, lon_grid, lat_grid,caxis=[0.001,10], label='Chl-a [mg/m$^{-3}$]', log10=True,save_image=os.path.splitext(file)[0],title=os.path.splitext(file)[0])
    f.close()  # from now on, we no longer need the hdf5 file, so we close it
