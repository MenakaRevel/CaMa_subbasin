#! /usr/bin/python 
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import sys
import matplotlib.ticker as ticker
import copy
from mpl_toolkits.basemap import Basemap
import matplotlib
# matplotlib.use('Agg')
# import matplotlib.backends as backend
import re
from numpy import ma

# Takashima & Prakat

### Plotting Subbasins
## For Merdians and Parallers (Subbasin Boundaries)
# ----------------------------------------------------------------------------------------------------------------------------------------------------#
def getGridLines (dat):
    lonsize = float (east - west ) / nx
    latsize = float (south - north ) / ny

    # Get meridians
    lons, lats = [], []
    lats_north = np.linspace (north, south, ny + 1 ) [: -1 ]
    for ix in  range (nx - 1 ):
        lats_inter = copy.copy (lats_north)
        lats_inter [dat [:, ix] == dat [:, ix + 1 ]] = np.nan
        lats_this = np.r_ [np.c_ [lats_north, lats_inter].reshape ( -1 ), south]
        lons_this = np.ones ((ny * 2 + 1 )) * (west + (ix + 1 ) * lonsize)
        lons.append (np.r_ [lons_this, np.nan])
        lats.append (np.r_ [lats_this, np.nan])
    meridians = (np.array (lons) .reshape (-1),np.array (lats) .reshape (-1))

    # Get parallels
    lons, lats = [], []
    lons_west = np.linspace (west, east, nx + 1 ) [: -1 ]
    for iy in  range (ny - 1 ):
        lons_inter = copy.copy (lons_west)
        lons_inter [dat [iy,:] == dat [iy + 1 ,: ]] = np.nan
        lons_this = np.r_ [np.c_ [lons_west, lons_inter].reshape ( -1 ), east]
        lats_this = np.ones ((nx * 2 + 1 )) * (north + (iy + 1 ) * latsize)
        lons.append (np.r_ [lons_this, np.nan])
        lats.append (np.r_ [lats_this, np.nan])
    parallels = (np.array (lons) .reshape ( -1 ),np.array (lats) .reshape ( -1 ))

    return meridians, parallels

# ----------------------------------------------------------------------------------------------------------------------------------------------------#
# For plotting the merdians, Paralles and the Colorplot and Boundaries
# ----------------------------------------------------------------------------------------------------------------------------------------------------#

def  draw (subbasin, maxbasin, meridians, parallels,llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,lonint = 5 , latint = 5 , miss = None , figname = None ):
    lonsize = float (east - west ) / nx
    latsize = float (south - north ) / ny
    x0 = int ( max (np.floor ((llcrnrlon-west) / lonsize), 0 ))
    x1 = int ( min (np.ceil ((urcrnrlon-west) / lonsize), nx))
    y0 = int ( max (np.floor ((urcrnrlat-north) / latsize), 0 ))
    y1 = int ( min (np.ceil ((llcrnrlat-north) / latsize), ny))
    west_ = west + lonsize * x0
    east_ = west + lonsize * x1
    north_ = north + latsize * y0
    south_ = north + latsize * y1

    if miss is  None  or miss is 'nan':
        dat_plt = dat [y0: y1, x0: x1]% 30 
    else :
        dat_plt = np.ma.masked_where (dat [y0: y1, x0: x1] == miss,dat [y0: y1, x0: x1]% 30 )                                   ## masking the miss value

    # Basemap with axis as ax , change ax to chnage the plot axis   
    m = Basemap (llcrnrlon = west_, llcrnrlat = south_,urcrnrlon = east_, urcrnrlat = north_, ax = ax)                             
 
    # Drawing fill colours
    im = m.imshow(ma.masked_greater(ma.masked_less(subbasin,0.0),maxbasin),interpolation="nearest",origin="upper",cmap="brg",zorder=100,vmin=1.0,vmax=maxbasin)
    m.colorbar(im,"right",size="2%")

    ##  Drawing Meridian and Parallel (with silver) with labels (Comment if doesnot want to plot)
    m.drawmeridians (np.arange ( -180 , 180 + 1 , lonint), linewidth = 0.3,color = 'grey' , labels = [ 0 , 0 , 0 , 1 ])
    m.drawparallels (np.arange ( -90 , 90 + 1 , latint), linewidth = 0.3,color = 'grey' , labels = [ 1 , 0 , 0 , 0 ])

    # for the colorplot, change cmap for color pattern
    # m.imshow (dat_plt, origin = 'upper' , interpolation = 'nearest' ,cmap = cm.Greens_r, alpha = 0.8)                         

    ## Subbasin Boundaries Plotting
    m.plot (meridians [ 0 ], meridians [ 1 ], linewidth = 0.3, color = 'k',zorder=101 )
    m.plot (parallels [ 0 ], parallels [ 1 ], linewidth = 0.3 , color = 'k' ,zorder=101)

    if figname is not None:
        plt.savefig(figname, bbox_inches='tight',
                    pad_inches=0.1, dpi=300)
    else:
        plt.show()

##Example
fig, ax = plt.subplots(1)
# mapname="glb_15min"
mapname="amz_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
lonint = 10 #40
latint = 10 #20
figname="./img/subbasin_boundry_"+mapname+".png"
# ny, nx = 350 , 500                                  ## Pixel Sizes
# west, east, south, north = -90 , -40 , -25, 10      ## Boundary of data
#----
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()

#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])


# f = '/cluster/data6/menaka/CaMa_subbasin/output/subbasin_reg_06min_srtm.bin'
f ='./output/subbasin_'+mapname+'.bin'
sbasin= np.fromfile(f, np.float32).reshape(ny,nx)       ## Subbasin data
sbasin[(sbasin>=2.000)|(sbasin<1.000)]=-999    ## Masking the subbasins not in amazon

# sbasin[(sbasin>=2.000)|(sbasin<1.000)]=-999    ## Masking the subbasins not in amazon



dat = sbasin     
meridians1, parallels1 = getGridLines (dat)
maxbasin=1.041
im1 = draw(sbasin,maxbasin,meridians1,parallels1,west,south,east,north,lonint=lonint,latint=latint,miss=-9999,figname=figname)  # lonint and latint spacing of meridians and paralles