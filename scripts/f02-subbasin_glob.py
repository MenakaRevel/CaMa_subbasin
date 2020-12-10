#!/opt/local/bin/pyuhon
#-*- coding: utf-8 -*-
 
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
from mpl_toolkits.basemap import Basemap
import os
from numpy import ma
import re
#---------
# CaMa_v394
#CaMa_dir="/work/a02/menaka/CaMa_package_v394_20190508/CaMa_package_v394_20190508/"
#CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/work/a01/modi/cama_v395"
#mapname="glb_06min"
mapname="reg_06min_srtm"
#nx,ny=(3600,1800)
#gsize=0.100
nextxy=CaMa_dir+"/map/"+mapname+"/nextxy.bin"
#----
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])
# #nextxy= "nextxy.bin"
# nextxy=np.fromfile(nextxy,np.int32).reshape(2,720,1440)
# nextx_v394 = nextxy[0]
# nexty_v394 = nextxy[1]
# # CaMa_v362
# nextxy="/work/a02/menaka/CaMa_etc/nextxy.bin"
# nextxy=np.fromfile(nextxy,np.int32).reshape(2,720,1440)
# nextx_v362 = nextxy[0]
# nexty_v362 = nextxy[1]
# subbasin
subbasin="../output/subbasin_"+mapname+".bin"
#subbasin=np.fromfile(subbasin,np.float32).reshape(ny,nx)#*1.0e3
subbasin=np.fromfile(subbasin,np.int32).reshape(ny,nx)

subcolor="../output/subbasin_"+mapname+".bin"
#subbasin=np.fromfile(subbasin,np.float32).reshape(ny,nx)#*1.0e3
subcolor=np.fromfile(subcolor,np.int32).reshape(ny,nx)
#--
# pos_v362=ma.masked_where(nextx_v362==-9999, nexty_v362*1440 + nextx_v362).filled(-9999)
# pos_v394=ma.masked_where(nextx_v394==-9999, nexty_v394*1440 + nextx_v394).filled(-9999)
#---------
land="#FFFFFF"
water="#C0C0C0"

lllat= south #-90.0
urlat= north # 90.0
lllon= west #-180.0
urlon= east # 180.0

londiff=(urlon-lllon)*int(1.0/gsize)
latdiff=(urlat-lllat)*int(1.0/gsize)

npix=int((90-urlat)*int(1.0/gsize))
spix=int((90-lllat)*int(1.0/gsize))
wpix=int((180+lllon)*int(1.0/gsize))
epix=int((180+urlon)*int(1.0/gsize))

#--subbasin----
plt.close()
cmap=cm.rainbow_r #viridis_r
cmap.set_under("w",alpha=0)
maxbasin=1.200
resol=1
plt.figure()#figsize=(7*resol,3*resol))
m = Basemap(projection='cyl',llcrnrlat=lllat+0.1,urcrnrlat=urlat-0.1,llcrnrlon=lllon-0.1,urcrnrlon=urlon+0.1, lat_ts=0,resolution='c')
#m.drawcoastlines( linewidth=0.3, color='k' )
m.fillcontinents(color=land,lake_color=water)
m.drawmapboundary(fill_color=water)
m.drawparallels(np.arange(lllat,urlat+0.1,20), labels = [1,0,0,0], fontsize=8,linewidth=0.1)
m.drawmeridians(np.arange(lllon,urlon+0.1,40), labels = [0,0,0,1], fontsize=8,linewidth=0.1)
#im = m.imshow(np.flipud(ratio),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
im = m.imshow(ma.masked_greater(ma.masked_less(subbasin,0.0),maxbasin),interpolation="nearest",origin="upper",cmap=cmap,zorder=100)#vmin=0,vmax=1440*720
m.colorbar(im,"right",size="2%")
plt.title("Subbasin "+mapname)#"nextx_v362")
figname="../img/subbasin_"+mapname+".png"#"nextx_v362.png"
plt.savefig(figname,dpi=800,bbox_inches="tight", pad_inches=0)

#plt.show()

##---pos_v394----
#plt.close()
#cmap=cm.viridis_r
#cmap.set_under("w",alpha=0)
#
#resol=1
#plt.figure()#figsize=(7*resol,3*resol))
#m = Basemap(projection='cyl',llcrnrlat=lllat+0.1,urcrnrlat=urlat-0.1,llcrnrlon=lllon-0.1,urcrnrlon=urlon+0.1, lat_ts=0,resolution='c')
##m.drawcoastlines( linewidth=0.3, color='k' )
#m.fillcontinents(color=land,lake_color=water)
#m.drawmapboundary(fill_color=water)
#m.drawparallels(np.arange(lllat,urlat+0.1,20), labels = [1,0,0,0], fontsize=8,linewidth=0.1)
#m.drawmeridians(np.arange(lllon,urlon+0.1,40), labels = [0,0,0,1], fontsize=8,linewidth=0.1)
##im = m.imshow(np.flipud(ratio),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
#im = m.imshow(ma.masked_less(pos_v394,0),vmin=0,vmax=1440*720,interpolation="nearest",origin="upper",cmap=cmap,zorder=100)
#m.colorbar(im,"right",size="2%")
#plt.title("pos_v394")#"nextx_v394")
#figname="pos_v394.png"#"nextx_v394.png"
#plt.savefig(figname,dpi=800,bbox_inches="tight", pad_inches=0)
#
#plt.show()
#
###---diff v394 v362----
#plt.close()
#cmap=cm.bwr
##cmap.set_under("w",alpha=0)
#
#resol=1
#plt.figure()#figsize=(7*resol,3*resol))
#m = Basemap(projection='cyl',llcrnrlat=lllat+0.1,urcrnrlat=urlat-0.1,llcrnrlon=lllon-0.1,urcrnrlon=urlon+0.1, lat_ts=0,resolution='c')
##m.drawcoastlines( linewidth=0.3, color='k' )
#m.fillcontinents(color=land,lake_color=water)
#m.drawmapboundary(fill_color=water)
#m.drawparallels(np.arange(lllat,urlat+0.1,20), labels = [1,0,0,0], fontsize=8,linewidth=0.1)
#m.drawmeridians(np.arange(lllon,urlon+0.1,40), labels = [0,0,0,1], fontsize=8,linewidth=0.1)
##im = m.imshow(np.flipud(ratio),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
#im = m.imshow(pos_v394-pos_v362,interpolation="nearest",origin="upper",cmap=cmap,zorder=100)
#m.colorbar(im,"right",size="2%")
#plt.title("diff_pos")#"nextx_v394")
#figname="diff_pos.png"#"nextx_v394.png"
#plt.savefig(figname,dpi=800,bbox_inches="tight", pad_inches=0)
#
#plt.show()