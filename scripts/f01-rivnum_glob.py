#! /usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib import colors
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import matplotlib as mpl
import sys
import os
import matplotlib.gridspec as gridspec
import string
import calendar
import errno
import re
import math
from numpy import ma 
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
#from slacker import Slacker
from multiprocessing import Pool
from multiprocessing import Process
####################
def mk_dir(sdir):
     try:
         os.makedirs(sdir)
     except:
         pass
####################
mk_dir("../img")
mapname="glb_06min"
nx,ny=(3600,1800)
gsize=0.100
rivnum = "../output/rivnum_"+mapname+".bin"
rivnum = np.fromfile(rivnum,np.int32).reshape(ny,nx)
#--major rivers and Ids
rivid={}
fname="../output/rivmth_"+mapname+".txt"
f = open(fname,"r")
lines = f.readlines()
f.close()
#---
for line in lines[1::]:
  line    = filter(None, re.split(" ",line))
  riverid = int(line[0])
  lon     = int(line[1])
  lat     = int(line[2])
  upar    = float(line[3])
  #print river
  rivid[riverid]=[lon,lat]
#------
boundsrivid=np.arange(0.5,31.5,1.0)
cmapS=cm.jet #_r #Set3_r #Dark2 #
cmapS.set_under("#000000",alpha=0)
cmapS.set_over("#FFFFFF",alpha=0)
cmapS.set_bad("#000000",alpha=0)
normriv=colors.BoundaryNorm(boundsrivid,256)
#--
land="#FFFFFF"#"grey"#
water="#C0C0C0"#"royalblue"#
#--
lllat, urlat, lllon, urlon = -90.0, 90.0, -180.0, 180.0
#-- meridians and parallels
meridians = 40.0
parallels = 20.0
plt.close() 
#--figure in A4 size
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure()#figsize=(11.69,8.27))
G = gridspec.GridSpec(1,1)
#  mtitle="Assimilation Index"
#  fig.suptitle(mtitle,fontsize=12,fontweight="bold")
ax1 = fig.add_subplot(G[0,0])
M  = Basemap(resolution="c", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon,ax=ax1) 
M.drawcoastlines(color="k",linewidth=0.5, zorder=101)
#M.fillcontinents(color=land,lake_color=water,zorder=99)
# M.drawmapboundary(fill_color=water,zorder=98)
#--
M.drawmeridians(np.arange(lllon,urlon+1, meridians), labels=[0, 0, 0, 1], fontsize=10, rotation=0,linewidth=0.5,zorder=104)
M.drawparallels(np.arange(lllat,urlat+1, parallels), labels=[1, 0, 0, 0], fontsize=10,linewidth=0.5,zorder=104)
im=M.imshow(ma.masked_less_equal(rivnum,0.0),interpolation="nearest",cmap=cmapS,origin="upper",alpha=1.0,zorder=100,norm=normriv)
for riverid in np.arange(1,30+1):
     # target pixel
     ix  = rivid[riverid][0]
     iy  = rivid[riverid][1]
     lon = -180.0 + (ix-1)*gsize
     lat = 90.0 - (iy-1)*gsize
     plt.scatter(lon,lat,s=60,marker="o",color="red",zorder=105)
     annotate_string="%d"%(riverid)
     plt.annotate(annotate_string,xy=(lon,lat),xycoords="data",horizontalalignment="left",verticalalignment="top",fontsize=12,zorder=106)
#----------
cb=plt.colorbar(im,orientation="horizontal",ticks=np.arange(1,30.1,5))#,cax=fig.add_axes([0.1,0.07,0.4,0.01]))
#  cb.set_label(label="Number of SWOT observations per cycle",size="10")
plt.title("Basin ID - Global", fontsize=10)
plt.savefig("../img/rivnum_"+mapname+".png",pad_inches=0,bbox_inches="tight",dpi=800)
#plt.show()