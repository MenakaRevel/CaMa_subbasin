#!/opt/local/bin/python
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
from mpl_toolkits.basemap import Basemap
import os
import re
from numpy import ma
import random
####################
#-----------------------------
def mk_dir(sdir):
     try:
         os.makedirs(sdir)
     except:
         pass
#-----------------------------
def latlon_river(rivername):
     global lllat, urlat, lllon, urlon
     if rivername=="LENA":
         lllat = 50.
         urlat = 80.
         lllon = 100.
         urlon = 145.
     if rivername=="NIGER":
         lllat = 0.
         urlat = 25.
         lllon = -15.
         urlon = 20.
     if rivername=="AMAZONAS":
         lllat = -20.
         urlat = 10.
         lllon = -80.
         urlon = -45.
     if rivername=="MEKONG":
         lllat = 10.
         urlat = 35.
         lllon = 90.
         urlon = 110.
     if rivername=="MISSISSIPPI":
         lllat = 25.
         urlat = 50.
         lllon = -115.
         urlon = -75.
     if rivername=="OB":
         lllat = 40.
         urlat = 70.
         lllon = 55.
         urlon = 95.
     if rivername=="CONGO":
         lllat = -15.
         urlat = 10.
         lllon = 10.
         urlon = 35.
     if rivername=="INDUS":
         lllat = 20.
         urlat = 40.
         lllon = 60.
         urlon = 80.
     return lllat, urlat, lllon, urlon
#-------------------------------------
def riveridname(rivername):
     river=rivername[0]+rivername[1::].lower()
     if rivername=="AMAZONAS":
         river="Amazon"
     if rivername=="ST._LAWRENCE":
         river="St Lawrence"
     if rivername=="BRAHMAPUTRA":
         river="Ganges-Brahamaputra"
     return river
#---------
mk_dir("./img")
#-------
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# mapname="glb_06min"
# mapname="glb_15min"
mapname="amz_06min"
# CaMa_dir="/work/a01/modi/cama_v395"
# mapname="reg_06min_srtm"
# nx,ny=(3600,1800)
# gsize=0.100
mk_dir("./img/subbasin_"+mapname)
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
# subbasin
subbasin="./output/subbasin_"+mapname+".bin"
subbasin=np.fromfile(subbasin,np.float32).reshape(ny,nx)#*1.0e3
# rivnum
rivnum = "./output/rivnum_"+mapname+".bin"
rivnum = np.fromfile(rivnum,np.int32).reshape(ny,nx)
#---------
w=0.15#*2
alpha=1
alpha1=1
width=0.5

land="#FFFFFF"
water="#C0C0C0"

lllat= -90.0
urlat=  90.0
lllon= -180.0
urlon=  180.0

londiff=(urlon-lllon)*4
latdiff=(urlat-lllat)*4
#********************
#major rivers and Idsi
rivid={}
fname="./output/river30_id.txt"
f = open(fname,"r")
lines = f.readlines()
f.close()
#---
for line in lines:
     line    = filter(None, re.split(",",line))
     riverid = int(line[0])
     river   = filter(None, re.split("\n",line[1]))[0].strip()
     #print river
     rivid[river]=riverid
#---------
# rivernames=["AMAZONAS"]
rivernames=["MISSISSIPPI"]
# rivernames=["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS"]
for rivername in rivernames:
     data=np.zeros([ny,nx],np.float32)
     c_nextx=(rivnum==rivid[riveridname(rivername)])*1.0
     latlon_river(rivername)
     data=subbasin*((rivnum==rivid[riveridname(rivername)])*1.0)
     data=np.array(data)
     #print np.amax(data)
     subs=int((np.amax(data)-rivid[riveridname(rivername)])*1e3) + 1
     #rd=(np.random.normal(500,,int(subs))+500).astype(int)
     #rd=random.randrange(0,1000)i
     #print subs
     rd=np.random.randint(0,1000,size=(int(subs)))
     #print len(rd)
     for i in np.arange(0,subs):
         basin=round(rivid[riveridname(rivername)])*1.0+int(i)*0.001
         index=np.where(subbasin==basin)
         #print basin, index[0], index[1]
         data[index]=rd[int(i)]
         #print rd[int(i)]
    #------------
     npix=int((90-urlat)*int(1.0/gsize))
     spix=int((90-lllat)*int(1.0/gsize))
     wpix=int((180+lllon)*int(1.0/gsize))
     epix=int((180+urlon)*int(1.0/gsize))

     #--subbasin----
     plt.close()
     cmap=cm.nipy_spectral #jet_r #rainbow_r #viridis_r
     cmap.set_under("w",alpha=0)

     resol=1
     plt.figure()#figsize=(7*resol,3*resol))
     M = Basemap(projection='cyl',llcrnrlat=lllat+0.1,urcrnrlat=urlat-0.1,llcrnrlon=lllon-0.1,urcrnrlon=urlon+0.1, lat_ts=0,resolution='c')
     #m.drawcoastlines( linewidth=0.3, color='k' )
     M.fillcontinents(color=land,lake_color=water)
     M.drawmapboundary(fill_color=water)
     M.drawparallels(np.arange(lllat,urlat+0.1,5.), labels = [1,0,0,0], fontsize=8,linewidth=0.1)
     M.drawmeridians(np.arange(lllon,urlon+0.1,5.), labels = [0,0,0,1], fontsize=8,linewidth=0.1)
     #im = M.imshow(np.flipud(ratio),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
     #im = M.imshow(ma.masked_greater(ma.masked_less(subbasin,1.000),1.999),interpolation="nearest",origin="upper",cmap=cmap,zorder=100)#vmin=0,vmax=1440*720
     im = M.imshow(ma.masked_less_equal(data[npix:spix,wpix:epix],0.0),interpolation="nearest",origin="upper",cmap=cmap,zorder=100)#
     #--
     box="%f %f %f %f"%(lllon,urlon,urlat,lllat) 
     #  os.system("./bin/txt_vector "+str(lllon)+str(urlon)+str(urlat)+str(lllat)+" > tmp.txt")
     os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > tmp.txt") 
     for LEVEL in range(3,10+1):
         os.system("./bin/print_rivvec tmp.txt 1 "+str(LEVEL)+" > tmp2.txt")
         width=float(LEVEL)*w
         #print width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
         # open tmp2.txt
         f = open("tmp2.txt","r")
         lines = f.readlines()
         f.close() 
         #--- 
         for line in lines:
             line    = filter(None, re.split(" ",line))
             lon1 = float(line[0])
             lat1 = float(line[1])
             lon2 = float(line[3]) 
             lat2 = float(line[4])

             iix = int((lon1 - west)*int(1.0/gsize)) 
             iiy = int((-lat1 + north)*int(1.0/gsize))

             if c_nextx[iiy,iix] <= 0:
                 continue
                  #print lon1,lat1,width
             x1,y1=M(lon1,lat1)
             x2,y2=M(lon2,lat2)
             M.plot([x1,x2],[y1,y2],color="#C0C0C0",linewidth=width,zorder=101,alpha=alpha)
     #M.colorbar(im,"right",size="2%")
     plt.title(rivername)
     figname="./img/subbasin_"+mapname+"/"+rivername+".png"
     print figname
     plt.savefig(figname,dpi=800,bbox_inches="tight", pad_inches=0)
os.system("rm -r tmp*.txt")