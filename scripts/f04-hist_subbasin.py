#! usr/python
import numpy  as np
import os
import errno
import datetime
import shutil
import numpy.random as rd
from multiprocessing import Pool
from multiprocessing import Process
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
#---
mapname="glb_06min"
fname="../output/outsub_"+mapname+".txt"
f = open(fname,"r")
lines = f.readlines()
f.close()
uparea=[]
print lines[0]
for line in lines[1::]:
    line    = filter(None, re.split(" ",line))
    #print line[0] #float(line[0])
    ID = float(line[0])
    if ID >= 600.0:
        continue
    print (ID, float(line[3]))
    uparea.append(float(line[3]))
#--
uparea=np.array(uparea)
#--
fig=plt.figure()
G = gridspec.GridSpec(1,1)
ax=fig.add_subplot(G[0,0])
#num_bins=100
#n, bins, patches = ax.hist(uparea, num_bins, facecolor='blue', alpha=0.5)
sns.distplot(uparea, hist=True, color="orange")
ax.set_xlim(xmin=0.0,xmax=1.0e12)#np.amax(uparea))
plt.savefig("../img/hist_uparea_"+mapname+".png")