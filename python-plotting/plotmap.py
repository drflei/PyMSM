import sys 
#sys.path.append("/home/flei/workspace/REST-SIM/extra_libs/C++/spenvis_csv")

#from SpenvisCSVFileHandler import SpenvisCSVFileHandler

#Read map file
###############

if len(sys.argv) is not 2 :  
    print "Wrong arguments specified."
    sys.exit(1)

file = sys.argv[1]
#hist = sys.argv[2]

from mpl_toolkits.basemap import Basemap
from matplotlib.mlab import griddata

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
# make up data.
xi=[]; yi=[]; zi=[]
nlon = 73
nlat = 37
ndata = nlon*nlat

with open(file, "r") as f:
    for line in f:
        entries = line.split()
        xi.append(float(entries[1]))
        yi.append(float(entries[0]))
        zi.append(float(entries[4]))
       
# define grid.
xi = np.array(xi).reshape(nlat,nlon)
yi = np.array(yi).reshape(nlat,nlon)
zi = np.array(zi).reshape(nlat,nlon)
xi=xi.transpose()
yi=yi.transpose()
zi=zi.transpose()

#print xi
#print yi

# create figure, axes instances.
fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
# create Basemap instance for Robinson projection.
# coastlines not used, so resolution set to None to skip
# continent processing (this speeds things up a bit)
m = Basemap(projection='robin',lon_0=0.0,resolution=None)
# compute map projection coordinates of grid.
#print zi.mean(), zi.min(), zi.max()

x, y = m(xi, yi)
# draw line around map projection limb.
# color background of map projection region.
# missing values over land will show up this color.
m.drawmapboundary(fill_color='0.1')

mycmap=mpl.cm.jet
#mycmap.set_clim(0.,1.5)
cmax = zi.max()
#if cmax > 1.0: cmax = 1.0
mynorm = mpl.colors.Normalize(vmin=0.,vmax=cmax)
im1 = m.pcolor(x,y,zi,shading='flat',cmap=mycmap, norm=mynorm)
#im2 = m.pcolor(x,y,ice,shading='flat',cmap=plt.cm.gist_gray)
# draw parallels and meridians
m.drawparallels(np.arange(-90,90,30),labels=[1,1,0,1])
m.drawmeridians(np.arange(-180,180.,90.),labels=[1,1,0,1] )
# add colorbar
cb = m.colorbar(im1,"bottom", size="5%", pad="6%")
# add a title.
ax.set_title(' Map for file: %s'%(file))
plt.show()



