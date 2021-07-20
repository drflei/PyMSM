import sys 

if len(sys.argv) is not 7 :  
    print "To run: plotcontour.py file1 file2 skip-lines x-col y-col z-ol"
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]
sl = int(sys.argv[3])
cx = int(sys.argv[4])
cy = int(sys.argv[5])
cz = int(sys.argv[6])

#from mpl_toolkits.basemap import Basemap
from matplotlib.mlab import griddata

#import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np

x,y,z = np.loadtxt(file1,skiprows=sl,usecols = (cx,cy,cz),unpack=True)
z1 = np.loadtxt(file2,skiprows=sl,usecols = (cz,),unpack=True)
z = z - z1

#with open(file, "r") as f:
#    for line in f:
#        entries = line.split()
#        x.append(float(entries[1]))
#        y.append(float(entries[0]))
#        z.append(float(entries[4]))
       
# define grid.
xi = np.linspace(-0.1,360.1,360)
yi = np.linspace(-90.1,90.1,180)
# grid the data.
zi = griddata(x,y,z,xi,yi)

# contour the gridded data, plotting dots at the nonuniform data points.
plt.figure(figsize=(11,7))
plt.subplots_adjust(right=1.0)
CS = plt.contour(xi,yi,zi,20,linewidths=0.5,colors='k')
plt.clabel(CS, inline=1, fontsize=10)
CS = plt.contourf(xi,yi,zi,20,cmap=plt.cm.rainbow,
                  vmax=abs(zi).max(), vmin=-abs(zi).max())
plt.colorbar() # draw colorbar
# plot data points.

#plt.scatter(x,y,marker='o',c='b',s=5,zorder=10)

#plt.xlim(-2,2)
#plt.ylim(-2,2)


plt.xlabel("Longitude [Deg]")
plt.ylabel("Latitude [Deg]")
plt.title(' File: ' + file1 + '-' + file2)
plt.show()
