import sys 

if len(sys.argv) is not 8 :  
    print "To run: difstat.py file1 skip-lines z-ol file2 skip-lines z-col rc"
    sys.exit(1)

file1 = sys.argv[1]
sl1 = int(sys.argv[2])
cz1 = int(sys.argv[3])
file2 = sys.argv[4]
sl2 = int(sys.argv[5])
cz2 = int(sys.argv[6])

rc = float(sys.argv[7])

import numpy as np

z = np.loadtxt(file1,skiprows=sl1,usecols = (cz1,),unpack=True)
z1 = np.loadtxt(file2,skiprows=sl2,usecols = (cz2,),unpack=True)
z0 = np.where(z1>rc,(z - z1)/z1,0.)
x = z0[z0!=0.]

mu = np.mean(x)
sigma = np.std(x)

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

num_bins = 50
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
# add a 'best fit' line
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.xlabel('Errors in Rc')
plt.ylabel('Probability')
plt.title(r'$\mu=$'+str(mu)+r' $\sigma=$'+str(sigma))

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()

