def plotmap(xi, yi, zi):
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    
    plt.figure(figsize=(11,7))
    plt.subplots_adjust(right=1.0)
        
    plt.subplot(aspect=1, title='Global_map', ylim=[-90, 90], xlim=[0, 360])
    pc = plt.pcolor(xi, yi, zi, norm=LogNorm(vmin=1e-1, vmax=20))
    cs = plt.contour(xi, yi, zi, np.logspace(-5, 2, 8), colors='orange')
    plt.clabel(cs, inline=1)

    plt.colorbar(pc)
    plt.show()

    #plt.savefig('plot_density_meridian.png')#
    
def plotmap_basemap(xi, yi, zi):

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import numpy as np

    # create figure, axes instances.
    fig = plt.figure(figsize=(11,7))
    ax = fig.add_axes([0.05,0.05,0.9,0.9])
        
    # create Basemap instance for Robinson projection.
    # coastlines not used, so resolution set to None to skip
    # continent processing
    m = Basemap(projection='robin',lon_0=-180,resolution=None)
    # compute map projection coordinates of grid.
    #
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
    im1 = m.pcolor(x,y,zi,cmap=mycmap, norm=mynorm)
    #im2 = m.pcolor(x,y,ice,shading='flat',cmap=plt.cm.gist_gray)
    # draw parallels and meridians
    m.drawparallels(np.arange(-90,90,30),labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360.,60.),labels=[1,1,0,1] )
    # add colorbar
    cb = m.colorbar(im1,"bottom", size="5%", pad="6%")
    # add a title.
    ax.set_title(' Map for file: %s'%("Global Map"))
    plt.show()

def plotmap_contour(xi, yi, zi,Title="Global Map"):    

    import matplotlib.pyplot as plt   
    from scipy.interpolate import griddata
    import numpy as np

#    grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
#    grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
#    grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
    # grid the data.
    grid_x, grid_y = np.mgrid[0:360:360j, -90:90:180j]
    grid_z = griddata((xi.flatten(),yi.flatten()),zi.flatten(),(grid_x,grid_y), method='cubic')

    # contour the gridded data, plotting dots at the nonuniform data points.
    plt.figure(figsize=(11,7))
    plt.subplots_adjust(right=1.0)
    CS = plt.contour(grid_x,grid_y,grid_z,30,linewidths=0.5,colors='k')
    plt.clabel(CS, inline=1, fontsize=10)
    CS = plt.contourf(grid_x,grid_y,grid_z,30,cmap=plt.cm.rainbow,
                  vmax=abs(grid_z).max(), vmin=-abs(grid_z).max())
    plt.colorbar() # draw colorbar
    plt.xlabel("Longitude [Deg]")
    plt.ylabel("Latitude [Deg]")
    plt.title(Title)
    plt.show()    
    #plt.savefig('rigidity_map.png')#
    
def plotmapfile(file):
    
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy as np
    
    x,y,z = np.loadtxt(file,skiprows=0,usecols = (1,0,5),unpack=True)
    
    # define grid.
#    xi = np.linspace(-0.1,360.1,360)
#    yi = np.linspace(-90.1,90.1,180)
    xg, yg = np.mgrid[0:360:360j, -90:90:180j]
    # grid the data.
    zg = griddata((x,y),z,(xg,yg),method='cubic')
    
    # contour the gridded data, plotting dots at the nonuniform data points.
    plt.figure(figsize=(11,7))
    plt.subplots_adjust(right=1.0)
    CS = plt.contour(xg,yg,zg,30,linewidths=0.5,colors='k')
    plt.clabel(CS, inline=1, fontsize=10)
    CS = plt.contourf(xg,yg,zg,30,cmap=plt.cm.rainbow,
                      vmax=abs(zg).max(), vmin=-abs(zg).max())
    plt.colorbar() # draw colorbar
    
    plt.xlabel("Longitude [Deg]")
    plt.ylabel("Latitude [Deg]")
    plt.title(' File: ' + file)
    plt.show()
    
def plotscatter(x,y,xtit='x-axis',ytit='y-axis',title='x-y scatter plot'):
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-whitegrid')
    plt.figure(figsize=(11,7))
    plt.plot(x, y, '.', color='black')
    plt.xlabel(xtit)
    plt.ylabel(ytit)
    plt.title(title)
    plt.xscale('log')
    plt.show()

def plot3D(hists, xgrid=None, ygrid=None, xtit='X', ytit='Y',ztit='Z', title='3D plot' ):
    
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig = plt.figure(figsize=(11,7))
    ax = fig.add_subplot(111, projection='3d')
    dm = hists.shape
    if xgrid is None: 
        hist = np.linspace(0,dm[1],dm[1])
    else:
        hist = xgrid
    if ygrid is None:
        nbins = dm[0]
        z_array = np.linspace(0,nbins*10,nbins)
    else:
        z_array = ygrid
        
    for a,z in zip(hists, z_array):
        ax.bar(hist, a, zs=z, zdir='y', alpha=0.8)
    
    ax.set_xlabel(xtit)
    ax.set_ylabel(ytit)
    ax.set_zlabel(ztit)
    ax.set_title(title)
    
    plt.show()

def plotbar3d(hists, colname = None, rowname = None, xtit='X', ytit='Y',ztit='Z', title='Bar3D plot'):
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig = plt.figure(figsize=(11,7))
    ax = fig.add_subplot(111, projection='3d')

    dm =hists.shape
    
    lx= dm[1]            # Work out matrix dimensions
    ly= dm[0]
    xpos = np.arange(0,lx,1)    # Set up a mesh of positions
    ypos = np.arange(0,ly,1)
    xpos, ypos = np.meshgrid(xpos+0.25, ypos+0.25)

    xpos = xpos.flatten()   # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)

    x = 0.5 * np.ones_like(zpos)
    y = x.copy()
    z = hists.flatten()

    ax.bar3d(xpos, ypos, zpos, x, y, z)
    
    if colname is not None: ax.w_xaxis.set_ticklabels(colname)
    if rowname is not None: ax.w_yaxis.set_ticklabels(rowname)
    ax.set_xlabel(xtit)
    ax.set_ylabel(ytit)
    ax.set_zlabel(ztit)
    ax.set_title(title)
    
    plt.show()
    