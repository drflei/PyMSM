"""Plotting utilities for PyMSM visualization."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata


def plotmap(xi, yi, zi):
    """Plot a 2D map with logarithmic color scale.
    
    Args:
        xi: X coordinates (longitude).
        yi: Y coordinates (latitude).
        zi: Z values to plot.
    """
    plt.figure(figsize=(11, 7))
    plt.subplots_adjust(right=1.0)

    plt.subplot(aspect=1, title='Global_map', ylim=[-90, 90], xlim=[0, 360])
    pc = plt.pcolor(xi, yi, zi, norm=LogNorm(vmin=1e-1, vmax=20))
    cs = plt.contour(xi, yi, zi, np.logspace(-5, 2, 8), colors='orange')
    plt.clabel(cs, inline=1)

    plt.colorbar(pc)
    plt.show()


def plotmap_basemap(xi, yi, zi):
    """Plot a map using Basemap projection.
    
    Args:
        xi: X coordinates (longitude).
        yi: Y coordinates (latitude).
        zi: Z values to plot.
    """
    import matplotlib as mpl
    from mpl_toolkits.basemap import Basemap

    # Create figure, axes instances
    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])

    # Create Basemap instance for Robinson projection
    # Coastlines not used, so resolution set to None to skip continent processing
    m = Basemap(projection='robin', lon_0=-180, resolution=None)
    # Compute map projection coordinates of grid
    x, y = m(xi, yi)

    # Draw line around map projection limb
    # Color background of map projection region
    # Missing values over land will show up this color
    m.drawmapboundary(fill_color='0.1')

    mycmap = mpl.cm.jet
    cmax = zi.max()
    mynorm = mpl.colors.Normalize(vmin=0., vmax=cmax)
    im1 = m.pcolor(x, y, zi, cmap=mycmap, norm=mynorm)
    # Draw parallels and meridians
    m.drawparallels(np.arange(-90, 90, 30), labels=[1, 1, 0, 1])
    m.drawmeridians(np.arange(0, 360., 60.), labels=[1, 1, 0, 1])
    # Add colorbar
    cb = m.colorbar(im1, "bottom", size="5%", pad="6%")
    # Add a title
    ax.set_title('Map for file: %s' % ("Global Map"))
    plt.show()


def plotmap_contour(xi, yi, zi, Title="Global Map"):
    """Plot a contour map with gridded data.
    
    Args:
        xi: X coordinates (longitude).
        yi: Y coordinates (latitude).
        zi: Z values to plot.
        Title: Plot title (default: "Global Map").
    """
    # Grid the data
    grid_x, grid_y = np.mgrid[0:360:360j, -90:90:180j]
    grid_z = griddata((xi.flatten(), yi.flatten()), zi.flatten(), (grid_x, grid_y), method='cubic')

    # Contour the gridded data
    plt.figure(figsize=(11, 7))
    plt.subplots_adjust(right=1.0)
    CS = plt.contour(grid_x, grid_y, grid_z, 30, linewidths=0.5, colors='k')
    plt.clabel(CS, inline=1, fontsize=10)
    CS = plt.contourf(grid_x, grid_y, grid_z, 30, cmap=plt.cm.rainbow,
                      vmax=abs(grid_z).max(), vmin=-abs(grid_z).max())
    plt.colorbar()  # Draw colorbar
    plt.xlabel("Longitude [Deg]")
    plt.ylabel("Latitude [Deg]")
    plt.title(Title)
    plt.show()


def plotmapfile(file):
    """Plot a map from a data file.
    
    Args:
        file: Path to the data file.
    """
    x, y, z = np.loadtxt(file, skiprows=0, usecols=(1, 0, 5), unpack=True)

    # Define grid
    xg, yg = np.mgrid[0:360:360j, -90:90:180j]
    # Grid the data
    zg = griddata((x, y), z, (xg, yg), method='cubic')

    # Contour the gridded data
    plt.figure(figsize=(11, 7))
    plt.subplots_adjust(right=1.0)
    CS = plt.contour(xg, yg, zg, 30, linewidths=0.5, colors='k')
    plt.clabel(CS, inline=1, fontsize=10)
    CS = plt.contourf(xg, yg, zg, 30, cmap=plt.cm.rainbow,
                      vmax=abs(zg).max(), vmin=-abs(zg).max())
    plt.colorbar()  # Draw colorbar

    plt.xlabel("Longitude [Deg]")
    plt.ylabel("Latitude [Deg]")
    plt.title('File: ' + file)
    plt.show()


def plotscatter(x, y, xtit='x-axis', ytit='y-axis', title='x-y scatter plot'):
    """Create a scatter plot with logarithmic x-axis.
    
    Args:
        x: X-axis data.
        y: Y-axis data.
        xtit: X-axis title (default: 'x-axis').
        ytit: Y-axis title (default: 'y-axis').
        title: Plot title (default: 'x-y scatter plot').
    """
    plt.style.use('seaborn-whitegrid')
    plt.figure(figsize=(11, 7))
    plt.plot(x, y, '.', color='black')
    plt.xlabel(xtit)
    plt.ylabel(ytit)
    plt.title(title)
    plt.xscale('log')
    plt.show()


def plot3D(hists, xgrid=None, ygrid=None, xtit='X', ytit='Y', ztit='Z', title='3D plot'):
    """Create a 3D bar plot.
    
    Args:
        hists: 2D array of histogram data.
        xgrid: Optional x-axis grid.
        ygrid: Optional y-axis grid.
        xtit: X-axis title (default: 'X').
        ytit: Y-axis title (default: 'Y').
        ztit: Z-axis title (default: 'Z').
        title: Plot title (default: '3D plot').
    """
    from mpl_toolkits.mplot3d import axes3d

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111, projection='3d')
    dm = hists.shape
    if xgrid is None:
        hist = np.linspace(0, dm[1], dm[1])
    else:
        hist = xgrid
    if ygrid is None:
        nbins = dm[0]
        z_array = np.linspace(0, nbins * 10, nbins)
    else:
        z_array = ygrid

    for a, z in zip(hists, z_array):
        ax.bar(hist, a, zs=z, zdir='y', alpha=0.8)

    ax.set_xlabel(xtit)
    ax.set_ylabel(ytit)
    ax.set_zlabel(ztit)
    ax.set_title(title)

    plt.show()


def plotbar3d(hists, colname=None, rowname=None, xtit='X', ytit='Y', ztit='Z', title='Bar3D plot'):
    """Create a 3D bar plot with optional axis labels.
    
    Args:
        hists: 2D array of histogram data.
        colname: Optional column names for x-axis.
        rowname: Optional row names for y-axis.
        xtit: X-axis title (default: 'X').
        ytit: Y-axis title (default: 'Y').
        ztit: Z-axis title (default: 'Z').
        title: Plot title (default: 'Bar3D plot').
    """
    from mpl_toolkits.mplot3d import axes3d

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111, projection='3d')

    dm = hists.shape

    lx = dm[1]  # Work out matrix dimensions
    ly = dm[0]
    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)

    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x = 0.5 * np.ones_like(zpos)
    y = x.copy()
    z = hists.flatten()

    ax.bar3d(xpos, ypos, zpos, x, y, z)

    if colname is not None:
        ax.w_xaxis.set_ticklabels(colname)
    if rowname is not None:
        ax.w_yaxis.set_ticklabels(rowname)
    ax.set_xlabel(xtit)
    ax.set_ylabel(ytit)
    ax.set_zlabel(ztit)
    ax.set_title(title)

    plt.show()
    