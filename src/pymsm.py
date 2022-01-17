import os, sys
from math import sqrt, acos, cos, pow, log10
import datetime
import numpy as np
from IRBEM import MagFields, Coords

# #from astropy.units import centiyear
# #from _pylief import NONE
# #from statsmodels.formula.api import wls
# #from astropy.wcs.docstrings import lat
# import spacepy.time as spt
# import spacepy.coordinates as spc
# import spacepy.irbempy as ib
# import spacepy.omni as om

# map size
nlon = 73
nlat = 37
ndata = nlon*nlat
# make up the grid
xi = np.linspace(0, 365, 73)    # X grid
yi = np.linspace(90, -95, 37)   # Y grid

    
class MapDB():
    # all instancs share the same map dict!
    maps = {}
    
    def __init__(self):
        '''
        Constructor
        '''
    
    def getMap(self, year, kp, ut):
        '''
        return the pre-calculated rigidity map for give year kp and ut
        
        inputs:
            year: string 
            kp: string
            ut: string
        '''
        tag = year+kp+ut
        if tag not in MapDB.maps:
            file = "../MAPS/"+year+"/AVKP"+kp+"T"+ut+".AVG"
            print(file)
            MapDB.maps[tag] = self.readMap(file)
        return MapDB.maps[tag]    
    
    def readMap(self, file):
        '''
        read in the map file
         
        '''
        try:
            lm, rc = np.loadtxt(file,skiprows=0,usecols = (3,5),unpack=True)
            lm = lm.reshape((nlat,nlon))
            rc = rc.reshape((nlat,nlon))             
#            return lm.transpose(),rc.transpose()
            return lm,rc,rc*lm**2
        except Exception as e:
            print(e, sys.stderr)
            
class PyMSM(object):
    def __init__(self, times, positions, kps, rc=None):
        '''
        
        Inputs:
         times: 1D numpy arrary of date and time in datetime.datetime format
         positions: 2D numpy array locations corresponding to the times in GDZ coordinates, e.g., [[alt0, lat0, lon0],[alt1, lat1, lon1],...]
         kps: 1D numpy array of Kp indices, corresponding to the times  
         rc: 1D numpy array of rigidity cut offs, in GV, for which the transmission factor to be calculated. This is optional, if not specified, the defaults are:
                         [0.1, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, \
                          6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13, 14., 15., 16., 17., \
                        20., 25., 20., 30., 40., 50., 60.]
         
         Note: the length of the first 3 inputs should match!
        
        '''
        self.times = times
        self.model = MagFields(options = [0,30,0,0,0], kext=4, verbose = False)
        self.cyears = ['1955','1960','1965','1970','1975','1980','1985','1990','1995','2000','2005','2010','2015','2020','2025']
        self.cuts = ['00','03','06','09','12','15','18','21']
        self.ckps = ['0','1','2','3','4','5','6','7','8','9','X']
        self.kps = kps
        try:
            if rc == None:
                # The rigidity values to be used for calculating the transmission function
                self.rc = [0.1, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, \
                           6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13, 14., 15., 16., 17., \
                           20., 25., 20., 30., 40., 50., 60.]
        except ValueError:
            self.rc = rc
        #  positions are in  GDZ 
        self.positions = Coords().coords_transform(times, positions, 'GDZ', 'GDZ')
        self.radius = Coords().coords_transform(times, self.positions, 'GDZ', 'RLL')[:,0]
        #
        self.lla = {}
        self.lla['x1'] = self.positions[:,0]
        self.lla['x2'] = self.positions[:,1]
        self.lla['x3'] = self.positions[:,2]
        self.lla['dateTime'] = times
        self.maginput = {'Kp':kps*10.}
        self.model.make_lstar(self.lla, self.maginput)
        # the (B,L) for the inputs times and positions
        self.lm = np.abs(self.model.make_lstar_output['Lm'])
        self.bm = np.abs(self.model.make_lstar_output['blocal'])        
        #
        self.dbMgr = MapDB()
    
    def getTransmissionFunctions(self):
        '''
        Main method to obtain the vertical rigidity for the given times and locations ...
        
        Return: 
         LMs: the   
        '''
        # first obtain the interpolated vertical cutoffs
        Rcv = self.getRc()
        # 2nd get the magnetic latitude, ether using the getEMLat() or calculateRInv method
        Mlats = self.getEMLat()
        # 3rd 
        TF = self.getTransfact(Mlats, Rcv)
        # 4th get the Earth shadowing factors       
        ES = self.facshadow(self.radius)  
        #
        return self.lm, self.bm, Mlats, Rcv, ES, TF
        
    def getRc(self):
        '''
               
        ''' 

        t_utc = self.lla['dateTime']

        rclist = []
        local = {}
        #
                     
        for i in range(len(t_utc)):
            year = t_utc[i].year
            if year < 1955: year = 1955
            if year > 2025: year = 2025
            iy = int((year - 1955)/5)
            ir = (year - 1955)%5
            cyear = self.cyears [iy]
            iu = int((t_utc[i].hour+t_utc[i].minute/60. + 1.5)/3.)
            # UT =1 corresponds to ut: 1.5 - 4.5 hrs
            if iu > 7: iu = 0
            cut = self.cuts[iu]
            ckp = self.ckps[self.kps[i]]
            self.dbMgr.getMap(cyear,ckp,cut)
            mkey = cyear+ckp+cut
            #  in first map 
            local['x1'] = 450.
            local['x2'] = self.lla['x2'][i]
            local['x3'] = self.lla['x3'][i]
            local['dateTime'] = datetime.datetime(int(cyear),1,1,int(cut)).isoformat()
            maginput = {'Kp':self.maginput['Kp'][i]}
            self.model.make_lstar(local, maginput)
            lm = np.abs(self.model.make_lstar_output['Lm'])
            rc = self.getRC450km(mkey, self.lla['x2'][i], self.lla['x3'][i],lm)
            #
            w =(ir + (t_utc[i].timetuple().tm_yday + (t_utc[i].hour + t_utc[i].minute/60.)/24.)/365.)/5. 
            # the next map at +5 years if required
            #
            if 1955 < year < 2025 and w < 1.:
                cyear = self.cyears [iy+1]
                self.dbMgr.getMap(cyear,ckp,cut)
                mkey = cyear+ckp+cut
                # in 2nd map  
                local['dateTime'] =  datetime.datetime(int(cyear),1,1,int(cut)).isoformat()  
                self.model.make_lstar(local, maginput)
                lm1 = np.abs(self.model.make_lstar_output['Lm'])
                rc1 = self.getRC450km(mkey,self.lla['x2'][i], self.lla['x3'][i],lm1)
                lm = lm*w + (1.- w)*lm1
                rc = rc*w + (1.- w)*rc1 
            #    
            #now apply altitude interpolation
            #
            lmr = self.lm[i] # for real time and altitude
            #
            rcr = rc*(lm/lmr)**2 # scaled by LM^2
            '''        
            Further Radial Distance Adjustment according to email from Don on 09/12/2013
             " If you examine the cutoff interpolation FORTRAN code in detail, 
               you may notice that there is an adjustment in the "L" value altitude
               interpolation process at the end of Subroutine LINT5X5 in a section
               labeled "Radial Distance Adjustment".  While the "L" interpolation
               equation has the basic form of L**-2, when this exact form is used to
               extend to geosynchronous altitude, the cutoff values extrapolated from
               the near Earth low altitudes are too high; approximately 0.3 GV at the
               magnetic equator.  (See figure 7 of Shea & Smart, JGR, 72, 3447, 1967)
               We incorporated a "patch" (actually an "ad-hoc" exponential function
               that makes adjustments so at 6.6 earth radii the vertical cutoff rigidity
               for local noon at the magnetic equator under extremely quiet magnetic
               conditions is about zero (or extremely small).
               This "ad-hoc" exponential function is not going to be reliable
               beyond geosynchronous distances."
            '''
            radist = self.radius[i]+1.                    
            rcorr = log10(radist*radist)/14.   # Don used 11. but 14. is better
            rcr -= rcorr 
            if rcr < 0.: rcr = 0.
            try: 
                rct = rcr[0]
            except:
                rct = rcr
            #
            rclist.append(rct)
            
        return np.array(rclist)          
            
    def getRC450km(self,mkey,lon,lat,lm):
        '''
        interpolation to obtain Rc for the given location at 450km 
        it should be called after the mkey map has been prepared, e.g., after use of dbMgr.getMap()
        
        inputs:
            string mkey:  the map key which is cyear+ckp+cut
            float lon, lat: longitude and latitude in degrees
            float lm: the L shell number of the position at 450km altitude 
            
        outputs:
            float rc: evrtical rigidity cutoff in (GV) at 450km altitude
             
        '''
        # get the left-top box corner idxs
        i, j = self.getGridIdx(lon, lat)
        #get the weights
        wl,wr,wt,wb = self.getWeights(lon,lat)
        # get the Lm and Rc from the maps
        # left-top corner
        rclm_LT = self.dbMgr.maps[mkey][2][i,j]
        # right-top corner
        rclm_RT = self.dbMgr.maps[mkey][2][i+1,j]
        # left-bot corner
        rclm_LB = self.dbMgr.maps[mkey][2][i,j+1]
        # right-bot corner
        rclm_RB = self.dbMgr.maps[mkey][2][i+1,j+1]
        # 
        # get the weights
        wl,wr, wt, wb = self.getWeights(lon,lat)
        # 
        rclm_l = wt*rclm_LT + wb*rclm_LB 
        rclm_r = wt*rclm_RT + wb*rclm_RB
        #
        rclm = wl*rclm_l + wr*rclm_r
        # 
        return rclm/lm**2
     
    def getWeights(self,lon,lat):
        
        # weights in longitude
        wl = (lon%5)/5. # left side of the box 
        wr = 1.0 - wl # right side of the box
        # weights in latitude
        wb = (-lat%5)/5. # bottom side of the box 
        wt = 1. - wb 
        return wl,wr, wt, wb
        
              
    def getGridIdx(self,lon, lat):
        ix = int(lon/5)
        iy = int(18 -lat/5)
        return ix, iy 

    def getEMLat(self):
        '''    
        calculate the equivalent magnetic latitude of the given locations 
    
            Get Corrected geomagnetic latitude (GMLATC) at sub-satellite point
            Then, get Invariant latitude at satellite position
                       INVARIANT LAT = ACOS(1.0/SQRT(L))
            Select the smaller value for magnetic latitud
                   We will always use absolute value of equivalent magnetic latitude
        Inputs:
            
        Returns:
            Emlats: np.array of the equivalent magnetic latitudes in radians
        ''' 
        # calculate the corrected magnetic latitude
        # need reset the altitudes = Re    

        #
        # GEO -> MAG ->sph
        mpos = Coords().coords_transform(self.times, self.positions, 'GEO', 'MAG')
        mpos = Coords().coords_transform(self.times, self.positions, 'MAG', 'RLL')
        
        # mpos[:,1] are the magnetic latitude in degrees. Note this is not the same as 
        # the corrected geomagnetic latitude, but the difference should be small
        gmlatcr = mpos[:,1]/57.2957795 # convert to radians
        emlats = []
        for i in range(len(self.lm)):
            glmdar = 0.0
            if (self.lm[i] > 1.): glmdar = acos(1.0/sqrt(self.lm[i]))
            if abs(gmlatcr[i])< abs(glmdar):  glmdar = abs(gmlatcr[i])
            emlats.append(glmdar)
        #
        return np.array(emlats)    
    

    def calculateRInv(self,B,L):
            '''
            calculte the invariant radial distance (R) and the magnetic latitude lambda based on the method of 
           
            "Roberts, C. S. (1964), Coordinates for the study of particles trapped in
            the Earth’s magnetic field: A method of converting from B, L to R, l
            coordinates, J. Geophys. Res., 69, 5089-- 5090."
            
            Not used at the moment, Need to compare lambda vs emlats above!
            
            Inputs: 
                B, L: in numpy 1D arrays  
            Ouputs:
                R, lambda: in numpy 1D arrays 
                
            '''
            a = [1.25992106, -0.19842592, -0.04686632, -0.01314096, -0.00308824, 0.00082777, -0.00105877, 0.00183142]
            Md = 31165.3  #nT*Re^3
            if (L<0. ): return -1., -1.
            p = np.pow(np.pow(L,3.)*B/Md,-1./3.)
            #
            s = 0.
            for i in range(8):
                s += a[i]*np.pow(p,i)
            ps = p*s 
            for i in np.nonzero(ps>1.):
                ps[i] = 1.
            #
            lamb = np.degrees(np.acos(np.sqrt(ps)))
            R = L*ps
            for i in np.nonzer(p<0. or p > 1.):
                R[i] = -1.
                lamb[i] = -1.
            return R, lamb
    
    def getTransfact(self, mlat, rcv):
        '''

         computes the angle-averaged transmmison function at the given RCs    
         averaging over arrival directions. 

         Inputs:
             mlat: the magnetic latitude in radians, in numpy 1D array or list
             rcv: the vertical cut-off, in numpy 1D array or list
        
         Returns:
             facs: the transmission function at the specified rigidities. In numpy 1D array or list   
         
        '''

        #
        N = len(self.rc)
        fac = np.empty(N)
        facs = np.empty(shape=(len(mlat),N))
        rcv = np.array(rcv)
        # 
        fac.fill(1.0)
        for i in range(len(rcv)):
            if (rcv[i] <= 0.1):
                facs[i] = fac   
            else:
                facs[i] = self.calcTF(mlat[i],rcv[i])
        
        return facs
            
    def calcTF(self, mlat, rcv):
        '''
        '''
        # table of One-Minus-Cos-Angles-Over-2 :
        omcao2 = np.array([0., .067, .146, .25, .5, .75, .854, .933, 1.])
        ang = np.array([0.01, .5236, .785, 1.047, 1.571, 2.094, 2.356, 2.618, 3.1416])
        nangle = 9
        #
        cosa = np.cos(ang)
        cosl = cos(mlat)
        cut = 4.*rcv/(1.0+np.sqrt(1.0+cosa*cosl**3))**2
        #
        N = len(self.rc)
        fac = np.empty(N)        
        for ir in range(N):
            fac[ir] = 0.  
            if (self.rc[ir] >= cut[0]):
                fac[ir] = 1.0  
                for ia in range(1,nangle):
                    # find the angular location where the cutoff goes over the rc:
                    if (self.rc[ir] <= cut[ia]):
                        fac[ir] = omcao2[ia-1] + (self.rc[ir]-cut[ia-1])*(omcao2[ia]-omcao2[ia-1]) \
                        /(cut[ia]-cut[ia-1])
                        break
        return fac
                        
    def facshadow(self, R):
        '''
    
       This is a correction factor for the earth's shadow on 
       the spacecraft according to simple geometrical optics.
       
       Inputs:
           R = radius in Re, in a list or 1d numpy array
       Return:
       
        '''
        fac = 1. - .5 * (1.-np.sqrt(R**2-1.)/R)
    
        return fac
                 
        
def plotmap(zi):
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    
    plt.subplot(aspect=1, title='Global_map', ylim=[-90, 90], xlim=[0, 360])
    pc = plt.pcolor(xi, yi, zi, norm=LogNorm(vmin=1e-4, vmax=1e2))
    cs = plt.contour(xi, yi, zi, np.logspace(-5, 2, 8), colors='pink', linewidth=0.5, linesytle='dashed')
    plt.clabel(cs, inline=1)

    plt.colorbar(pc)
    plt.show()

    #plt.savefig('plot_density_meridian.png')#
    
def plotmap_b(zi):
    import sys 
    from mpl_toolkits.basemap import Basemap
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    #from matplotlib.mlab import griddata

    # create figure, axes instances.
    fig = plt.figure()
    ax = fig.add_axes([0.05,0.05,0.9,0.9])
    # create Basemap instance for Robinson projection.
    # coastlines not used, so resolution set to None to skip
    # continent processing (this speeds things up a bit)
    m = Basemap(projection='robin',lon_0=0.0,resolution=None)
    # compute map projection coordinates of grid.
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
    m.drawmeridians(np.arange(-180,180.,90.),labels=[1,1,0,1] )
    # add colorbar
    cb = m.colorbar(im1,"bottom", size="5%", pad="6%")
    # add a title.
    ax.set_title(' Map for file: %s'%("Global Map"))
    plt.show()

def plotmap_c(zi,Title="Global Map"):    

    import matplotlib.pyplot as plt
    
    # contour the gridded data, plotting dots at the nonuniform data points.
    plt.figure(figsize=(11,7))
    plt.subplots_adjust(right=1.0)
    CS = plt.contour(xi,yi,zi,20,linewidths=0.5,colors='k')
    plt.clabel(CS, inline=1, fontsize=10)
    CS = plt.contourf(xi,yi,zi,20,cmap=plt.cm.rainbow,
                  vmax=abs(zi).max(), vmin=-abs(zi).max())
    plt.colorbar() # draw colorbar
    plt.xlabel("Longitude [Deg]")
    plt.ylabel("Latitude [Deg]")
    plt.title(Title)
    plt.show()    
    #plt.savefig('rigidity_map.png')#
    
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

 
def main():
    '''

    '''
    N = 19*37
    res = 10  # grid size
#
    mapMgr = MapDB()
#
    kps = np.empty(N,dtype=int)
    kps.fill(1)
    # rc = [0.1, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, \
    #   6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13, 14., 15., 16., 17., \
    #   20., 25., 20., 30., 40., 50., 60.] 
#
    times = np.empty(N,dtype='object')   
    #dat = datetime.datetime(2012,7,1,12,0,0)
    dat = datetime.datetime.fromisoformat('2012-07-01T12:00:00')
    
    times.fill(dat)
#      
    coord = []
    alti = 500.
    for i in range(90,-100,-res): 
        for j in range(0,370,res):
            coord.append([alti,i,j])  # [alti, lati, longi]

#  
    pm = PyMSM(times,np.array(coord),kps)


    print (pm.getTransmissionFunctions())
    print ('completed')
    
#    lm, rc, rclm2 = mapMgr.getMap('2000','9','12')
#    plotscatter(lm,rc)
#    plotmap_c(rc)

if __name__ == '__main__':
    main() 