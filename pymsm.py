import os, sys

currdir = os.getcwd() + "/"
# i       = currdir.rfind("/UNILIB/")
# if ( i == -1 ) :
#   print("The example cannot be run from outside of the UNLIB directory")
#   print("or a subdirectory")
# else :
#   libdir = currdir[:i+10] + "lib"
#   sys.path.insert(0, libdir)
libdir = currdir + "UNILIB/lib"
sys.path.insert(0, libdir)
print("Library location is ", libdir)

# sys.path.insert(0,'./UNILIB/lib')

from PyUNILIB import PyUT990, PyUM510, PyUT540, PyUM520, PyUM522, PyUM524, PyUM536, PyUM530, PyUL220, PyUL225
from PyUNILIB import PyZDAT, PyZGEO, PyZVEC

from math import sqrt, acos, cos, pow, log10
import datetime
import numpy as np
 
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
            file = "./MAPS/"+year+"/AVKP"+kp+"T"+ut+".AVG"
            print(file)
            MapDB.maps[tag] = self.readMap(file)
        return MapDB.maps[tag]    
    
    def readMap(self, file):
        '''
        read in the map file
         
        '''
        try:
            lm, rc = np.loadtxt(file,skiprows=0,usecols = (3,5),unpack=True)
            lm = lm.reshape((37,73))  # maps are in 5x5 degs grid size and 37x73 points 
            rc = rc.reshape((37,73))             
#            return lm.transpose(),rc.transpose()
            rc = np.where(lm == 99.99, 0., rc) 
            return lm,rc,rc*lm**2
        except Exception as e:
            print(e, sys.stderr)
            
class PyMSM(object):
    def __init__(self, times, positions, kps, rc=None):
        '''
        main method to obtain the vertical rigidity for a given time and location
        
        Inputs:
         times: Numpy array of datetime - datetime.datetime()
         positions: Numpy array of geographic location [alti, lati, long] (km, deg)
         kps:Numpy array of kps
         rc: rigidity cut offs for which the transmission factor to be calculated. In a list or numpy 1D array 
         
         Note: the lengthes of the first 3 inputs should match
        
        '''
        self.Re = 6371.2  # Earth's radius in km as in UNILIB
        self.cyears = ['1955','1960','1965','1970','1975','1980','1985','1990','1995','2000','2005','2010','2015']
        self.cuts = ['00','03','06','09','12','15','18','21']
        self.ckps = ['0','1','2','3','4','5','6','7','8','9','X']
        try:
            if rc == None:
                # The rigidity values to be used for calculating the transmission function
                self.rc = [0.1, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, \
                           6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13, 14., 15., 16., 17., \
                           20., 25., 20., 30., 40., 50., 60.]
        except ValueError:
            self.rc = rc
        #
        #
        # Initialize the library; the output will be to "fort.1", but you can open a
        # suitable unit with number kunit.
        #
        kunit  = 1
        kinit  = 1
        ifail  = PyUT990(kunit, kinit)
        #  print ("Result from UT990 call = ", ifail)
        if ( ifail < 0 ) :
            print("Error in call to PyUT990.  ifail = ",ifail)
            return
        #
        #
        # Set the geomagnetic field model
        # (kint=0, DGRF/IGRF)
        #
        kint   = 0
        year   = 2015.0        # It can be reset by UM513 but it is pythonised yet. Perhaps call the again when the epock changes.
        lbint, ifail  = PyUM510(kint, year, kunit)
        #  print ("Result from UM510 call = ", lbint, ifail)
        if ( ifail < 0 ) :
            print("Error in call to PyUM510.  ifail = ",ifail)
            return
#    
        self.times = times
        self.coords = positions
        self.kps = kps
                
        #
        self.dbMgr = MapDB()
    
    def getTransmissionFunctions(self):
        '''
        Return: 
         TF: the  
        '''
        tf = []
        # first obtain the interpolated vertical cutoffs
        bm, lm, la, rcv = self.getRc()

#         # 2nd get the Invariant magnetic latitude, ether using the getInvLat() or calculateInvLat() method
#         invlats = self.getInvLat()
        # 3rd  get Transmission factors at the specificed rigidities
        mla = np.array(la)/57.2957795 # convert to radians
        tf = self.getTransfact(mla, rcv)
        # 4th get the Earth shadowing factors
        es = self.facshadow((self.coords[:,0]+self.Re)/self.Re)
        #
        return lm, bm, la, rcv, es, tf
        
    def getRc(self):
        '''
               
        ''' 

        rclist = []
        lmlist = []
        bmlist = []
        lalist = []
        #
        mkeyold = ''
        mdate = PyZDAT()
        mgeod = PyZGEO()
        kint = 0
        kunit = 1
        alpha = [90.,]
        nfbm = 1 
        # set the month and day for the map which are 01/01
        mdate.imonth  = 1
        mdate.iday    = 1       
        for i in range(len(self.times)):
            year = self.times[i].year
            if year < 1955: year = 1955
            if year > 2015: year = 2015
            iy = int((year - 1955)/5)
            ir = (year - 1955)%5
            cyear = self.cyears[iy]
            iu = int((self.times[i].hour+self.times[i].minute/60. + 1.5)/3.)
            # UT =1 corresponds to ut: 1.5 - 4.5 hrs
            if iu > 7: iu = 0
            cut = self.cuts[iu]
            ckp = self.ckps[self.kps[i]]
            self.dbMgr.getMap(cyear,ckp,cut)
            mkey = cyear+ckp+cut
            # set the year and ut of the map
            mdate.iyear   = int(cyear)
            mdate.ihour   = int(cut)
            #
            if mkey != mkeyold: 
                # initialise the internal field
                lbint, ifail  = PyUM510(kint, year, kunit)
                if ( ifail < 0 ) :
                    print("Error in call to PyUM510.  ifail = ",ifail)
                    return
                mkeyold = mkey
            # prepare the external field parameters, needed by GetBmLm()
            self.param = [self.kps[i], -30.0, 0.0, 25.0, 300.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            #
            # Set the geographic location
#           mgeod.radius  = self.coords[i][0]
            mgeod.colat   = 90 - self.coords[i][1] # note UNILIB requires the colatitude!
            mgeod.elong   = self.coords[i][2]

            mgeod.radius = 450.0 + self.Re  # at 450 km is the altitude of the map and coords are in GEO
        
            #
            dbm, dlm = self.GetBmLm(mdate, mgeod)
            lm450 = dlm[0]
            rc450 = self.getRC450km(mkey,self.coords[i][1],self.coords[i][2],lm450)
            #
            w =(ir + (self.times[i].timetuple().tm_yday -1 + self.times[i].hour/24.)/365.)/5.  # DoY starts from 1 
            # the next map at +5 years if required
            if 1955 < year < 2015 and w > 0.001:
                cyear = self.cyears [iy+1]
                self.dbMgr.getMap(cyear,ckp,cut)
                mkey = cyear+ckp+cut
                mdate.iyear   = int(cyear)
                if mkey != mkeyold: 
                    # initialise the interfield
                    lbint, ifail  = PyUM510(kint, int(cyear), kunit)
                    if ( ifail < 0 ) :
                        print("Error in call to PyUM510.  ifail = ",ifail)
                        return
                    mkeyold = mkey
                dbm, dlm = self.GetBmLm(mdate, mgeod)
                lm1 = dlm[0]                          
                rc1 = self.getRC450km(mkey,self.coords[i][1],self.coords[i][2],lm1)
                lm450 = lm450*(1-w) + w*lm1
                rc450 = rc450*(1-w) + w*rc1 
            #
            #calculate the (B,L) of the actual location at the given date and time
            # set the date and time
            mdate.iyear   = self.times[i].year
            mdate.imonth  = self.times[i].month
            mdate.iday    = self.times[i].day
            mdate.ihour   = self.times[i].hour
            mdate.imin    = self.times[i].minute
            mdate.isec    = self.times[i].second

            # initialise the interfield
            lbint, ifail  = PyUM510(kint, mdate.iyear, kunit)
            if ( ifail < 0 ) :
                print("Error in call to PyUM510.  ifail = ",ifail)
                return
            #
            # Set the geographic location
            mgeod.colat   = 90 - self.coords[i][1]
            mgeod.elong   = self.coords[i][2] # it could have been moved 
            mgeod.radius  = self.coords[i][0] + self.Re 
            dbm, dlm = self.GetBmLm(mdate, mgeod)
            # get the invariant latitude
            if dlm[0] == 99.99:
                frd=[1.0]
                fla=[75.]
            else:
                frd, fla, ifail = PyUL225(dbm, dlm, nfbm)
            if ( ifail < 0 ) :
                print("Error in call to PyUL225.  ifail = ",ifail)
                return
            bm = dbm[0]
            lm = dlm[0]
            rd = frd[0]
            la = fla[0] # in deg
            
# #     Write the result but to stdout, not kunit in this case.
# #
#             print("              magetic shell Lm   = %22.15E " %lm)
#             print("              field strength Bm  = %22.15E nT" %(bm*100000.))
#             print("              invariant radius   = %22.15E Re" %rd)
#             print("              invariant latitude = %22.15E deg" %la)
#             print(" -------------------------")
   
            #now apply altitude interpolation
            #
            rc = rc450*(lm450/lm)**2 # scaled by LM^2
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
            radist = mgeod.radius/self.Re  # in Re                    
            rcorr = log10(radist*radist)/14.   # Don used 11. but 14. is better
            rc -= rcorr 
            if rc < 0.: rc = 0.
            #
            rclist.append(rc)
            lmlist.append(lm)
            bmlist.append(bm*100000.)    # either B at the location or at the mirror point. Mirror point for UNILIB
            lalist.append(la)
        #           
        return bmlist, lmlist, lalist, rclist          
            
    def getRC450km(self,mkey,lat,lon,lm):
        '''
        interpolation to obtain Rc for the given location at 450km 
        it should be called after the mkey map has been prepared, e.g., after use of dbMgr.getMap()
        
        inputs:
            string mkey:  the map key which is cyear+ckp+cut
            float lat, lon: latitude and longitude in degrees
            float lm: the L shell number of the position at 450km altitude 
            
        outputs:
            float rc: evrtical rigidity cutoff in (GV) at 450km altitude
             
        '''
        # get the left-top box corner idxs
        i, j = self.getGridIdx(lat, lon)
        ii = (i+1 if i < 36 else 0)
        jj = (j+1 if j < 72 else 0)   

        # get the Rc*Lm^2 from the maps
        # left-top corner
        rclm_LT = self.dbMgr.maps[mkey][2][i,j]
        # right-top corner
        rclm_RT = self.dbMgr.maps[mkey][2][ii,j]
        # left-bot corner
        rclm_LB = self.dbMgr.maps[mkey][2][i,jj]
        # right-bot corner
        rclm_RB = self.dbMgr.maps[mkey][2][ii,jj]
        # 
        # get the weights
        wl,wr, wt, wb = self.getWeights(lat,lon)
        # 
        rclm_l = wt*rclm_LT + wb*rclm_LB 
        rclm_r = wt*rclm_RT + wb*rclm_RB
        #
        rclm = wl*rclm_l + wr*rclm_r
        # 
        return rclm/lm**2
    
    def GetBmLm(self, mdate, mgeod):
        kunit = 1
        # Compute the modified julian day
        # (based on January 1, 1950)
        #
        mdate1 = PyUT540(mdate)
        amjd  = mdate1.amjd
        #
        #
        # Set the external magnetic field model
        # (kext= 4, Tsyganenko 89c )
        #
        kext  = 4
        # the field parameters
        lbext, ifail  = PyUM520(kext, amjd, self.param, kunit)
        if ( ifail < 0 ) :
            print("Error in call to PyUM520.  ifail = ",ifail)
            return
        # position the Sun    
        PyUM522 (amjd, kunit)
        # convert SM coordinates
        PyUM524 (kunit)
        
        # Convert from geodetic to geocentric
        mpos = PyUM536 (mgeod)
        # get the L-shell value
        dbm, dlm, dkm, dsm, fbeq, ds, ifail = PyUL220(mpos, [90.,], 1)
        d = 0
        while ( ifail < 0 and d < 5 ) :
            #print("Error in call to PyUL220.  ifail = ",ifail, mpos.colat, mpos.elong)
            mgeod.elong += 1. # move the longitude east by 1 deg and recalculate
#             mgeod.colat += 1.
            d+=1
            mpos = PyUM536 (mgeod)
            dbm, dlm, dkm, dsm, fbeq, ds, ifail = PyUL220(mpos, [90.,], 1)
        if ifail < 0: dlm = [99.99,]
        return dbm, dlm
     
    def getWeights(self,lat, lon):
        
        # weights in longitude
        wr = (lon%5)/5. # left side of the box 
        wl = 1.0 - wr # right side of the box
        # weights in latitude
        wb = (-lat%5)/5. # bottom side of the box 
        wt = 1. - wb 
        return wl,wr, wt, wb
        
              
    def getGridIdx(self,lat, lon):
        ix = int(lon/5)
        iy = int(18 -lat/5)
        return iy, ix 

    def getInvLat(self):
        '''    
        calculate the invariant magnetic latitude of the given locations using irben
    
        Inputs:
            
        Returns:
            Emlats: list of the invariant magnetic latitudes in radians
        ''' 
        # calculate the corrected magnetic latitude
        # to to reset the altitudes = Re    
        radia = np.empty(len(self.times))
        radia.fill (self.Re)
        radia_old = self.coords[:,0] 
        self.coords[:,0] = radia
        #
        # GDZ -> MAG
        mpos = ib.coord_trans(self.coords,'MAG','sph')
        # restore the radi in coords
        self.coords.radi = radi_old      
        
        #
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
        return emlats    
            
#
    def calculateInvLat(self,B,L):
            '''
            calculte the invariant radial distance (R) and the magnetic latitude lambda based on the method of 
           
            "Roberts, C. S. (1964), Coordinates for the study of particles trapped in
            the Earth’s magnetic field: A method of converting from B, L to R, l
            coordinates, J. Geophys. Res., 69, 5089-- 5090."
            
            Inputs: 
                B, L:  B,L at given locations in numpy 1D arrays  
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
             mlat: the invariant magnetic latitude in radians, in numpy 1D array or list
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
        
        return 1. - facs
            
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
                 
        
def plotmap(xi,yi,zi):
    
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
    
def plotmap_b(xi,yi,zi):
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

def plotmap_c(xi,yi,zi,Title="Global Map"):    

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
    # map size
    nlon = 73
    nlat = 37
    res = 5  # grid size
# make up the grid
    xi = np.linspace(0, 365, nlon)    # X grid
    yi = np.linspace(90, -95, nlat)   # Y grid

    N = nlat*nlon 
    mapMgr = MapDB()
    kps = np.empty(N,dtype=np.int)
    kps.fill(1)    
    rc = [0.1, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, \
          6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13, 14., 15., 16., 17., \
          20., 25., 20., 30., 40., 50., 60.]  
#   
    dat = datetime.datetime(1962,7,1,12,0,0)  
    coords = []
    times = []
    alti = 500.
    for i in range(90,-95,-res): 
        for j in range(0,365,res):
            coords.append([alti,i+0.5,j+0.5])  # [alti, lati, longi]
            times.append(dat)
#  
    pm = PyMSM(np.array(times),np.array(coords),kps, np.array(rc))
    lm, bm, la, rcv, es, tf = pm.getTransmissionFunctions()
#
#    print (pm.getTransmissionFunctions())
    print ('completed')
    
#    lm, rc, rclm2 = mapMgr.getMap('2000','9','12')
    plotscatter(rc,np.sum(tf,axis=0)/N)
    plotmap_c(xi,yi,np.array(rcv).reshape((nlat,nlon)))

if __name__ == '__main__':
    main() 
