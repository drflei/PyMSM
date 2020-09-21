'''
This version included the handling of spenvis csv files, but it looks like the spenvis module copied from CIRSOS
is not working with the oribital file ...
Need to test the spenvis module and add csv file out put functions
 
'''

import os, sys
from math import sqrt, acos, cos, pow
import numpy as np
#from astropy.units import centiyear
#from _pylief import NONE
#from statsmodels.formula.api import wls
#from astropy.wcs.docstrings import lat
import spacepy.time as spt
import spacepy.coordinates as spc
import spacepy.irbempy as ib

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
            lm = lm.reshape((nlat,nlon))
            rc = rc.reshape((nlat,nlon))             
#            return lm.transpose(),rc.transpose()
            return lm,rc,rc*lm**2
        except Exception as e:
            print(e, sys.stderr)
            
class Interpolations(object):
    def __init__(self, times, positions, kps):
        '''
        main method to obtain the vertical rigidity for a given time and location
        
        Inputs:
         times: a SpacyPy Ticktock instance
         positions: a SpacePy Coords instance
         kps: a list of kps
         
         Note: the lengthes of the inputs should match
        
        Outputs:
        '''
        self.cyears = ['1955','1960','1965','1970','1975','1980','1985','1990','1995','2000','2005','2010','2015']
        self.cuts = ['00','03','06','09','12','15','18','21']
        self.ckps = ['0','1','2','3','4','5','6','7','8','9','X']
        self.times = times
        positions.ticks = self.times  # need to set the ticks before convert
        self.coords = positions.convert('GDZ','sph')
        # get the Lms at actual positions, they wil be used for altitude scaling 
        t_dic = ib.get_Lm(self.times,self.coords,90,'T89') # alternatively alpha=[90]   
        self.lm = t_dic['lm'].flatten()
        t_dic = ib.get_Bfield(self.times,self.coords,extMag='T89')
        self.bm = t_dic['Blocal'].flatten()
        self.kps = kps
        
    def getRc(self):
        '''
               
        ''' 

        t_utc = self.times.UTC

        rclist = []
        #
                     
        for i in range(len(self.kps)):
            year = t_utc[i].year
            if year < 1955: year = 1955
            if year > 2015: year = 2015
            iy = (year - 1955)/5
            ir = (year - 1955)%5
            cyear = self.cyears [iy]
            iu = int((t_utc[i].hour+t_utc[i].minute/60. + 1.5)/3.)
            # UT =1 corresponds to ut: 1.5 - 4.5 hrs
            if iu > 7: iu = 0
            cut = self.cuts[iu]
            ckp = self.ckps[kp]
            self.dbMgr.getMap(cyear,ckp,cut)
            mkey = cyear1+ckp+cut
            #
            isotime = [cyear+'-01-01T'+cut+':00:00']  # time-date of the map
            #atime = spt.Ticktock(['2002-02-02T12:00:00'],'ISO')
            atime = spt.Ticktock(isotime,'ISO')
            #y = spc.Coords([[3,0,0],[2,0,0],[1,0,0]],'GEO','car')
            aposi = self.positions[i]
            aposi.radi = [aposi.Re+450.]
            t_dic = ib.get_Lm(atime,aposi,[90.],'T89')
            lm = t_dic['lm'].flatten()[0]
            lon = aposi.long
            lat = aposi.lati 
            rc = getRC450km(mkey,lon,lat,lm)
            #
            w =(ir + (self.times[i].DOY + t_utc[i].hour/24.)/365.)/5. 
            # the next map at +5 years if required
            #
            if 1955 < year < 2015 and w < 1.:
                cyear = self.cyears [iy+1]
                self.dbMgr.getMap(cyear,ckp,cut)
                mkey = cyear+ckp+cut
                isotime = [cyear+'-01-01T'+cut+':00:00']  # time-date of the map
                atime = spt.Ticktock(isotime,'ISO')
                t_dic = ib.get_Lm(atime,aposi,[90.],'T89')
                lm1 = t_dic['lm'].flatten()[0]
                rc1 = getRC450km(mkey,lon,lat,lm1)
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
            radist = aposi.radi/aposi.Re                    
            rcorr = log(radist*radist)/14.   # Don used 11. but 14. is better
            rcr -= rcorr 
            if rcr < 0.: rcr = 0.
            #
            rclist.appeend(rcr)
            
        return rclist          
            
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
        i, j = getGridIdx(lon, lat)
        #get the weights
        wl,wr,wt,wb = getWeights(lon,lat)
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
        wl,wr, wt, wb = getWeights(lon,lat)
        # 
        rclm_l = wt*rclm_LT + wb*rclm_LB 
        rcml_r = wt*rclm_RT + wb*rclm_RB
        #
        rclm = wl*rclm_l + wr*rclm_r
        # 
        return rclm/lm**2
     
    def getWeights(lon,lat):
        
        # weights in longitude
        wl = (lon%5)/5. # left side of the box 
        wr = 1.0 - wl # right side of the box
        # weights in latitude
        wb = (-lat%5)/5. # bottom side of the box 
        wt = 1. - wb 
        return wl,wr, wt, wb
        
              
    def getGridIdx(self,lon, lat):
        ix = int(lon)/5
        iy = 18 + int(-lat)/5
        return ix, iy 

    def getEMLat():
        '''    
        calculate the equivalent magnetic latitude of the given locations using SpacePy
    
            Get Corrected geomagnetic latitude (GMLATC) at sub-satellite point
            Then, get Invariant latitude at satellite position
                       INVARIANT LAT = ACOS(1.0/SQRT(L))
            Select the smaller value for magnetic latitud
                   We will always use absolute value of equivalent magnetic latitude
        Inputs:
            times:
            positions:
        
        Returns:
        
        ''' 
        # calculate the corrected magnetic latitude
        # to to reset the altitudes = Re    
        radia = np.empty(len(self.times))
        radia.fill (self.coords.Re)
        positions=self.coords
        positions.radi = radia
        #
        # GDZ -> MAG
        mpos = ib.coord_trans(positions,'MAG','sph')
        #
        # mpos.lati are the magnetic latitude in degrees. Note this is not the same as 
        # the corrected geomagnetic latitude, but the difference should be small
        gmlatcr = mpos.lati/57.2957795 # convert to radians
        emlats = []
        for i in range(len(posistions)):
            glmdar = 0.0
            if (self.lm[i] > 1.): glmdar = acos(1.0/sqrt(self.lm[i]))
            if abs(gmlatcr)< abs(glmdar):  glmdar = abs(gmlatcr)
            emlats.append(glmdar)
        #
        return emlats    
    
import Spenvis
#
class SpenvisCSVFileHandler():
    '''
    classdocs
    This class allows to read/write Spenvis CSV file
    It make use of the C++_to_python interface Spenvis.so 
    '''
     
    def __init__(self):
        '''
        Constructor
        '''
        self.list_blocks=[]
         
    def ReadFile(self,file_name):
        theCSVCollection=Spenvis.SpenvisCSVCollection()
        self.list_blocks=[]
        if (not os.path.exists(file_name)):
            print("The CSV file that you want to open does not exist!") 
            return False
        list_CSV=theCSVCollection.ReadFile(file_name)
        for name in list_CSV:
            self.list_blocks+=[SpenvisCSVBlock(list_CSV[name])]
        return True
     
    def ReadSpenvisOrbitFile(self,file_name):
        aBool=self.ReadFile(file_name)
        if (not aBool):
            print("Problem when reading the file ", file_name)
            return None,None
        aBool=False
        if (len(self.list_blocks)==1):
            meta_var_bl1=self.list_blocks[0].MetaVariables
            if "MOD_ABB" in meta_var_bl1:
                if (meta_var_bl1["MOD_ABB"]['values'][0].replace("'","")=="ORB"):
                    aBool=True
        if (not aBool):
            print("The file ", file_name, "is not a SPENVIS orbit file")
            return None,None
        orbit_data=self.list_blocks[0].Variables
        return orbit_data
     
class SpenvisCSVBlock():
    '''
    This class is copied from CIRSOS
      
    '''
     
    def __init__(self,aSpenvisBlock):
        '''
        Constructor
        ''' 
         
        #Get the variable names and characteristic
        ######################
        self.Comments=[]
        for i in range(aSpenvisBlock.GetNumTextLines()):
            self.Comments+=[aSpenvisBlock.GetComment(i)]
        #Get the meta variable
        ######################
        self.MetaVariables=dict()
         
        for i in range(aSpenvisBlock.GetNumMetaVariables()):
            name=aSpenvisBlock.GetMetaVariableName(i)
            line=aSpenvisBlock.GetMetaVarLine(name)
            words=line.split(",")
            new_words=[]
            ii=0
            while ii<len(words):
                word=words[ii]
                if (word[0] =="'" and word[-1]!="'"):
                    word+=","+words[ii+1]
                    ii+=1
                new_words+=[word]   
                ii+=1    
            ndim=int(new_words[1])
            values=[]
            if ndim >0:
                for i in range(ndim):
                    values+=[float(new_words[2+i])]    
            else:
                for i in range(abs(ndim)):
                    values+=[new_words[2+i]]        
                 
            self.MetaVariables[name]={"dim":ndim,"values":np.array(values)}
             
        #Get the variable names and characteristic
        ######################
        self.Variables=dict()
         
        nrows=aSpenvisBlock.GetNumDataLines()
        ncol=aSpenvisBlock.GetNumDataColumns()
        data=np.array([])
        if (nrows >0 and ncol>0):
            data=np.zeros((nrows,ncol))
        for i in range(aSpenvisBlock.GetNumVariables()):           
            line=aSpenvisBlock.GetVariableLine(i) 
            words=line.split(",")
            new_words=[]
            ii=0
            while ii<len(words):
                word=words[ii]
                if (word[0] =="'" and word[-1]!="'"):
                    word+=","+words[i+1]
                    ii+=1
                new_words+=[word]   
                ii+=1    
            ndim=int(new_words[2]) 
            name=  new_words[0].replace("'","") 
            data=np.array([])
            if (nrows >0 and ndim>0):
                data=np.zeros((nrows,ndim))    
            for j in range(nrows):
                data[j,:]=np.array(aSpenvisBlock.GetDataRecordForPython(name,j))
            if (ndim ==1):
                data=np.reshape(data,nrows)
            self.Variables[name]={"dim":ndim,"desc":new_words[3],"unit":new_words[1],"data":data}
             
#
def calculateRInv(B,L):
        '''
        calculte the invariant radial distance (R) and the magnetic latitude lambda based on the method of 
       
        "Roberts, C. S. (1964), Coordinates for the study of particles trapped in
        the Earth’s magnetic field: A method of converting from B, L to R, l
        coordinates, J. Geophys. Res., 69, 5089-- 5090."
        
        Inputs: 
            B, L
        Ouputs:
            R, lambda
            
        '''
        a = [1.25992106, -0.19842592, -0.04686632, -0.01314096, -0.00308824, 0.00082777, -0.00105877, 0.00183142]
        Md = 31165.3  #nT*Re^3
        if (L<0. ): return -1., -1.
        p = pow(pow(L,3.)*B/Md,-1./3.)
        if (p<0. or p > 1.): return -1., -1.
        
        s = 0.
        for i in range(8):
            s += a[i]*pow(p,i)
        ps = p*s 
        if ps > 1.:
            # print ("ps = ",ps)
            ps = 1.
        lamb = degrees(acos(sqrt(ps)))
        R = L*ps
        return R, lamb


def getExpfact(rcv, mlat, rc):
    '''
C----------------------------------------------------------------------
c     SUBROUTINE EXPFAC(RADIUS, MLAT, E, Z, A, F, N, ISTORM)
C     computes an exposure factor at B,L for particles Z,A 
C     arriving from interplanetary space with energies E(N)
C     by averaging over arrival directions. 
C     Exposure factor returned in F(N).
C     First implementation E.J. Daly ESA/ESTEC/WMA  9/87
c     Updated by F. Lei 03/14
c     1) change to rigidity based in stead of energy, so no need of Z and A
c     2) the Rcv is given as input, no need for RADIUS and ISTORM
c     3) use omcao2
c
    '''
    # table of One-Minus-Cos-Angles-Over-2 :
    omcao2 = np.array([0., .067, .146, .25, .5, .75, .854, .933, 1.])
    ang = np.array([0.01, .5236, .785, 1.047, 1.571, 2.094, 2.356, 2.618, 3.1416])
    nangle = 9
    #
    N = len(rc)
    fac = np.empty(N)
    cut = np.empty[nangle]
    #
    if (rcv <= 0.1): 
        fac.fill(1.0)
    else:
        for ia in range(nangle):
            cutf = cutfac(mlat,ang[ia])
            cut[ia] =4.*rcv*cutf
        for ir in range(N):
            fac[ir] = 0.
            # is there an value in the array above the lowest (Eastward, from West!) cutoff? 
            if (rc[ir] >= cut[0]):
                fac[ir] = 1.0  
                for ia in range(1,nangle):
                    # find the angular location where the cutoff goes over the rc:
                    if (rc[ir] <= cut[ia]):
                        fac[ir] = omcao2[ia-1] + (rc[ir]-cut[ia-1])*(omcao2[ia]-omcao2[ia-1]) \
                        *(cut[ia]-cut[ia-1])
                        break
            
    return 
                    
def facshadow(R):
    '''

   This is a correction factor for the earth's shadow on 
   the spacecraft according to simple geometrical optics.
   
   Inputs:
       R = dadius in 
   Return:
   
    '''
    fac = 1. - .5 * (1.-sqrt(R**2-1.)/R)

    return fac
      

def cutfac (fmlat, a):
    '''
    '''
    cosa = cos(a)
    cosl = cos(fmlat)
    cutfac = 1.0/(1.0+sqrt(1.0+cosa*cosl**3))**2
    return cutfac        
        
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
    csvMgr = SpenvisCSVFileHandler()
    orbdata = csvMgr.ReadSpenvisOrbitFile('test_orbit.txt')
    print (orbtata)
#    csvMgr.ReadFile('test2.csv')
#     csvMgr.list_blocks
    mapMgr = MapDB()
#    lm, rc, rclm2 = mapMgr.getMap('2000','9','12')
#    plotscatter(lm,rc)
#    plotmap_c(rc)

if __name__ == '__main__':
    main() 