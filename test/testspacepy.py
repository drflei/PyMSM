import spacepy.time as spt
import spacepy.coordinates as spc
import spacepy.irbempy as ib
import spacepy.omni as om
import datetime

import numpy as np

def testspacepy():
    '''
    This is a test of the spacepy module

    '''
    cyear ='2005'
    cut = '12'
#    isotime = [cyear+'-01-01T'+cut+':00:00.000000']  # time-date of the map
    # case 1 - iso 
    try: 
        atime = spt.Ticktock(['2002-02-02T12:00:00.000000'],'ISO')
        print('case1 - ',atime.getDOY())
    except:
        print('case1 - failed ')
        pass
    # case 2 - d-m-y
    try:
        atime = spt.Ticktock(['01-01-2013'], lambda x: datetime.datetime.strptime(x, '%d-%m-%Y'))
        print('case2 - ',atime.getDOY())
    except:
        print('case2 - failed ')
        pass
    # case 3 - mjd 
    try: 
        atime = spt.Ticktock([55100.2], 'MJD')
        print('case3 - ',atime.getDOY())
    except:
        print ('case3 - failed!')
        pass

#    atime = spt.Ticktock(isotime)
    aposi = spc.Coords([450,22,351],'GDZ','sph')
    aomni = om.get_omni(atime)
    aomni['Kp'] =[3.]
    try:
        t_dic = ib.get_Lm(atime,aposi,90,'T89',omnivals=bomni)
    except:
        t_dic = ib.get_Lm(atime,aposi,90,'T89')
    lm = abs(t_dic['Lm'].flatten()[0])
    print(lm,'Test1 on Lm completed!')
    #
    N = 11*1 
    times = np.empty(N,dtype='object')
    kps = np.empty(N,dtype=int)
    mjd = 55000.5
    for i in range(N):
        mjd += 0.00000001*i
        times[i] = mjd  
#    times.fill('2019-02-02T12:00:00')
    kps.fill(1)      
    coords =[]
    alti = 1000.
    for i in range(-90,91,18): 
        for j in range(0,360,360):
            coords.append([alti, i, j ]) 
    t = spt.Ticktock(times,'MJD')
    y = spc.Coords(coords,'GDZ','sph') 
    aomni = om.get_omni(t)
    #aomni['Kp'] =kps
    t_dic = ib.get_Lm(t,y,90,'T89',omnivals=aomni)
#    t_dic = ib.get_Lm(t,y,90,'T89')
    print (np.nan_to_num(np.abs(t_dic.pop('Lm')),nan=30.)) 

     

if __name__ == '__main__':
    testspacepy() 
