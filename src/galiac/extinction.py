'''


@author: jmolgo
'''

from os import path
import numpy as np
import math
import scipy.ndimage

class Drimmel(object):
    '''
    Drimmel extinction model
    '''
    
    delta_la_d=(0.2, 0.2, 0.02)
    delta_la_s=(0.2, 0.2, 0.02)
    delta_la_o=(0.05, 0.05, 0.02)
    
    delta_lo_d=(0.05, 0.05, 0.02)
    delta_lo_o=(0.02, 0.02, 0.02)
    
    moddir = path.dirname(__file__)
    avdisk_fname = path.join (moddir, 'avdisk.npy') 
    avspir_fname = path.join (moddir, 'avspir.npy')
    avori_fname = path.join (moddir, 'avori.npy')
    avori2_fname = path.join (moddir, 'avori2.npy')
    avdloc_fname = path.join (moddir, 'avdloc.npy')
    ncomp_fname = path.join (moddir, 'ncomp.npy')
    glon_fname = path.join (moddir, 'glon.npy')
    glat_fname = path.join (moddir, 'glat.npy')
    rfact_fname = path.join (moddir, 'rfact.npy')
    
    avdisk = np.load(avdisk_fname)
    avspir = np.load(avspir_fname)
    avori = np.load(avori_fname)
    avori2 = np.load(avori2_fname)
    avdloc = np.load(avdloc_fname)
    ncomp = np.load(ncomp_fname)
    glon = np.load(glon_fname)
    glat = np.load(glat_fname)
    rfact = np.load(rfact_fname)
    
    rfact_d = np.ones (rfact.size)
    rfact_s = np.ones (rfact.size)
    rfact_o = np.ones (rfact.size)
    
    ind_d = np.where (ncomp==1)
    rfact_d[ind_d] = rfact[ind_d]
    ind_s = np.where (ncomp==2)
    rfact_s[ind_s] = rfact[ind_s]
    ind_o = np.where (ncomp==3)
    rfact_o[ind_o] = rfact[ind_o]
    
    Al_Av = {'U':    1.531, 
             'B':    1.324,
             'V':    1.000,
             'R':    0.748,
             'I':    0.482,
             'J':    0.282,
             'H':    0.175,
             'K':    0.112,
             'L':    0.058,
             'M':    0.023,
             'N':    0.052,
             '8.0':  0.020,
             '8.5':  0.043,
             '9.0':  0.074,
             '9.5':  0.087,
             '10.0': 0.083,
             '10.5': 0.074,
             '11.0': 0.060,
             '11.5': 0.047,
             '12.0': 0.037,
             '12.5': 0.030,
             '13.0': 0.027}

    def __init__(self, Rsun, Zsun):
        '''
        Constructor
        '''
        self.Rsun = Rsun
        self.Zsun = Zsun
    
    def galToGlobalXY (self, r, l, b):
        return (-(self.Rsun-r*np.cos(b)*np.cos(l)), 
                r*np.cos(b)*np.sin(l), 
                r*np.sin(b)+self.Zsun)
    
    def galToLocalXY (self, r, l, b):
        return (r*np.cos(b)*np.cos(l), 
                r*np.cos(b)*np.sin(l), 
                r*np.sin(b)+self.Zsun)
    
    def xyzToGrid (self, x, y, z, delta, nele, avori=False):
        xfac = 2.5 if avori else 0.5
        return (np.minimum(x/delta[0] + xfac*(nele[0]-1), nele[0]-1.0001), 
                np.minimum(y/delta[1] + 0.5*(nele[1]-1), nele[1]-1.0001), 
                np.minimum(z/delta[2] + 0.5*(nele[2]-1), nele[2]-1.0001))
    
    def lin3dInterp (self, grid, i, j, k):
        return scipy.ndimage.map_coordinates( grid, np.array([i,j,k]), order=1)
    
    def adjustDistance (self, r, lon, lat, ori=False):
        if not ori:
            d1=0.5/abs(math.sin(lat))-self.Zsun/math.sin(lat)
            d2=15.0/abs(math.cos(lon))+self.Rsun/math.cos(lon)
            d2=d2/abs(math.cos(lat))
            d3=15.0/abs(math.sin(lon))
            d3=d3/abs(math.cos(lat))
        else:
            d1=0.5/abs(math.sin(lat))-self.Zsun/math.sin(lat)
            d2=3.75/abs(math.sin(lon))
            d2=d2/abs(math.cos(lat))
            d3=2.375/abs(math.cos(lon)) if (math.cos(lon)>0) else 1.375/abs(math.cos(lon))
            d3=d3/abs(math.cos(lat))
    
        mind = min(d1, d2, d3)
        return np.where (r<mind, r, mind)
    
    def lb2pix (self, l, b):
        x, y, z = self.ll2uv (l, b)
        xe, ye, ze = self.gal2ecliptic (x, y, z)
        return self.pixNumber (xe, ye, ze)

    def ll2uv (self, l, b):
        return math.cos(b) * math.cos(l), math.cos(b) * math.sin(l), math.sin(b)

    def gal2ecliptic (self, x, y, z):
        matr2 = np.array([[-0.054882486, -0.993821033, -0.096476249], [0.494116468, -0.110993846,  0.86228144], [-0.867661702, -0.000346354,  0.497154957]])
        matr = np.transpose(matr2)
        ov = np.empty (3)
        iv = np.array([x,y,z])
        for i in range (3):
            ov[i] = 0.0
            for j in range (3):
                ov[i] = ov[i]+iv[j]*matr[i][j]
        return ov[0], ov[1], ov[2]
    
    def pixNumber (self, v1, v2, v3, resolution=9):
        ix = np.zeros(129)
        iy=np.zeros(129)
        for i in range (1,129):
            j = i - 1
            k = 0
            ip = 1
            while (j != 0):
                idd = j%2
                j = j/2
                k = ip * idd + k
                ip = ip * 4
            ix[i] = k
            iy[i] = 2*k
        face2, x, y = self.axisxy(np.array([v1,v2,v3]))
        face4 = face2
        i = int(2.**14 * x)
        j = int (2.**14 * y)
        if (i > 16383):
            i = 16383
        if (j > 16383):
            j = 16383
        ih = i>>7 
        ih1 = ih    
        il = i - (ih1<<7)
        jh = j>>7
        jh1= jh    
        jl = j - (jh1<<7)
        pixel = (int(face4)<<28) + ix[il+1] + iy[jl+1] + ((int(ix[ih+1] + iy[jh+1]))<<14)
    
        #'Pixel' now contains the pixel number for a resolution of 15.  To
        #convert to the desired resolution, (integer) divide by 4 to the power
        #of the difference between 15 and the given resolution:
        
        pixel = int(pixel)>>(30-2*resolution)  
        return pixel

    def axisxy(self, v):
        v1=v[0]
        v2=v[1]
        v3=v[2]
        av1=abs(v[0])
        av2=abs(v[1])
        av3=abs(v[2])
        if (av3 > av2):
            if (av3 > av1):
                if (v3 > 0):
                    nface = 0
                    eta = -v1/v3
                    xi = v2/v3
                else:
                    nface = 5
                    eta=-v1/v3
                    xi=-v2/v3
            else:
                if (v1 > 0):
                    nface=1
                    xi=v2/v1
                    eta=v3/v1
                else:
                    nface=3
                    eta=-v3/v1
                    xi=v2/v1
        else:
            if(av2 > av1):
                if (v2 > 0):
                    nface=2
                    eta=v3/v2
                    xi=-v1/v2
                else:
                    nface=4
                    eta=-v3/v2
                    xi=-v1/v2
            else:
                if (v1 > 0):
                    nface=1
                    xi=v2/v1
                    eta=v3/v1
                else:
                    nface=3
                    eta=-v3/v1
                    xi=v2/v1
        x, y = self.incube (xi, eta)
        x = (x+1)/2
        y = (y+1)/2
        return nface, x, y
    
    def incube(self, beta, BETA):
        gstar=1.37484847732 
        G=-0.13161671474
        M=0.004869491981
        W1=-0.159596235474
        C00=0.141189631152
        C10=0.0809701286525
        C01=-0.281528535557
        C11=0.15384112876
        C20=-0.178251207466
        C02=0.106959469314
        D0=0.0759196200467
        D1=-0.0217762490699
        R0=0.577350269
        AA=beta**2
        BB=BETA**2
        A4=AA**2
        B4=BB**2
        ONMAA=1.-AA
        ONMBB=1.-BB
        
        gstar_1 = 1. - gstar
        m_g     = M - G
        c_comb  = C00 + C11*AA*BB
        
        X = beta * (gstar + AA * gstar_1 + ONMAA * (BB * (G + (m_g)*AA + ONMBB * (c_comb + C10*AA + C01*BB +C20*A4 + C02*B4)) +AA * (W1 - ONMAA*(D0 + D1*AA))))
    
        Y = BETA * (gstar + BB * gstar_1 + ONMBB * (AA * (G + (m_g)*BB +ONMAA * (c_comb + C10*BB + C01*AA +C20*B4+C02*A4)) + BB * (W1 - ONMBB*(D0 + D1*BB))))
        return X, Y    

    def compute_extinction (self, coords, filt='K', params={'useCobe':True, 'cobePix':None}):
        pix = params['cobePix']
        useCobe = params['useCobe']
        l = coords.l
        b = coords.b
        r = coords.r
        if pix is None and useCobe is True:
            pix = self.lb2pix (l, b)
        if useCobe:
            fd = self.rfact_d[pix]
            fs = self.rfact_s[pix] 
            fo = self.rfact_o[pix]
        else:
            fd = 1.0
            fs = 1.0
            fo = 1.0
                
        rr = self.adjustDistance (r, l, b)    
        x, y, z = self.galToGlobalXY (rr, l, b)
        
        #disk
        i,j,k = self.xyzToGrid (x, y, z, self.delta_la_d, self.avdisk.shape)
        Ad = self.lin3dInterp (self.avdisk, i, j, k)
    
        #spiral
        i,j,k = self.xyzToGrid (x, y, z, self.delta_la_s, self.avspir.shape)
        As = self.lin3dInterp (self.avspir, i, j, k)
    
        #orion arm        
        rr = self.adjustDistance (r, l, b, True)
        x, y, z = self.galToGlobalXY (rr, l, b)
        i,j,k = self.xyzToGrid (x, y, z, self.delta_la_o, self.avori.shape, True)
        Ao = self.lin3dInterp (self.avori, i, j, k)
    
        #correct for wavelength
        wavelengthCorrection = self.Al_Av[filt]
    
        return wavelengthCorrection*(fd*Ad + fs*As + fo*Ao)
    
