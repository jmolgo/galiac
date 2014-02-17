'''

@author: jmolgo
'''

import numpy as np
import math 
import pickle
from scipy import integrate
from os import path
import extinction

class Coordinates(object):
    '''
    This class encapsules the concept of galactic
    spatial coordinates.
    '''
    
    Rsun = 8.0
    
    def __init__(self, l, b, r, in_radians=False):
        ''' 
        Constructor.
        
        @param l: galactic longitude
        @param b: galactic latitude
        @param r: heliocentric distance, in kpc
        @param radians: l and b are given in radians if this parameter 
                        is True, and in degrees otherwise
        ''' 
        
        self._r = r
        self._l = np.empty (r.size)
        self._b = np.empty (r.size)
        if in_radians:
            self.l = l
            self.b = b
            self._l.fill (l)
            self._b.fill (b)
        else:
            self.l = math.radians(l)
            self.b = math.radians(b)
            self._l.fill (np.radians(l))
            self._b.fill (np.radians(b))
        
        #create the derived attributes, initializing them
        #to None. A lazy initialization will be performed
        #by the properties
        self._R = None
        self._theta = None
        self._x = None
        self._y = None
        self._z = None   
    
    @property
    def r (self):
        return self._r
    
    @property
    def x (self):
        ''' Cartesian galactocentric X coordinates ''' 
        if self._x is None:
            self._x = self.Rsun - self._r*np.cos(self._b)*np.cos(self._l)
        return self._x
    
    @property
    def y (self):
        ''' Cartesian galactocentric Y coordinates '''
        if self._y is None:
            self._y = self._r*np.cos(self._b)*np.sin(self._l)
        return self._y
    
    @property
    def z (self):
        ''' Cartesian galactocentric Z coordinates '''
        if self._z is None:
            self._z = self._r*np.sin(self._b)
        return self._z
    
    @property
    def R (self):
        ''' Cylindrical galactocentric radius '''
        if self._R is None:
            self._R = (self.x**2 + self.y**2)**0.5
        return self._R
    
    @property
    def theta(self):
        ''' Cylindrical galactocentric azimuth angle '''
        if self._theta is None:
            self._theta = np.arctan2(self.y, self.x)
        return self._theta
    
    def r_pc(self):
        ''' Returns the heliocentric distances, in parsecs '''
        return self._r*1e3
    
    
class ModelParameters(object):
    '''
    This class provides a repository for all the galactic
    model parameters.
    
    The model_parameters class implements the functionalities
    of reading or writing the model parameters as attributes, as a 
    dictionary, as files, etc.
    '''
    
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
    
    def update (self, new_parms):
        self.__dict__.update (new_parms)
    
    def from_file (self, file_name):
        f = open(file_name, 'r')
        parms = pickle.load(f)
        self.update (parms)
    
    def to_file (self, file_name):
        with open(file_name, 'w') as f:
            pickle.dump(self.__dict__, f)

class GalacticComponent(object):
    '''
    This is the base class for all the galactic component classes
    '''
    
    def __init__(self, params):
        self._lf = {}
        self._cumulated_lf = None        
        self.current_passband = ''
        self.params = params
        self.default_parameters = {}
        self.enabled = True
    
    def density (self, coords):
        return np.zeros(coords.shape)
    
    def set_passband (self, band):
        if band != self.current_passband:
            self._cumulated_lf = None
            self.luminosity_function = self._lf[band]
            self.current_passband = band
    
    @property
    def cum_luminosity_function(self):
        if self._cumulated_lf is None:
            self._cumulated_lf = integrate.cumtrapz (self.luminosity_function[1], self.luminosity_function[0])
        return self._cumulated_lf
    
    @property
    def param_names(self):
        return self.default_parameters.keys()

class ThinDisk(GalacticComponent):
    '''
    This class implements the Thin Disk galactic component, with flare and warp
    '''

    param_names = ['hr_thin', 'hz_thin_sun', 'hr_hole', 'hz_hole_lin', 'hz_hole_sq', 'Rflare', 'Cwarp', 'Ewarp', 'theta_warp']
    def __init__(self, params):
        '''
        Constructor.
        
        @param params: The parameters list object
        '''
         
        super(ThinDisk, self).__init__(params)
        self.default_parameters = {
                                   'hr_thin': 1.97
                                   ,'hz_thin_sun': 0.285
                                   ,'hr_hole': 3.74
                                   ,'hz_hole_lin': 0.21
                                   ,'hz_hole_sq': 0.056 
                                   ,'Rflare': 15.0
                                   ,'Zsun': 15.0
                                   ,'Cwarp': 2.1e-19 #2.1x10^(-19) pc
                                   ,'Ewarp': 5.25
                                   ,'theta_warp': math.radians(-5.0)
                                   }
        
        
        self.load_luminosity_functions()
    
    def density (self, coords):
        '''
        This function returns the stellar disk density in the given coordinates
        
        @param coords: The Coordinates object with the galactic locations where
                       the density shall be computed 
        '''
        
        p = self.params
        R = coords.R
        z = coords.z
        
        #compute scale height
        hz = self.disk_hz (coords)
        
        #compute warp
        zw = self.disk_warp (coords)
        
        #compute density
        return (p.A
               *(p.hz_thin_sun/hz)
               *np.exp(p.Rsun/p.hr_thin + p.hr_hole/p.Rsun)
               *np.exp(-R/p.hr_thin-p.hr_hole/R)
               *np.exp(-abs(z-zw)/hz))
        
    def disk_hz (self, coords):
        '''
        Returns the scale height of the disc in the given coordinates
         
        @param coords: The Coordinates object with the galactic locations where
                       the density shall be computed 
        '''
        p = self.params
        R = coords.R
        return p.hz_thin_sun*(1+p.hz_hole_lin*(R-p.Rflare)+p.hz_hole_sq*(R-p.Rflare)**2)
    
    def disk_warp (self, coords):
        #saturate for large R
        p = self.params
        theta = coords.theta
        R = np.where (coords.R<15.0, coords.R, 15.0)
        return (p.Cwarp*((R*1e3)**p.Ewarp)*np.sin(theta-p.theta_warp)+p.Zsun)*1e-3
    
    def load_luminosity_functions (self, lf_files={}):
        if len(lf_files) == 0:
            bands = ['B', 'V', 'J', 'H', 'K', '12', '25']
            moddir = path.dirname(__file__)
            lf_files = dict([ (band, path.join(moddir, 'disk_lf_'+band+'.npy')) for band in bands ])
        self._lf = {}
        for band in bands:
            self._lf[band] = np.load(lf_files[band])
        self.set_passband(bands[0])

class ThickDisk(ThinDisk):
    '''
    This class implements the Thick Disk galactic component, with flare and warp
    '''
    
    param_names = ['hr_thick', 'hz_thick_sun', 'hz_hole_lin_th', 'hz_hole_sq_th', 'RflareThick']

    def __init__(self, params):
        '''
        Constructor.
        
        @param params: The parameters list object
        '''
         
        super(ThickDisk, self).__init__(params)
        self.default_parameters = {
                                   'hr_thick': 3.0
                                   ,'hz_thick_sun': 0.631
                                   ,'hz_hole_lin_th': 0.21
                                   ,'hz_hole_sq_th': 0.056 
                                   ,'RflareThick': 8.0
                                   ,'Zsun': 15.0
                                   ,'Cwarp': 2.1e-19 #2.1x10^(-19) pc
                                   ,'Ewarp': 5.25
                                   ,'theta_warp': math.radians(-5.0)
                                   ,'thickDiskFactor': 0.09
                                   }
        
    def density (self, coords):
        '''
        This function returns the stellar disk density in the given coordinates
        
        @param coords: The Coordinates object with the galactic locations where
                       the density shall be computed 
        '''
        
        p = self.params
        R = coords.R
        z = coords.z
        
        #compute scale height
        hz = self.disk_hz (coords)
        
        #compute warp
        zw = self.disk_warp (coords)
        
        #compute density
        return (p.thickDiskFactor
               *p.A
               *(p.hz_thick_sun/hz)
               *np.exp(p.Rsun/p.hr_thick + p.hr_hole/p.Rsun)
               *np.exp(-R/p.hr_thick-p.hr_hole/R)
               *np.exp(-abs(z-zw)/hz))
        
    def disk_hz (self, coords):
        '''
        Returns the scale height of the disc in the given coordinates
         
        @param coords: The Coordinates object with the galactic locations where
                       the density shall be computed 
        '''
        p = self.params
        R = coords.R
        return np.where (R>p.RflareThick, 
                         p.hz_thick_sun*(1+p.hz_hole_lin_th*(R-p.RflareThick)+p.hz_hole_sq_th*(R-p.RflareThick)**2), 
                         p.hz_thick_sun)
    

class Bulge(GalacticComponent):
    '''
    This class implements the Bulge galactic component
    '''

    def __init__(self, params):
        '''
        Constructor.
        
        @param params: The parameters list object
        '''
         
        super(Bulge, self).__init__(params)
        self.default_parameters = {
                                    'bulgeAngle': math.radians(27.0)
                                   ,'bulgeAxialRatio': np.array([1.0, 0.49, 0.40])
                                   ,'bulgeCentralDensity': 6.6e9 #6.6 stars pc-3
                                   ,'bulgeScaleLength': 0.740*0.866 # scale length of 740 parsecs, plus
                                                                    # correction factor from t_ellipsoid
                                   }
        
        self.load_luminosity_functions()
    
    def density (self, coords):
        '''
        This function returns the bulge stellar density in the given coordinates
        
        @param coords: The Coordinates object with the galactic locations where
                       the density shall be computed 
        '''
        p = self.params
        alpha = p.bulgeAngle
        
        x = coords.x
        y = coords.y
        z = coords.z
         
        #rotate the coordinates to the reference frame defined by the bulge elipsoid axes        
        x1 =  math.cos(alpha)*x + math.sin(alpha)*y
        x2 = -math.sin(alpha)*x + math.cos(alpha)*y
        x3 = z #third axis is unrotated, perpendicular to the plane

        #compute the star density, in parametric form
        #use the BOXY model of the bulge
        t = (x1**4 + (x2/p.bulgeAxialRatio[1])**4 + (x3/p.bulgeAxialRatio[2])**4)**0.25
        return p.bulgeCentralDensity*np.exp(-t/p.bulgeScaleLength)
            
    def load_luminosity_functions (self, lf_files={}):
        if len(lf_files) == 0:
            bands = ['B', 'V', 'J', 'H', 'K', '12', '25']
            moddir = path.dirname(__file__)
            lf_files = dict([ (band, path.join(moddir, 'bulge_lf_'+band+'.npy')) for band in bands ])
        self._lf = {}
        for band in bands:
            self._lf[band] = np.load(lf_files[band])
        self.set_passband(bands[0])

class Bar(GalacticComponent):
    '''
    This class implements the Bar galactic component
    '''

    def __init__(self, params):
        '''
        Constructor.
        
        @param params: The parameters list object
        '''
         
        super(Bar, self).__init__(params)
        self.default_parameters = {
                                    'barAxialRatio': np.array([1.0, 0.10, 0.03]) #[Lop 2007]
                                   ,'barMajorSemiaxis': 3.9 #major axis is 7.2 kpc [Lop 2007]
                                   ,'barAngle': math.radians(43.0) # [Lop 2007]
                                   ,'barCentralDensity': 3.0e9 #3.5 stars/pc-3 [Lop 2007]
                                   }
        
        self.load_luminosity_functions()
    
    def density (self, coords):
        '''
        This function returns the bar stellar density in the given coordinates
        
        @param coords: The Coordinates object with the galactic locations where
                       the density shall be computed 
        '''
        p = self.params
        x = coords.x
        y = coords.y
        z = coords.z
        
        #rotate the coordinates to the reference frame defined by the bar major axis
        alpha = p.barAngle
        xb =  math.cos(alpha)*x + math.sin(alpha)*y
        yb = -math.sin(alpha)*x + math.cos(alpha)*y
        zb = z #third axis is unrotated, perpendicular to the plane

        #bar from [Lop 2007]
        rhobar = p.barCentralDensity*np.exp(-yb**2/(2*(p.barAxialRatio[1]*p.barMajorSemiaxis)**2))*np.exp(-zb**2/(2.*(p.barAxialRatio[2]*p.barMajorSemiaxis)**2))
        #truncate the bar
        rhobar[np.where(abs(xb) > p.barMajorSemiaxis)] = 0.0
        rhobar[np.where(abs(yb) > 3*p.barAxialRatio[1]*p.barMajorSemiaxis)] = 0.0
        rhobar[np.where(abs(zb) > 3*p.barAxialRatio[2]*p.barMajorSemiaxis)] = 0.0
        return rhobar
            
    def load_luminosity_functions (self, lf_files={}):
        if len(lf_files) == 0:
            bands = ['B', 'V', 'J', 'H', 'K', '12', '25']
            moddir = path.dirname(__file__)
            lf_files = dict([ (band, path.join(moddir, 'bar_lf_'+band+'.npy')) for band in bands ])
        self._lf = {}
        for band in bands:
            self._lf[band] = np.load(lf_files[band])
        self.set_passband(bands[0])

class Halo (GalacticComponent):
    '''
    This class implements the Halo galactic component
    '''

    def __init__(self, params):
        '''
        Constructor.
        
        @param params: The parameters list object
        '''
         
        super(Halo, self).__init__(params)
        self.default_parameters = {
                                    'haloFactor': 1.4e-3
                                   ,'haloExp1': 10.093
                                   ,'haloExp2': 7.0/8
                                   }
        
        self.load_luminosity_functions()
    
    def density (self, coords):
        '''
        This function returns the halo stellar density in the given coordinates
        
        @param coords: The Coordinates object with the galactic locations where
                       the density shall be computed 
        '''
        p = self.params
        R = coords.R
        z = coords.z
        
        Rsp = (R**2 + 2.52*z**2)**0.5
        Reff = np.where (Rsp>5.0, Rsp, 5.0)
        return p.haloFactor*p.A*np.exp(p.haloExp1*(1-(Reff/p.Rsun)**0.25))/((Reff/p.Rsun)**(p.haloExp2))
                    
    def load_luminosity_functions (self, lf_files={}):
        if len(lf_files) == 0:
            bands = ['B', 'V', 'J', 'H', 'K', '12', '25']
            moddir = path.dirname(__file__)
            lf_files = dict([ (band, path.join(moddir, 'halo_lf_'+band+'.npy')) for band in bands ])
        self._lf = {}
        for band in bands:
            self._lf[band] = np.load(lf_files[band])
        self.set_passband(bands[0])

class SpiralArms (GalacticComponent):
    '''
    This class implements the Spiral Arms galactic component
    '''

    def __init__(self, params):
        '''
        Constructor.
        
        @param params: The parameters list object
        '''
         
        super(SpiralArms, self).__init__(params)
        self.default_parameters = {
                                    'hz_arms_sun': 0.09
                                   ,'armInitialRadius': 2.2 #the one of the bar
                                   ,'armDensityFact': 0.083
                                   ,'sigmaArm': 0.5/2.3548 #[Wainscoat 1992]
                                   ,'armCentralCuttoff': 2.0
                                   ,'armPitch': math.radians(5.8685)
                                   ,'armInitialAngle': math.radians(80+180)
                                   }
        
        self.load_luminosity_functions()
    
    def density (self, coords):
        '''
        This function returns the arms stellar density in the given coordinates
        
        @param coords: The Coordinates object with the galactic locations where
                       the density shall be computed 
        '''
        p = self.params
        z = coords.z
        R = coords.R
        theta = coords.theta
        
        #parameters of the spiral
        b = math.tan(p.armPitch)
        a = p.armInitialRadius*math.exp(-b*p.armInitialAngle)
        a2 = p.armInitialRadius*math.exp(-b*(p.armInitialAngle+math.pi))
    
        # calculate the floating point approximation for n
        n1 = (np.log(R/a)/b - theta)/(2.0*math.pi)
        n2 = (np.log(R/a2)/b - theta)/(2.0*math.pi)
    
        # find the two possible radii for the closest point
        upper_r1 = a * np.exp(b * (theta + 2.0*math.pi*np.ceil(n1)))
        lower_r1 = a * np.exp(b * (theta + 2.0*math.pi*np.floor(n1)))
        upper_r2 = a2 * np.exp(b * (theta + 2.0*math.pi*np.ceil(n2)))
        lower_r2 = a2 * np.exp(b * (theta + 2.0*math.pi*np.floor(n2)))
    
        newR = np.where (np.abs(upper_r1 - R)<p.sigmaArm, upper_r1, 0.0)
        ind = np.abs(R - lower_r1)<p.sigmaArm
        newR[ind] = lower_r1[ind]
        ind = np.abs(upper_r2 - R)<p.sigmaArm
        newR[ind] = upper_r2[ind]
        ind = np.abs(R - lower_r2)<p.sigmaArm
        newR[ind] = lower_r2[ind]
        newR[R>17.0] = 0.0
    
        ind = newR>p.armInitialRadius
        armPointsR = newR[ind]
        armPointsZ = z[ind]
        armPointsTheta = theta[ind]
        rhoArmPoints = self.spiralArmsCentralDensity (p, armPointsR, armPointsTheta, armPointsZ)
        rhoarm = np.zeros(R.size)
        rhoarm[ind] = rhoArmPoints
        return rhoarm
    
    def spiralArmsCentralDensity (self, params, R, theta, z):

        hz = params.hz_arms_sun
        res = params.armDensityFact*params.A*(params.hz_arms_sun/hz)*np.exp((params.Rsun-R)/params.hr_thin)*np.exp(-abs(z)/hz)
        ind = R<params.Rsun
        res[ind] *= np.exp(params.hr_hole/params.Rsun)*np.exp(- params.hr_hole/R[ind])
        return res
                        
    def load_luminosity_functions (self, lf_files={}):
        if len(lf_files) == 0:
            bands = ['B', 'V', 'J', 'H', 'K', '12', '25']
            moddir = path.dirname(__file__)
            lf_files = dict([ (band, path.join(moddir, 'arms_lf_'+band+'.npy')) for band in bands ])
        self._lf = {}
        for band in bands:
            self._lf[band] = np.load(lf_files[band])
        self.set_passband(bands[0])


class Model(object):
    '''
    This is the main class of the Galiac package, which implements
    the module itself
    '''


    def __init__(self, components_list=None, params=None):
        '''
        Constructor
        
        @param components_list: The components that will build up the model. 
                                If it is None, a default list of components 
                                is used
        '''        
        
        self.params = ModelParameters()
                
        if components_list is not None:
            self.components = components_list
        else:
            self.components = self.create_default_components()
        
        if params is None:
            self.initialize_default_parameters()
        else:
            self.set_parameters (params)
            
        self.extinction_model = extinction.Drimmel(self.params.Rsun, self.params.Zsun*1e-3)
        self.set_passband('B')    
    
    def star_counts (self, l, b, w, m, dr=0.1, in_radians=False):
        '''
        Function that gives a prediction for the star counts in a given sky 
        location.
        
        This is the main method of the model. It computes the integrated star 
        counts along the line of sight, given a galactic coordinate, solid 
        angle and limiting magnitude.
        
        @param l: galactic longitude
        @param b: galactic latitude
        @param w: solid angle
        @param m: limit in apparent magnitude
        @param dr: heliocentric radius integration step, in kpc
        @param in_radians: if it is True, the l, b and w parameters are given in 
                           radians, or in degrees otherwise  
        '''
        #build the coordinates
        r = np.linspace (dr, self.params.rmax, self.params.rmax/dr)
        coords = Coordinates (l, b, r, in_radians)
        if not in_radians:
            w = w*(math.pi/180.0)**2
        
        #compute the limiting absolute magnitudes
        M = self.app_mag_to_absolute (m, coords)
        
        #compute star counts for each model component
        counts = 0.0
        for component in self.components.values():   
            if component.enabled:
                #compute component density in each integration point
                rho = component.density (coords) 
            
                #compute the integrated luminosity function
                iLF = self.integrated_luminosity_function (component, M)
                
                #compute the integrand
                dA = w*(r**2)*rho*iLF
                counts += integrate.trapz(dA, r)
        
        return counts
                         
    def set_area (self, area):
        self._w = area
    
    def set_magnitude_limit(self, mag):
        self._mag = mag
        
    def star_counts_cp(self, coords):
        l = coords[0]
        b = coords[1]
        return self.star_counts(l, b, self._w, self._mag)
            
    def app_mag_to_absolute (self, m, coords):
        '''
        Computes the absolute magnitude from the apparent magnitude,
        coordinates and extinction parameters
        '''
        if self.params.use_extinction:
            ext = self.extinction_model.compute_extinction (coords, self.current_band, self.params.extinction_params)
        else:
            ext = 0.0
        
        return m + 5 - 5*np.log10(coords.r_pc()) - ext 
    
    
    def integrated_luminosity_function(self, component, Mlim):
        ''' 
        For each absolute magnitude in the Mlim vector, it
        computes the integral of the luminosity function of
        the component, from an absolute magnitude of -infinity
        to Mlim
        '''
         
        intlf = component.cum_luminosity_function
        lf = component.luminosity_function
        ind = np.searchsorted (lf[0], Mlim)
        ind_before_interval = np.where(Mlim<=lf[0][0])
        ind_after_interval = np.where(ind >= lf[0].size)
        ind[ind_after_interval] = -1
        res = np.where (ind >=2, intlf[ind-2], 0.0)
        res += (Mlim-lf[0][ind-1])*(lf[1][ind-1]+0.5*(lf[1][ind]-lf[1][ind-1])*(Mlim-lf[0][ind-1])/(lf[0][ind]-lf[0][ind-1]))
        res[ind_before_interval] = 0.0
        res[ind_after_interval] = 1.0

        return res
    
    def set_passband (self, band):
        self.current_band = band
        for component in self.components.values():
            component.set_passband(band)
            
    def create_default_components(self):
        components = { 
                       'ThinDisk': ThinDisk(self.params)
                      ,'ThickDisk': ThickDisk(self.params)
                      ,'Bulge': Bulge(self.params)
                      ,'Bar': Bar(self.params)
                      ,'Halo': Halo(self.params)
                      ,'SpiralArms': SpiralArms(self.params)
                      }
        return components
    
    def initialize_default_parameters(self):
        self.params.update ({
                              'Rsun': 8.0 #galactocentric distance of the Sun (kpc)
                              ,'Zsun': 15.0
                             ,'A': 0.055e9 #stellar density in the Sun vicinity 
                             ,'rmax': 40.0 #maximum integration distance
                             ,'use_extinction': True
                             ,'extinction_params': {'useCobe':True, 'cobePix':None}
                             })
        for component in self.components.values():
            self.params.update(component.default_parameters)
                
    
    def set_parameters(self, new_params):
        if type(new_params) is dict:
            np = new_params
        else:
            np = vars(new_params)
        self.params.update(np)
        