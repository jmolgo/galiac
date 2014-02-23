'''
Handlers to use several catalogs with Galiac

@author: jmolgo
'''

import numpy as np
from os import path

class Catalog(object):
    '''
    Base class for the catalogs
    '''
    
    catalog = None
    band = None
    

class TwoMASS(Catalog):
    '''
    Handler for the 2MASS catalog
    '''

    _lower_mag = 7
    _upper_mag = 15
    _ind_mag = 2
    band = 'K'

    def __init__(self, nfile = path.join (path.dirname(__file__),'twomass.npy')):
        '''
        Constructor
        '''
        self.catalog = np.load (nfile)        
        
    def select (self, l=None, b=None, mag=14): 
        '''
        Method to perform a selection on the catalog
        '''
        
        ind_mag = self.colmag (mag)
        ind = np.repeat (True, self.catalog.shape[0])
        
        #select using longitude condition, if it is present
        if l is not None:
            #check if it is an interval
            try:
                linit = l[0]
                lend = l[1]
                ind &= (self.catalog[:,0] >= linit) & (self.catalog[:,0] <= lend)
            except (TypeError, IndexError):   
                #it is not an interval, but a specific longitude
                ind &= self.catalog[:,0] == l
                
        #select using latitude condition, if it is present
        if b is not None:
            #check if it is an interval
            try:
                binit = b[0]
                bend = b[1]
                ind &= (self.catalog[:,1] >= binit) & (self.catalog[:,1] <= bend)
            except (TypeError, IndexError):   
                #it is not an interval, but a specific longitude
                ind &= self.catalog[:,1] == b
        
        #return the counts selected sky points of the given magnitude
        filtered = self.catalog[ind]
        return filtered[:, [0, 1, ind_mag]]
    
    def colmag (self, mag):
        return mag-self._lower_mag+self._ind_mag
    
