'''
This module contains the implementation
of several methods to compute the Galactic Model
using different kinds of parallelism  

@author: jmolgo
'''

import galiac.core as gc
import multiprocessing
from multiprocessing import Pool
import copy_reg
import types
import numpy as np
from functools import partial

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    cls_name = ''
    if func_name.startswith('__') and not func_name.endswith('__'):
        cls_name = cls.__name__.lstrip('_')
    if cls_name:
        func_name = '_' + cls_name + func_name
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

    
class ModCalc_Multiprocessing(object):
    def __init__(self, model=None, nprocs=multiprocessing.cpu_count()):
        
        if model is None:
            self.model = gc.Model()
        else:
            self.model = model
        self.pool = Pool(processes=nprocs)
 
    
    
    def compute_region (self, l, b, area, mag, band):
        coords=np.empty((l.size, 2))
        coords[:, 0] = l
        coords[:, 1] = b  

        self.model.set_passband(band)
        self.model.set_area (area)
        self.model.set_magnitude_limit(mag)
        f = partial (gc.Model.star_counts_cp, self.model)
        return self.pool.map (f, coords)
