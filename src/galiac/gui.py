'''


@author: jmolgo
'''
import Tkinter as tk
import ttk
import tkFileDialog

import numpy as np
from matplotlib import pyplot as plt

import galiac.core as gc 
import galiac.parallel as gpar 
import galiac.catalog as gcat



class GaliacGUI(ttk.Frame): 
    '''
    This class implements a grafical user interface for
    the Galiac model, using python TKinter
    '''

    def __init__(self,  isapp=True, name='galiac'):        

        self.model = gc.Model()
        self.modcalc = gpar.ModCalc_Multiprocessing(model=self.model)
        self.catalog = gcat.TwoMASS()
        
        #st = ttk.Style()
        #st.theme_use ('alt')
        
        ttk.Frame.__init__(self, name=name)
        self.pack(expand=tk.Y, fill=tk.BOTH)
        self.master.title(name)
        self.isapp = isapp
        
        self.top_panel = ttk.Frame(self, name=name)
        self.top_panel.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.Y)
        self.create_general_panel()
        self.nb = ttk.Notebook(self.top_panel, name=name)
        self.nb.pack(fill=tk.BOTH, expand=tk.Y, padx=2, pady=3)
        
        self.panels = [ ComponentPanel(self.nb, self.model, 'Thin disk', self.model.components['ThinDisk'])
                       ,ComponentPanel(self.nb, self.model, 'Thick disk', self.model.components['ThickDisk'])
                       ,ComponentPanel(self.nb, self.model, 'Bulge', self.model.components['Bulge'])
                       ,ComponentPanel(self.nb, self.model, 'Bar', self.model.components['Bar'])
                       ,ComponentPanel(self.nb, self.model, 'Spiral Arms', self.model.components['SpiralArms'])
                       ,ComponentPanel(self.nb, self.model, 'Halo', self.model.components['Halo'])]
        
    def create_general_panel(self):
        #create the panel itself
        panel = ttk.Frame(self.top_panel, name='general')
        self.general_panel = panel
        
        #create the variables
        self.init_l_var = tk.DoubleVar (value=-180.0)
        self.end_l_var = tk.DoubleVar (value=180.0)
        self.b_var = tk.DoubleVar (value=40)
        self.area_var = tk.DoubleVar (value=1.0)
        self.mag_var = tk.DoubleVar(value=13.0)
        
        #create the widgets
        row=0
        label = ttk.Label (panel, text='l (from)')
        entry = ttk.Entry (panel, textvariable=self.init_l_var)
        label.grid (row=row, column=0, sticky=tk.W+tk.N)
        entry.grid (row=row, column=1, sticky=tk.W+tk.N)
        row += 1
        
        label = ttk.Label (panel, text='l (to)')
        entry = ttk.Entry (panel, textvariable=self.end_l_var)
        label.grid (row=row, column=0, sticky=tk.W+tk.N)
        entry.grid (row=row, column=1, sticky=tk.W+tk.N)
        row += 1

        label = ttk.Label (panel, text='b')
        entry = ttk.Entry (panel, textvariable=self.b_var)
        label.grid (row=row, column=0, sticky=tk.W+tk.N)
        entry.grid (row=row, column=1, sticky=tk.W+tk.N)
        row += 1

        label = ttk.Label (panel, text='area')
        entry = ttk.Entry (panel, textvariable=self.area_var)
        label.grid (row=row, column=0, sticky=tk.W+tk.N)
        entry.grid (row=row, column=1, sticky=tk.W+tk.N)
        row += 1

        label = ttk.Label (panel, text='magnitude')
        entry = ttk.Entry (panel, textvariable=self.mag_var)
        label.grid (row=row, column=0, sticky=tk.W+tk.N)
        entry.grid (row=row, column=1, sticky=tk.W+tk.N)
        row += 1
        
        ttk.Button (panel, text='Plot!', command = self.plot).grid(row=row)
        row += 1
        ttk.Button (panel, text='Load parameters...', command = self.load_params).grid(row=row)
        row += 1
        ttk.Button (panel, text='Save parameters...', command = self.save_params).grid(row=row)

        # add tab to notebook
        #self.nb.add(panel, text='General', underline=0, padding=20)
        panel.pack (side=tk.LEFT)
    
    def plot(self):
        self.set_model_parameters()
        linit = self.init_l_var.get()
        lend =  self.end_l_var.get()
        l = np.linspace (linit, lend, lend-linit+1)
        l = np.where (l==0, 1e-10, l)
        b = np.repeat(self.b_var.get(), l.size)
        b = np.where (b==0, 1e-10, b)
        mag = self.mag_var.get()
        area = self.area_var.get()
        
        cat_counts = self.catalog.select(l=(linit, lend), b=b[0], mag=mag)        
        cts_mod = self.modcalc.compute_region(l, b, area, mag, 'K')
        
        plt.clf()
        plt.plot (cat_counts[:,0], cat_counts[:,2], 'g.')            
        plt.plot(l, cts_mod, 'r')
        plt.grid()
        plt.minorticks_on()
        plt.gca().invert_xaxis()
        plt.title ('b='+str(b[0]))
        plt.xlabel ('l (deg)')
        plt.ylabel ('$\\rm{estrellas}/\\rm{grado}^2$')
        plt.show()
    
    def load_params(self):
        fname = tkFileDialog.askopenfilename(parent=self,title='Parameters file')
        if fname != None:
            self.model.params.from_file(fname)
            self.get_model_parameters()
    
    def save_params(self):
        fname = tkFileDialog.asksaveasfilename(parent=self,title='Save model parameters...')
        if fname != None:
            self.set_model_parameters()
            self.model.params.to_file(fname)
        
    def set_model_parameters(self):
        params = {}
        for panel in self.panels:
            for prop in panel.prop_editors:                
                if prop[2].__class__ is tk.StringVar:
                    arritems = prop[2].get().replace("[", "").replace("]","")
                    params[prop[0]]=np.fromstring(arritems, 
                                                  sep=' ', 
                                                  dtype=float)
                else:
                    params[prop[0]]=prop[2].get()
        print params
        self.model.set_parameters(params)
    
    def get_model_parameters(self):
        params = vars(self.model.params)
        for panel in self.panels:
            panel.display_model_parameters(params)
            
    
    
class ComponentPanel(object):
    '''
    Subpanel to display and setup the properties
    of a given model component
    '''
    
    def __init__(self, nb, model, name, component):
        self.root = nb
        self.model = model
        self.component = component
        
        self.frame = ttk.Frame(nb, name=name.lower())
        
        #create dynamically the editors for the 
        #component properties
        self.create_property_editors()
        
        # add tab to notebook
        nb.add(self.frame, text=name, underline=0, padding=20)
    
    def create_property_editors(self):
        self.prop_editors = []
        row = 0
        properties = self.component.param_names
        params = self.component.params        
        for prop in properties:
            if type(vars(params)[prop]) is np.ndarray:
                variable = tk.StringVar(value=vars(params)[prop])
            else:
                variable = tk.DoubleVar(value=vars(params)[prop])
            entry = ttk.Entry(self.frame, textvariable=variable)
            label = ttk.Label (self.frame, text=prop+": ")
            label.grid (row=row, column=0, sticky=tk.W+tk.N)
            entry.grid (row=row, column=1, sticky=tk.W+tk.N)
            self.prop_editors.append((prop, label, variable, entry))
            row += 1
                        
    def display_model_parameters(self, params):
        for prop in self.prop_editors:
            propname = prop[0]
            propvar = prop[2]
            if params.has_key(propname):
                propvar.set (params[propname])

if __name__ == '__main__':
    gui = GaliacGUI()
    gui.mainloop()
    
