'''


@author: jmolgo
'''
import Tkinter as tk
import ttk


import galiac.core as gc
import galiac.parallel as gpar



class GaliacGUI(ttk.Frame):
    '''
    This class implements a grafical user interface for
    the Galiac model, using python TKinter
    '''

    def __init__(self,  isapp=True, name='galiac'):        

        self.model = gc.Model()
        self.modcalc = gpar.ModCalc_Multiprocessing(model=self.model)
        
        ttk.Frame.__init__(self, name=name)
        self.pack(expand=tk.Y, fill=tk.BOTH)
        self.master.title(name)
        self.isapp = isapp
        
        self.top_panel = ttk.Frame(self, name=name)
        self.top_panel.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.Y)
        self.nb = ttk.Notebook(self.top_panel, name=name)
        self.nb.pack(fill=tk.BOTH, expand=tk.Y, padx=2, pady=3)
        
        self.create_general_panel()
        
        
        self.panels = [ ComponentPanel(self.nb, self.model, 'Thin disk', self.model.components['ThinDisk'])
                       ,ComponentPanel(self.nb, self.model, 'Thick disk', self.model.components['ThickDisk'])
                       ,ComponentPanel(self.nb, self.model, 'Bulge', self.model.components['Bulge'])
                       ,ComponentPanel(self.nb, self.model, 'Bar', self.model.components['Bar'])
                       ,ComponentPanel(self.nb, self.model, 'Spiral Arms', self.model.components['SpiralArms'])
                       ,ComponentPanel(self.nb, self.model, 'Halo', self.model.components['Halo'])]
        
    def create_general_panel(self):
        #create the panel itself
        panel = ttk.Frame(self.nb, name='general')
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

        # add tab to notebook
        self.nb.add(panel, text='General', underline=0, padding=20)
        
    
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
            variable = tk.DoubleVar(value=vars(params)[prop])
            entry = ttk.Entry(self.frame, textvariable=variable)
            label = ttk.Label (self.frame, text=prop+": ")
            label.grid (row=row, column=0, sticky=tk.W+tk.N)
            entry.grid (row=row, column=1, sticky=tk.W+tk.N)
            self.prop_editors.append((label, variable, entry))
            row += 1            

if __name__ == '__main__':
    gui = GaliacGUI()
    gui.mainloop()
    
