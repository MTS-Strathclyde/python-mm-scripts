# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 13:17:48 2014

@author: max
"""

import numpy as np
import statsmodels.api as sm

import solvation_database as sd
from db_interface import create_session
from sqlalchemy import and_
import matplotlib.pyplot as plt


ses = create_session()
del_points_f = open('deleted_points.txt', 'a')
del_mol_f = open('deleted_molecules.txt', 'a')


class PointBrowser:
    """
    Delete stuff    
    """
    def __init__(self, fig, ax, temps, values, line, title):
        self.fig = fig
        self.line = line
        self.ax = ax
        self.temps = temps
        self.values = values
        self.title = title
        
        X = create_fit_m(temps)
        e_fit_ob = sm.OLS(values, X).fit()
        self.r2 = e_fit_ob.rsquared

        
        self.lastind = 0
        self.text = ax.text(0.05, 0.95, 'selected: none\nR2:{}'.format(self.r2),
                            transform=ax.transAxes, va='top')
        self.selected, = ax.plot(temps[self.lastind], values[self.lastind], 
                                  'o', ms=12, alpha=0.4,
                                  color='yellow', visible=False)
                                  
        self.to_del = []
        self.del_molecule = False
        
        #        print 'fig',id(self.fig)
        #        print 'ax',id(ax)
     
                                  
    def onpick(self, event):
       if event.artist!=self.line: return True
    
       N = len(event.ind)
       if not N: return True
    
       # the click locations
       x = event.mouseevent.xdata
       y = event.mouseevent.ydata
    
    
       distances = np.hypot(x-self.temps[event.ind], y-self.values[event.ind])
       indmin = distances.argmin()
       dataind = event.ind[indmin]
    
       self.lastind = dataind
       self.update()

    def onpress(self, event):
        if self.lastind is None: return
        if event.key == 'd':
            self.delete(self.lastind)
            self.lastind = 0
            self.update()
        elif event.key == 'x':
            self.del_molecule = not self.del_molecule
            self.text.set_text("Molecule will be deleted={}".format(self.del_molecule))
            self.fig.canvas.draw()
        else:
            return
    
    
    def delete(self, ind):
        self.to_del.append((self.temps[ind], self.values[ind]))
        
        self.ax.cla()
        self.ax.set_title(self.title)
        self.temps = np.delete(self.temps, ind)
        self.values = np.delete(self.values, ind)

        dum_range = np.linspace(273, 373, 200)
        X = create_fit_m(self.temps)
        e_fit_ob = sm.OLS(self.values, X).fit()
        self.r2 = e_fit_ob.rsquared
        self.ax.plot(dum_range, e_fit_ob.predict(create_fit_m(dum_range)))

        self.line, = self.ax.plot(self.temps, self.values, 'o', picker=5)  # 5 points tolerance


        self.text = self.ax.text(0.05, 0.95, 'selected: none\nR2:{}'.format(self.r2),
                            transform=self.ax.transAxes, va='top')
        
        

    def update(self):
        if self.lastind is None: return

        dataind = self.lastind

        self.selected.set_visible(True)
        self.selected.set_data(self.temps[dataind], self.values[dataind])
        
        t = self.temps[dataind]
        v = self.values[dataind]

        potential_temps = np.delete(self.temps, dataind)
        potential_values = np.delete(self.values, dataind)
        X = create_fit_m(potential_temps)
        e_fit_ob = sm.OLS(potential_values, X).fit()
        
        resid = e_fit_ob.resid[dataind]
        r2 = e_fit_ob.rsquared

        self.text.set_text("""T: {}\nG: {}\nDeviation: {}
On delete R2 change: {}"""\
                            .format(round(t,1), round(v,2), round(resid,3),
                                    round(r2 - self.r2, 4)))
        self.fig.canvas.draw()


class Molecule(object):
    def __init__(self, db_mol):
        self.db_mol = db_mol
        self.deleted = False
        self.experiments = []
        self.e_fit = None
        self.t_en_e = None
        self.get_experiments()
        self.get_data()
        
    def get_data(self):
        e_fit_ob = build_experimental_ols(self.db_mol)
        self.e_fit = e_fit_ob.params
        self.fit_er = e_fit_ob.bse
        self.t_en_e = create_exp_t_energy_m(self.db_mol)
        
    def get_experiments(self):
        for exp in self.db_mol.Experiments:
            self.experiments.append(Experiment(self, exp))

    def delete(self):
        ses.delete(self.db_mol)
        self.deleted = True
    
    def __repr__(self):
        return self.db_mol.__repr__()

        
class Experiment(object):
    def __init__(self, mol_ob, db_exp):
        self.mol = mol_ob
        self.db_exp = db_exp
        self.deleted = False
        
        self.T = self.db_exp.Temperature
        self.value = self.db_exp.SolvEnergy
        self.method = self.db_exp.Method
        self.source = self.db_exp.Source
        
    def delete(self):
        ses.delete(self.db_exp)
        self.deleted = True


def create_exp_t_energy_m(molecule):
    """Create matrix (T, Exp_E)"""
    t_energy_m = []
    for exp in molecule.Experiments:
        t_energy_m.append((float(exp.Temperature),float(exp.SolvEnergy)))
    t_energy_m.sort()
    return np.array(t_energy_m)


def create_fit_m(t):
    """ Given t create (ones, delta T, deltaT - T*ln T/T0) matrix """
    dt = t - 298.15
    lt = dt - t*np.log(t/298.15)
    return np.vstack((np.ones(t.shape), dt, lt)).T


def build_experimental_ols(molecule):
    """Create statmodel result object describing fit of experimental data to curve."""
    t_energy = create_exp_t_energy_m(molecule)
    X = create_fit_m(t_energy[:, 0])
    return sm.OLS(t_energy[:, 1], X).fit()


def investigate_mol(mol):

    temps = mol.t_en_e[:, 0]
    values = mol.t_en_e[:, 1]

    fig, ax = plt.subplots()
    title = mol.db_mol.IUPACName
    ax.set_title(title)
    
    dum_range = np.linspace(273, 373, 200)
    X = create_fit_m(mol.t_en_e[:, 0])
    e_fit_ob = sm.OLS(mol.t_en_e[:, 1], X).fit()
    ax.plot(dum_range, e_fit_ob.predict(create_fit_m(dum_range)))
    
    
    line, = ax.plot(temps, values, 'o', picker=5)  # 5 points tolerance
    
    browser = PointBrowser(fig, ax, temps, values, line, title)

    fig.canvas.mpl_connect('pick_event', browser.onpick)
    fig.canvas.mpl_connect('key_press_event', browser.onpress)

    plt.show()
    for t, v in browser.to_del:
        for exp in mol.experiments:
            if abs(exp.T - t) < 1e-6 and abs(exp.value - v) < 1e-6:
                exp.delete()
                descrip = '{} {} {} {}\n'.format(exp.T, exp.value, exp.method, exp.source)
                #print descrip
                del_points_f.write(descrip)
                break
        continue
    del_points_f.flush()
    print len(browser.to_del),'points removed'
    if browser.del_molecule:
        for exp in mol.experiments:
            exp.delete()
        ses.commit()
        mol.delete()
        del_mol_f.write('{}\n'.format(mol.db_mol.IUPACName))
        del_mol_f.flush()
        print mol,'deleted'
    ses.flush()
    ses.commit()
    

def main():    

    molecules = ses.query(sd.Molecule).filter(sd.Molecule.Experiments != None).all()

    mol_objs = []
    for mol in molecules:
        mol_objs.append(Molecule(mol))    
    
    for mol in mol_objs:
        investigate_mol(mol)
        
        

try:
    main()
finally:
    del_points_f.close()
    del_mol_f.close()
    ses.close()
    