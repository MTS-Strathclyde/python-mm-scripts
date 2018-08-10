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
avgd_points = open('averaged_points.txt', 'a')
interp_points = open('interpoltd_points.txt', 'a')


class PointBrowser:
    """
    Delete stuff
    Press d to delete point
    Press x to delete molecule (deprecated)
    Close figure to continue
    """
    def __init__(self, fig, ax, int_t, int_v, avg_t, avg_v,
                 line, title, true_t, true_v):
        self.fig = fig
        self.line = line
        self.ax = ax
        self.int_t = np.array(int_t)
        self.int_v = np.array(int_v)
        self.avg_t = avg_t
        self.avg_v = avg_v
        self.true_t = true_t
        self.true_v = true_v
        self.title = title
        
        X = create_fit_m(true_t)
        e_fit_ob = sm.OLS(true_v, X).fit()
        self.r2 = e_fit_ob.rsquared

        
        self.lastind = 0
        self.text = ax.text(0.05, 0.95, 'selected: none\nR2:{}'.format(self.r2),
                            transform=ax.transAxes, va='top')
        self.selected, = ax.plot(int_t[self.lastind], int_v[self.lastind], 
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
        
        distances = np.hypot(x-self.int_t[event.ind], y-self.int_v[event.ind])
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
#        elif event.key == 'x':
#            self.del_molecule = not self.del_molecule
#            self.text.set_text("Molecule will be deleted={}".format(self.del_molecule))
#            self.fig.canvas.draw()
        else:
            return
    
    
    def delete(self, ind):
        self.to_del.append((self.int_t[ind], self.int_v[ind]))
        
        self.ax.cla()
        self.ax.set_title(self.title)
        self.int_t = np.delete(self.int_t, ind)
        self.int_v = np.delete(self.int_v, ind)

        dum_range = np.linspace(273, 373, 200)
        X = create_fit_m(self.true_t)
        e_fit_ob = sm.OLS(self.true_v, X).fit()
        self.ax.plot(self.true_t, self.true_v, 'x')
        self.ax.plot(self.avg_t, self.avg_v, 'rx')        
        self.ax.plot(dum_range, e_fit_ob.predict(create_fit_m(dum_range)))

        self.line, = self.ax.plot(self.int_t, self.int_v, 'ro', picker=5)  # 5 points tolerance
        
        

    def update(self):
        if self.lastind is None: return

        dataind = self.lastind

        self.selected.set_visible(True)
        self.selected.set_data(self.int_t[dataind], self.int_v[dataind])
        
        t = self.int_t[dataind]
        v = self.int_v[dataind]

        self.text.set_text("""T: {}\nG"""\
                            .format(round(t,1), round(v,2)))
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



    
def eval_params(params, t):
    dt = t - 298.15
    lt = dt - t*np.log(t/298.15)
    
    return params[0] + dt*params[1] + lt*params[2]


def interpolate(mol):
    print 'mol'
    #new_temps = [0, 10, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
    new_temps = [25]
    mol.all_exp_temps = sorted(list(mol.t_en_e[:, 0]))
    mol.new_temp_dic = {}
    for t in new_temps:
        mol.new_temp_dic[t + 273.0] = []
    for exp in mol.experiments:
        if round(exp.T - 273.0) in new_temps:
            mol.new_temp_dic[round(exp.T)].append((exp.T, exp))
    mol.averaged = []
    mol.interpolated = []
    for k, v in mol.new_temp_dic.iteritems():
        if len(v) == 1:
            mol.averaged.append([k, v[0][1].value, 'AVG', v[0][1].source])
        elif len(v) > 1:
            avg_value = np.average(np.array([i[1].value for i in v], dtype=float))
            total_source = ' _AND_ '.join([i[1].source for i in v])
            mol.averaged.append([k, avg_value, 'AVG', total_source])            
        elif len(v) == 0:
            point_above = False
            point_belove = False
            for p in mol.all_exp_temps:
                if 10 >= k-p > 0:   #modified for 25 interpolation
                    point_belove = True
                elif 0 > k-p >= -10:
                    point_above = True
            if point_above and point_belove:
                interp_v = eval_params(mol.e_fit, k)
                mol.interpolated.append([k, interp_v, 'INTERP', ''])
                
    if mol.interpolated:
        avg_t = np.array([i[0] for i in mol.averaged])
        avg_v = np.array([i[1] for i in mol.averaged])
        
        int_t = np.array([i[0] for i in mol.interpolated])
        int_v = np.array([i[1] for i in mol.interpolated])
    
        true_t = mol.t_en_e[:, 0]
        true_v = mol.t_en_e[:, 1]

        fig, ax = plt.subplots()
        title = mol.db_mol.IUPACName
        ax.set_title(title)
        
        
        ax.plot(mol.t_en_e[:, 0], mol.t_en_e[:, 1], 'x')
        ax.plot(avg_t, avg_v, 'rx')        
        dum_range = np.linspace(273, 373, 200)
        X = create_fit_m(mol.t_en_e[:, 0])
        e_fit_ob = sm.OLS(mol.t_en_e[:, 1], X).fit()
        ax.plot(dum_range, e_fit_ob.predict(create_fit_m(dum_range)))
        
        
        line, = ax.plot(int_t, int_v, 'ro', picker=5)  # 5 points tolerance
        
        browser = PointBrowser(fig, ax, int_t, int_v, avg_t, avg_v, 
                               line, title, true_t, true_v)
    
        fig.canvas.mpl_connect('pick_event', browser.onpick)
        fig.canvas.mpl_connect('key_press_event', browser.onpress)
    
        plt.show()
        for t, v in browser.to_del:
            for i, entry in enumerate(mol.interpolated):
                if len(entry) == 4:
                    T, val, _, _ = entry
                    if abs(T - t) < 1e-6 and abs(val - v) < 1e-6:
                        mol.interpolated[i] = 'd'
                        break
        
    added = 0
    interpolated_count = 0
    for t, v, method, source in mol.averaged:
        db_exp = sd.Experiment(Temperature=t, SolvEnergy=v,
                               Method=method, Source=source)
        mol.db_mol.Experiments.append(db_exp)
        added += 1
        descrip = '{} {} {} {}\n'.format(t, v, method, source)
        avgd_points.write(descrip)

    for entry in mol.interpolated:
        if len(entry) == 4:
            t, v, method, source = entry
            db_exp = sd.Experiment(Temperature=t, SolvEnergy=v,
                                   Method=method, Source=source)
            mol.db_mol.Experiments.append(db_exp)
            added += 1
            interpolated_count += 1
            descrip = '{}, {}, {}, {}, {}\n'.format(mol.db_mol.ID, t, v, method, source)
            interp_points.write(descrip)
            interp_points.flush()
            
    print added,'points added'
    print interpolated_count,'of them interpolated'

    ses.flush()
    ses.commit()
        
    

def main():    

    molecules = ses.query(sd.Molecule).filter(sd.Molecule.Experiments != None).all()

    mol_objs = []
    for mol in molecules:
        mol_objs.append(Molecule(mol))    
    
    for mol in mol_objs:
        interpolate(mol)
        
        

try:
    main()
finally:
    avgd_points.close()
    interp_points.close()
    ses.close()
    