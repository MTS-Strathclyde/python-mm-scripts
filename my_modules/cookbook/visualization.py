# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 07:16:00 2015

@author: max
"""

import matplotlib.pyplot as plt
import numpy as np


class PosterFig(object):
    def __init__(self, fig_width=7.824, fig_height=5.956,
                 alpha=.6, var_name=r'\Delta G', method_name='HNC/ISc',
                 units='kcal/mol', labelsize=35):
        """
        New poster fig class.
        
        Good colours:
        [i/255. for i in (63,74,77)],
        [i/255. for i in (95,158,209)],
        [i/255. for i in (105,89,189)]
        
        Parameters
        ----------
        x_data : list
            List of list-like arrays of points.
        
        y_data : list
            List of list-like arrays of points.
        """
        # general
        self.back_color = [i/255. for i in (249,250,255)]
        # box positions
        self.box_pos_x = .47
        self.box_pos_y = .62
        # text parameters
        xlabel = '$\mathrm{' + var_name + '_{exp}}$'
        ylabel = '$\mathrm{' + var_name + '_{' + method_name + '}}$'
        xlabel = xlabel + units
        ylabel = ylabel + units
        # initalize figure
        self.fig, self.ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
        self.ax.set_ylabel(ylabel, fontsize=labelsize)
        self.ax.set_xlabel(xlabel, fontsize=labelsize)
        

    def plot_scatter(self, x_data, y_data, colour=(0.247,0.29,0.302), 
                     marker='o', size=35, alpha=.6, label=None):
#        if len(self.x_data) == 1:
#            self.textstr = "RMSE={}\nSD={}\nbias={}".format(*self.get_errors())
        self.ax.scatter(x_data, y_data, c=colour, marker=marker, s=size, 
                        alpha=alpha, label=label)

    def draw_diag_line(self, alpha=.6, colour='b'):
        xmin, xmax = self.ax.xlim()
        ymin, ymax = self.ax.ylim()
        lmin = min(xmin, ymin)
        lmax = max(xmax, ymax)
        self.ax.plot([lmin, lmax], [lmin, lmax], '-', alpha=alpha, c=colour)
    
    def add_error_box(self, x_data, y_data, x_pos=.47, y_pos=.62,
                      fontsize=35):
        def rmse(errs): return np.sqrt(np.average(errs**2))
        def sepb(errs): return np.std(errs, ddof=1)
        def bias(errs): return np.average(errs)        
        mod_err = y_data - x_data
        textstr = "RMSE={:.2f}\nSD={:.2f}\nbias={:.2f}".format(rmse(mod_err),
                                                 sepb(mod_err), bias(mod_err))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in lower left in axes coords
        self.ax.text(x_pos, y_pos, textstr, transform=self.ax.transAxes, 
                     fontsize=fontsize, ha='right', va='bottom', bbox=props)
        
    def set_limits(self,limits=None, expand=.01):
        if limits:
            xmin, xmax, ymin, ymax = limits
        else:
            xmin, xmax = self.ax.xlim()
            ymin, ymax = self.ax.ylim()
            lmin = min(xmin, ymin)
            lmax = max(xmax, ymax)
            ax_border = (lmax - lmin)*expand
            xmin -= ax_border
            xmax += ax_border
            ymin -= ax_border
            ymax += ax_border
        self.ax.xlim(xmin, xmax)
        self.ax.ylim(xmin, xmax)
    
            
fig = PosterFig()
fig.plot_scatter([1,2,3],[.1, .3, .1])
plt.show()
    
