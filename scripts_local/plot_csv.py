#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 09:13:41 2016

@author: max
"""

import bokeh.io as bokio
import bokeh.plotting as bokpl
import bokeh.models as bokml

import pandas as pd
import sys
import numpy as np

def rmse(x,y):
    return np.sqrt(np.nanmean((y-x)**2))

def bias(x,y):
    return np.nanmean(y-x)

def sepb(x,y):
    return np.nanstd(y-x, ddof=0)

def mape(x,y):
    return np.mean(np.abs((x - y) / x)) * 100


def create_plot(x, y, df, source=None, title=None, plot_size=500, circ_size=8):
    if not title:
        title=y
    if not source:
        source = bokpl.ColumnDataSource(df)
    TOOLS = ['resize,pan,wheel_zoom,box_zoom,reset']
    plot_options = dict(width=plot_size, plot_height=plot_size)
    figure = bokpl.figure(tools=TOOLS + \
                        [bokml.HoverTool(tooltips=[("(exp,pred)", "(@{}, @{})".format(x,y)),
                                               ("name", "@Name")])],
                        **plot_options)
    figure.title = title
    limits = (np.amin(np.r_[df[x], df[y]]), np.amax(np.r_[df[x], df[y]]))
    figure.line(limits, limits, line_color='black', line_width=2)
    figure.circle(x, y, source=source, size=circ_size)
    plot_range = limits[1] - limits[0]
    box_l = limits[0] + plot_range*0.55
    box_b = limits[0] + plot_range*0.25
    textstrs = ["RMSE={:.3f}".format(rmse(df[x], df[y])),
                'SD={:.3f}'.format(sepb(df[x], df[y])),
                'bias={:.3f}'.format( bias(df[x], df[y]))]
    for t in textstrs:
        box_b = box_b - plot_range*0.07
        figure.text(x=box_l, y=box_b, text_align="left", text=[t],
                    text_font_size="12pt")
    return figure

def main(argv):
    bokio.output_file('plot.html')
    csv = argv[0]
    df = pd.read_csv(csv)
    df = df.dropna()
    col1 = argv[1]
    col2 = argv[2]
    p = create_plot(col1, col2, df)
    bokio.show(p)
    

if __name__ == '__main__':
    main(sys.argv[1:])    
    

    