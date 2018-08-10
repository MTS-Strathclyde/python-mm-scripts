# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 15:52:18 2014

@author: max
"""

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-1, 2), ylim=(-1, 2))
polygon = plt.Polygon([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
ax.add_patch(polygon)

# set the picker to True, so that pick events are registered
polygon.set_picker(True)

# create a function to be bound to pick events: here the event has an
# attribute `artist` which points to the object which was clicked
class A(object):
    def on_pick(self, event):
        event.artist.set_facecolor(np.random.random(3))
        fig.canvas.draw()

a = A()
# bind pick events to our on_pick function
fig.canvas.mpl_connect('pick_event', a.on_pick)
plt.show()