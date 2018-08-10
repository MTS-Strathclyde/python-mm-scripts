import numpy as np
import matplotlib.pyplot as plt


class PointBrowser:
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """
    def __init__(self, fig, ax, ax2, line, xs, ys, X):
        self.lastind = 0

        self.fig = fig
        self.ax2 = ax2
        self.xs = xs
        self.ys = ys
        self.line = line

        self.text = ax.text(0.05, 0.95, 'selected: none',
                            transform=ax.transAxes, va='top')
        self.selected,  = ax.plot([self.xs[0]], [self.ys[0]], 'o', ms=12, alpha=0.4,
                                color='yellow', visible=False)

    def onpress(self, event):
        print event
        if self.lastind is None: return
        if event.key not in ('n', 'p'): return
        if event.key=='n': inc = 1
        else:  inc = -1


        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(self.xs)-1)
        self.update()

    def onpick(self, event):
       print event
       if event.artist!=self.line: return True

       N = len(event.ind)
       if not N: return True

       # the click locations
       x = event.mouseevent.xdata
       y = event.mouseevent.ydata


       distances = np.hypot(x-self.xs[event.ind], y-self.ys[event.ind])
       indmin = distances.argmin()
       dataind = event.ind[indmin]

       self.lastind = dataind
       self.update()

    def update(self):
        if self.lastind is None: return

        dataind = self.lastind

        self.ax2.cla()
        self.ax2.plot(self.X[dataind])

        self.ax2.text(0.05, 0.9, 'mu=%1.3f\nsigma=%1.3f'%(self.xs[dataind], self.ys[dataind]),
                 transform=self.ax2.transAxes, va='top')
        self.ax2.set_ylim(-0.5, 1.5)
        self.selected.set_visible(True)
        self.selected.set_data(self.xs[dataind], self.ys[dataind])

        self.text.set_text('selected: %d'%dataind)
        self.fig.canvas.draw()


def main():
    X = np.random.rand(100, 200)
    xs = np.mean(X, axis=1)
    ys = np.std(X, axis=1)

    fig, (ax, ax2) = plt.subplots(2, 1)
    ax.set_title('click on point to plot time series')
    line, = ax.plot(xs, ys, 'o', picker=5)  # 5 points tolerance


    browser = PointBrowser(fig, ax, ax2, line, xs, ys, X)

    fig.canvas.mpl_connect('pick_event', browser.onpick)
    print fig
    fig.canvas.mpl_connect('key_press_event', browser.onpress)

    fig.show()

main()