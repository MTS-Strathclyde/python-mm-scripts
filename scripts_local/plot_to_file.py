#!/usr/bin/env python2
import numpy as np
import matplotlib.pyplot as plt
import sys

mat = np.loadtxt(sys.argv[1])

for i in range(1, mat.shape[-1]):
    plt.plot(mat[:,0], mat[:,i])
    plt.savefig(sys.argv[1] + '.png', dpi=300)

