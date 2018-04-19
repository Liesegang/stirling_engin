#!/usr/bin/env python
# -*- coding: utf-8 -*-

import stirling
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

cm3 = 10 ** -9

if __name__ == '__main__':
    # write your code here
    x = np.linspace(1, 360, 360)

    s = stirling.Stirling_Engin(19.847 * cm3, 21.991 * cm3, 3.848 * cm3, 600, 100, x, 1000)
    plt.plot(x, s.W_i())
    plt.show()
