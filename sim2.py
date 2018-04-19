#!/usr/bin/env python
# -*- coding: utf-8 -*-

import stirling
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

cm3 = 10 ** -6

if __name__ == '__main__':
    # write your code here
    x = np.linspace(0, 1, 1000)[1:-1]

    s = stirling.Stirling_Engin(40 * x * cm3, 40 * (1-x) * cm3, 10 * cm3, 600, 100, 90, 1000)
    plt.plot(x, s.W_i())
    plt.show()
