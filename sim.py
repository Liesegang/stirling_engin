#!/usr/bin/env python
# -*- coding: utf-8 -*-

import stirling
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

cm3 = 10 ** -9

if __name__ == '__main__':
    # write your code here
    x = np.linspace(0, 10, 101)[1:-1]
    y = np.linspace(0, 180, 1000)
    X, Y = np.meshgrid(x, y)

    s = stirling.Stirling_Engin(1 * cm3, 1 * cm3, X * cm3, 400, 50, Y)
    Z = s.W_i()

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_wireframe(X,Y,Z)
    ax.set_xlabel("体積比")
    ax.set_ylabel("無効体積(0-1)")
    ax.set_zlabel("仕事(W_i)")
    plt.show()
