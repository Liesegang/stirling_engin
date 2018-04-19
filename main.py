#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import *
import matplotlib.pyplot as plt
import numpy as np

base = 273.15
cm3 = 10 ** (-6)
kpa = 10 ** 3
R = 289
rho = 1.293


#V_SE = 0.628 * cm3
#V_SC = 0.628 * cm3
#V_D = 1.5 * cm3

T_E = 600 + base
T_C = 50 + base

#m() = 1.9268 * 10 ** (-6)

alpha = pi / 2


tau = T_C / T_E
T_R = (T_E + T_C) / 2

def X():
    return V_D / V_SE

def kappa():
    return V_SC / V_SE

def S():
    return tau + 4 * tau * X() / (1 + tau) + kappa()

def B():
    return sqrt(tau ** 2 + 2 * tau * kappa() * cos(alpha) + kappa() ** 2)

def delta():
    return B() / S()

def V_E(theta):
    return (1/2) * V_SE * (1 - cos(theta))

def V_C(theta):
    return (1/2) * V_SC * (1-cos(theta - alpha))

def V(theta):
    return V_E(theta) + V_D + V_C(theta)

def phi():
    res = atan(kappa() * sin(alpha) / (tau + kappa() * cos(alpha)))
    if res >= 0:
        return res
    else:
        return atan(kappa() * sin(alpha) / (tau + kappa() * cos(alpha)) + pi)

def m():
    return rho * V(pi/2)

def P(theta):
    return 2 * m() * R * T_C / (V_SE * (S() - B() * cos(theta - phi())))

def P_mean():
    return 2 * m() * R * T_C / (V_SE * sqrt(S()**2 - B()**2))

def m_from_pmean(pmean):
    return V_SE * sqrt(S() ** 2 - B() ** 2) * pmean / (2 * R * T_C)

def W_E():
    return P_mean() * V_SE * pi * delta() * sin(phi()) / (1 + sqrt(1 - delta() ** 2))

def W_C():
    return -1 * P_mean() * V_SE * pi * delta() * tau * sin(phi()) / (1 + sqrt(1 - delta() ** 2))

def W_i():
    print(V_SE)
    return W_E() + W_C()

def plot():
    x  = np.arange(0, 2 * pi, 0.1)
    p  = list(map(lambda x : P(x) / kpa, x))
    v  = list(map(lambda x : V(x) / cm3, x))
    ve = list(map(lambda x : V_E(x) / cm3, x))
    vc = list(map(lambda x : V_C(x) / cm3, x))

    plt.title("P-V線図")
    plt.xlabel("容積 V(cm3)")
    plt.ylabel("圧力 P(kPa)")
    plt.plot(v,p, "-o")
    plt.show()

def print_detail():
    print(f"高温部行程容積(V_SE):\t{V_SE/cm3} cm3")
    print(f"低温部行程容積(V_SC):\t{V_SC/cm3} cm3")
    print(f"全無効容積(V_D):\t{V_D/cm3} cm3")
    print()

    print(f"高温部温度(T_E):\t{T_E - base} ℃")
    print(f"低温部温度(T_C):\t{T_C - base} ℃")
    print()

    print(f"エンジン内ガス質量(m):\t{m() * 10 ** 6} mg")
    print()

    print(f"位相角(alpha):\t\t{alpha * 360 / 2 / pi} degree")
    print()

    print(f"理論図示仕事(W_i):\t{W_i()} W")

def plot_kappa():
    W = []
    k = []

    V_all = 1.5 * cm3
    global V_D
    V_D = 1.5 * cm3

    for t in np.arange(0.3, 1, 0.0001):
        global V_SE
        global V_SC

        V_SE = V_all * t
        V_SC = V_all * (1-t)

        k.append(kappa())
        W.append(W_i())

    plt.plot(k, W)
    plt.show()

def plot_v():
    W = []
    V = []

    for v in np.arange(1, 100, 0.1):
        v *= cm3

        global V_D
        global V_SE
        global V_SC

        V_D = v
        V_SE = v * 0.5
        V_SC = v * 0.5 

        V.append(v)
        W.append(W_i())

    plt.plot(V, W)
    plt.show()


if __name__ == '__main__':
    plot_v()
