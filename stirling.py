#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
import numpy



R = 289
Rho = 1.293

cm3 = 10 ** -6
kPa = 10 ** 3



def cercius2kelvin(x):
    return x + 273.15

def deg2rad(x):
    return x / 180 * numpy.pi

def rad2deg(x):
    return x / numpy.pi * 180


class Stirling_Engin:
    def __init__(self, V_SE, V_SC, V_D, T_E, T_C, alpha, N):
        self.V_SE = V_SE            # 膨張部(加熱部)の有効体積
        self.V_SC = V_SC            # 収縮部(冷却部)の有効堆積
        self.V_D = V_D              # 全無効体積
        self.T_E = T_E              # 膨張部(加熱部)温度
        self.T_C = T_C              # 収縮部(冷却部)温度
        self.alpha = deg2rad(alpha) # 偏角
        self.N = N                  # 回転数

    def T_R(self):
        return cercius2kelvin((self.T_E + self.T_C) / 2)

    def X(self):
        return self.V_D / self.V_SE

    def kappa(self):
        return self.V_SC / self.V_SE

    def S(self):
        return self.tau() + 4 * self.tau() * self.X() / (1 + self.tau()) + self.kappa()

    def B(self):
        return numpy.sqrt(self.tau() ** 2 + 2 * self.tau() * self.kappa() * numpy.cos(self.alpha) + self.kappa() ** 2)

    def tau(self):
        return cercius2kelvin(self.T_C) / cercius2kelvin(self.T_E)

    def delta(self):
        return self.B() / self.S()

    def V_E(self, theta):
        return 1/2 * self.V_SE * (1 - numpy.cos(theta))

    def V_C(self, theta):
        return 1/2 * self.V_SC * (1 - numpy.cos(theta - self.alpha))

    def V(self, theta):
        return self.V_E(theta) + self.V_D + self.V_C(theta)

    def phi(self):
        t = numpy.arctan(self.kappa() * numpy.sin(self.alpha) / (self.tau() + self.kappa() * numpy.cos(self.alpha)))
        return numpy.where(t >= 0, t, t + numpy.pi)

    def m(self):
        return Rho * self.V(numpy.pi / 2)

    def P(self, theta):
        return 2 * self.m() * R * cercius2kelvin(self.T_C) / (self.V_SE * (self.S() - self.B() * numpy.cos(theta - self.phi())))

    def P_mean(self):
        return 2 * self.m() * R * cercius2kelvin(self.T_C) / (self.V_SE * numpy.sqrt(self.S() ** 2 - self.B() ** 2))

    def W_E(self):
        return self.P_mean() * self.V_SE * numpy.pi * self.delta() * numpy.sin(self.phi()) / (1 + numpy.sqrt(1 - self. delta() ** 2))

    def W_C(self):
        return -1 * self.P_mean() * self.V_SE * numpy.pi * self.delta() * self.tau() * numpy.sin(self.phi()) / (1 + numpy.sqrt(1 - self.delta() ** 2))

    def W_i(self):
        return self.W_E() + self.W_C()

    def m_from_pmean(self, pmean):
        return self.V_SE * numpy.sqrt(self.S() ** 2 - self.B() ** 2) * pmean / (2 * R * cercius2kelvin(self.T_C))

    def L_i(self):
        return self.W_i() * self.N / 60

    def L(self):
        return self.L_i() * 0.6 * 0.5

    def show_details(self):
        print(f"高温部行程容積(V_SE):\t{self.V_SE/cm3} cm3")
        print(f"低温部行程容積(V_SC):\t{self.V_SC/cm3} cm3")
        print(f"全無効容積(V_D):\t{self.V_D/cm3} cm3")
        print()

        print(f"高温部温度(T_E):\t{self.T_E} ℃")
        print(f"低温部温度(T_C):\t{self.T_C} ℃")
        print()

        print(f"エンジン内ガス質量(m):\t{self.m() * 10 ** 6} mg")
        print()

        print(f"位相角(alpha):\t\t{rad2deg(self.alpha)} degree")
        print()

        print(f"理論図示仕事(W_i):\t{self.W_i()} W")


if __name__ == '__main__':
    # write your code here
    pass
