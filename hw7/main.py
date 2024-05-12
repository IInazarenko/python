from sympy import *
from sympy.vector import cross
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Set the backend before importing pyplot
import matplotlib.pyplot as plt


if __name__ == '__main__':
    t = symbols('t')
    theta = symbols('theta', cls=Function)(t)
    x = symbols('x', cls=Function)(t)
    k = symbols('k')
    m2 = symbols('m2')
    g = symbols('g')
    m1 = symbols('m1')
    length = symbols('l')
    # xdot = symbols('xdot')
    # thetadot = symbols('thetadot')

    Pc = k*x**2 / 2
    Pp = -m2*g*cos(theta)*length
    P = Pc + Pp

    Kc = 1/2*m1*diff(x, t)**2
    Kp = 1/2*m2*(diff(x, t)**2 + 2*diff(x, t)*length*diff(theta, t)*cos(theta)+diff(theta, t)**2*length**2)
    K = Kc + Kp

    L = K - P

    dL_dx = diff(L, x)
    xdot = diff(x, t)
    dL_dXdot = diff(L, xdot)
    d_dt_dLd_dXdot = diff(dL_dXdot, t)
    Lagrangian_x = d_dt_dLd_dXdot - dL_dx

    dL_dtheta = diff(L, theta)
    thetadot = diff(theta, t)
    dL_dThetadot = diff(L, thetadot)
    d_dt_dLd_dThetadot = diff(dL_dThetadot, t)
    Lagrangian_theta = d_dt_dLd_dThetadot - dL_dtheta




    print(Lagrangian_x.expand())
    print(Lagrangian_theta.expand())
