from scipy.integrate import ode
from math import sin, cos, pi, sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


def soft_sign(x):
    eps = 1e-2
    return x / (abs(x) + eps)


def dynamics1(t, state):
    theta, thetadot = state
    thetaR = pi + 0.3 * sin(2*t)
    thetaRdot = 0.6*cos(2*t)
    thetaRdotdot = -1.2*sin(2*t)
    M = 1 - 1/2*cos(theta)**2
    C = 1/2* cos(theta)*sin(theta)*thetadot
    Greal = 4.5*sin(theta)
    Gcontrol = 4*sin(theta)
    dM = cos(theta)*sin(theta)*thetadot
    Kp = 1
    Kd = 1
    u = C*thetadot+Gcontrol+M*thetaRdotdot - 1/2 * dM * (thetadot - thetaRdot) - Kp * (theta - thetaR) - Kd * (thetadot - thetaRdot)
    thetadotdot = (u - thetadot*C - Greal)/M
    return [thetadot, thetadotdot]

def dynamics2(t, state):
    theta, thetadot = state
    thetaR = pi + 0.3 * sin(2*t)
    thetaRdot = 0.6*cos(2*t)
    thetaRdotdot = -1.2*sin(2*t)
    M = 1 - 1/2*cos(theta)**2
    C = 1/2 * cos(theta)*sin(theta)*thetadot
    Greal = 4.5 *sin(theta)
    Gcontrol = 4*sin(theta)
    Kp = 1
    Kd = 2
    A = np.array([
        [0, 1],
        [-Kp, -Kd]
    ])
    B = np.array([
        [0],
        [1]
    ])
    E = np.array([
        [theta - thetaR],
        [thetadot - thetaRdot]])
    P = np.array([
        [3/2, 1/2],
        [1/2, 1/2]])
    z = B.transpose()@P@E
    if z == 0:
        w = 0
    else:
        w = - 2 * soft_sign(z[0, 0])
    v = thetaRdotdot - Kp*(theta - thetaR) - Kd*(thetadot - thetaRdot) + w
    u = M*v + C*thetadot+Gcontrol
    thetadotdot = (u - thetadot*C - Greal)/M
    return [thetadot, thetadotdot]

def dynamics3(t, state):
    theta, thetadot, est1, est2, est3 = state
    thetaR = pi + 0.3 * sin(2*t)
    thetaRdot = 0.6*cos(2*t)
    thetaRdotdot = -1.2*sin(2*t)
    M = 1 - 1/2*cos(theta)**2
    C = 1/2* cos(theta)*sin(theta)*thetadot
    Greal = 4.5 * sin(theta)
    Gcontrol = est3* sin(theta)
    Kp = 1
    Kd = 2
    A = np.array([
        [0, 1],
        [-Kp, -Kd]
    ])
    B = np.array([
        [0],
        [1]
    ])
    E = np.array([
        [theta - thetaR],
        [thetadot - thetaRdot]])
    P = np.array([
        [3 / 2, 1 / 2],
        [1 / 2, 1 / 2]])
    Ge = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])

    v = thetaRdotdot - Kp*(theta - thetaR) - Kd*(thetadot - thetaRdot)
    u = M*v + C*thetadot+Gcontrol
    thetadotdot = (u - thetadot*C - Greal)/M
    F = np.array([[thetadotdot*(1 - 1/2*(cos(theta))**2), thetadot**2 *(1/2*cos(theta)*sin(theta)), sin(theta)]])
    estimationdot = - Ge@F.transpose()@B.transpose()@P@E
    return [thetadot, thetadotdot, estimationdot[0, 0], estimationdot[1, 0], estimationdot[2, 0]]


def simulation(dynamics):
    solver = ode(dynamics).set_integrator('dopri5', atol=1e-8, rtol=1e-8, max_step=1e-3, nsteps=1000)
    solver = solver.set_initial_value([0, 0, 1, 1, 4])
    theta = [0]
    thetadot = [0]
    estimation1 = [1]
    estimation2 = [1]
    estimation3 = [4]

    for t in zip(np.arange(0, 30, 0.001).tolist()[1:]):
            th, dth, est1, est2, est3 = solver.integrate(t)
            theta.append(th)
            thetadot.append(dth)
            estimation1.append(est1)
            estimation2.append(est2)
            estimation3.append(est3)


    t_array = np.arange(0, 30, 0.001)
    ref_array = (pi + 0.3 * np.sin(2*t_array)).tolist()
    t_array = t_array.tolist()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t_array, ref_array, "-g")
    ax.plot(t_array, theta, "-r")
    ax.set_xlabel("t")
    ax.set_ylabel("Theta")
    plt.show()


# simulation(dynamics1)
# simulation(dynamics2)
simulation(dynamics3)