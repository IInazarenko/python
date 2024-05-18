from scipy.integrate import ode
from math import sin, cos, pi
import numpy as np
import matplotlib.pyplot as plt
def dynamics(t, state):
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


solver = ode(dynamics).set_integrator('dopri5', atol=1e-8, rtol=1e-8, max_step=1e-3, nsteps=1000)
solver = solver.set_initial_value([0, 0])
theta = [0]
thetadot = [0]

for t in zip(np.arange(0, 15, 0.001).tolist()[1:]):
            th, dth = solver.integrate(t)
            theta.append(th)
            thetadot.append(dth)


t_array = np.arange(0, 15, 0.001)
ref_array = (pi + 0.3 * np.sin(2*t_array)).tolist()
t_array = t_array.tolist()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t_array, ref_array, "-g")
ax.plot(t_array, theta, "-r")
ax.set_xlabel("t")
ax.set_ylabel("Theta")
plt.show()