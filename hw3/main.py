import numpy as np
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from math import sin, cos

s = 0
dt = 0.001
T = 3
rad = 0.5

r = R.from_matrix([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1]])

state = [s, r]


# def omega(t):
#     return np.array([[1],
#                     [ 0],
#                     [ 0]])

def omega(t):
    return np.array([[cos(t)],
                    [ sin(t)],
                    [ 0]])

# def omega(t):
#     return np.array([[cos(2 * t) ** 3],
#                      [sin(2 * t) ** 3],
#                      [0]])


def skew(om):
    return np.array([[0, -om[2, 0], om[1, 0]], [om[2, 0], 0, -om[0, 0]], [-om[1, 0], om[0, 0], 0]])


def solver(t, state):
    s, r = state
    skew_om = skew(omega(t))
    dr = skew_om @ r.as_matrix()
    r = R.from_matrix(r.as_matrix() + dr * dt)

    vec = np.array([[0],
                    [0],
                    [rad]])

    vel = skew_om @ vec
    s += np.linalg.norm(vel * dt)

    state = [s, r]

    return state


for i in range(int(T / dt)):
    state = solver(i * dt, state)

sk, r = state

omega_new = np.transpose(r.as_matrix()) @ omega(3)
print(omega_new)




