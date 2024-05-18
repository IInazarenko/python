from sympy import *
from sympy.vector import cross
import numpy as np
from scipy import linalg
import matplotlib
matplotlib.use('TkAgg')
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
    f = symbols('F')
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


    leq1 = Lagrangian_x.expand()
    leq2 = Lagrangian_theta.expand()

    print(leq1)
    print(leq2, "\n")

    system_of_equations = [Eq(leq1, f), Eq(leq2, 0)]
    solutions = solve(system_of_equations, (diff(x,t,2), diff(theta,t,2)))
    # print(solutions)
    xdotdot     = solutions[list(solutions.keys())[1]]
    thetadotdot = solutions[list(solutions.keys())[0]]

    print_latex(thetadotdot.simplify())
    print_latex(xdotdot.simplify())
    print("\n\n\n\n")

    a31 = diff(xdotdot, x)
    a32 = diff(xdotdot, theta)
    a33 = diff(xdotdot, xdot)
    a34 = diff(xdotdot, thetadot)
    u3 =  diff(xdotdot, f)

    a41 = diff(thetadotdot, x)
    a42 = diff(thetadotdot, theta)
    a43 = diff(thetadotdot, xdot)
    a44 = diff(thetadotdot, thetadot)
    u4 =  diff(thetadotdot, f)


    A = Matrix([[0,0,1,0],[0,0,0,1],[a31, a32, a33, a34],[a41, a42, a43, a44]])
    A4 = A
    A = A.subs({x:0, theta: pi, xdot: 0, thetadot:0}).expand()
    B = Matrix([[0],[0],[u3],[u4]]).subs({x:0, theta: pi, xdot: 0, thetadot:0})

    print_latex(A)
    print_latex(B)

    control_mat = B.row_join(A@B).row_join(A@A@B).row_join(A@A@A@B).expand()
    print_latex(control_mat)

    Anp = np.array(A.subs({length: 0.8, k: 25, m1: 1, m2: 0.5, g: 9.8})).astype(np.float64)
    Bnp = np.array(B.subs({length: 0.8, k: 25, m1: 1, m2: 0.5, g: 9.8})).astype(np.float64)
    print(Anp)
    print(Bnp)

    R = 1
    Q = np.eye(4)
    K =  - Bnp.T @ linalg.solve_continuous_are(Anp, Bnp, Q, R)

    print(K)

    good_list = []
    bad_list = []

    for ttheta in np.arange(0,2*pi, pi/5):
        for tx in np.arange(-1,1, 1/5):
            for tdtheta in np.arange(-1,1, 1/5):
                At = A4.subs({x:tx, theta:ttheta, thetadot:tdtheta, xdot: 0, length: 0.8, k: 25, m1: 1, m2: 0.5, g: 9.8, f:0})
                Atnp = np.array(At).astype(np.float64)
                M = Atnp + Bnp @ K
                val = np.linalg.eigvals(M).real
                val = val < 0
                if val.prod():
                    good_list.append([tx, ttheta, tdtheta])
                else:
                    bad_list.append([tx, ttheta, tdtheta])

    bad_list = zip(*bad_list)
    good_list = zip(*good_list)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.plot(*bad_list, ".r")
    ax.plot(*good_list, ".g")
    ax.set_xlabel("X")
    ax.set_ylabel("Theta")
    ax.set_zlabel("dTheta/dt")
    plt.show()
