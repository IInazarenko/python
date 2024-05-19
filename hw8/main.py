import sympy as sp

if __name__ == '__main__':
    t = sp.symbols('t')
    e1 = sp.symbols('e1', cls=sp.Function)(t)
    e2 = sp.symbols('e2', cls=sp.Function)(t)
    x1 = sp.symbols('x1', cls=sp.Function)(t)
    x2 = sp.symbols('x2', cls=sp.Function)(t)
    g = sp.symbols('g')
    L = sp.symbols('L')
    m = sp.symbols('m')
    M = sp.symbols('M')
    u1 = sp.symbols('u1', cls=sp.Function)(t)
    u2 = sp.symbols('u2', cls=sp.Function)(t)

    H01 = sp.Matrix([[sp.cos(sp.pi/2), -sp.sin(sp.pi/2), 0, 0], [sp.sin(sp.pi/2), sp.cos(sp.pi/2), 0, 0],
                     [0, 0, 1, 0], [0, 0, 0, 1]]) @\
          sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, x1], [0, 0, 0, 1]]) @\
          sp.Matrix([[1, 0, 0, 0], [0, sp.cos(sp.pi / 2), -sp.sin(sp.pi / 2), 0],
                     [0, sp.sin(sp.pi / 2), sp.cos(sp.pi / 2), 0], [0, 0, 0, 1]])

    H12 = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, x2], [0, 0, 0, 1]]) @ \
          sp.Matrix([[1, 0, 0, 0], [0, sp.cos(-sp.pi / 2), -sp.sin(-sp.pi / 2), 0],
                     [0, sp.sin(-sp.pi / 2), sp.cos(-sp.pi / 2), 0], [0, 0, 0, 1]])

    H23 = sp.Matrix([[sp.cos(e1 - sp.pi/2), -sp.sin(e1 - sp.pi/2), 0, 0], [sp.sin(e1 - sp.pi/2), sp.cos(e1 - sp.pi/2), 0, 0],
                     [0, 0, 1, 0], [0, 0, 0, 1]]) @ \
          sp.Matrix([[1, 0, 0, 0], [0, sp.cos(-sp.pi / 2), -sp.sin(-sp.pi / 2), 0],
                     [0, sp.sin(-sp.pi / 2), sp.cos(-sp.pi / 2), 0], [0, 0, 0, 1]])

    H34 = sp.Matrix([[sp.cos(e2), -sp.sin(e2), 0, 0], [sp.sin(e2), sp.cos(e2), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) @\
          sp.Matrix([[1, 0, 0, L], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

    H = H01 @ H12 @ H23 @ H34

    rm = H @ sp.Matrix([[0], [0], [0], [1]])
    print(rm)
    # print(H)

    # print(rm)
    print(rm[0, 0])
    print(rm[1, 0])
    print(rm[2, 0])

    v = sp.diff(rm[0, 0], t)**2 + sp.diff(rm[1, 0], t)**2 + sp.diff(rm[2, 0], t)**2
    # print (v.simplify().expand())

    P = m * g * rm[1, 0]

    KM = 1 / 2 * M * (sp.diff(x1, t)**2 + sp.diff(x2, t)**2)
    Km = 1 / 2 * m * v;
    K = KM + Km

    L = K - P

    dL_dx1 = sp.diff(L, x1)
    x1dot = sp.diff(x1, t)
    x1dotdot = sp.diff(x1, t, 2)
    dL_dX1dot = sp.diff(L, x1dot)
    d_dt_dLd_dX1dot = sp.diff(dL_dX1dot, t)
    Lagrangian_x1 = d_dt_dLd_dX1dot - dL_dx1
    print("x1 ")
    print( Lagrangian_x1.simplify().expand())

    dL_de1 = sp.diff(L, e1)
    e1dot = sp.diff(e1, t)
    e1dotdot = sp.diff(e1, t, 2)
    dL_de1dot = sp.diff(L, e1dot)
    d_dt_dLd_de1dot = sp.diff(dL_de1dot, t)
    Lagrangian_e1 = d_dt_dLd_de1dot - dL_de1
    print("e1 ")
    print(Lagrangian_e1.simplify().expand())

    dL_dx2 = sp.diff(L, x2)
    x2dot = sp.diff(x2, t)
    x2dotdot = sp.diff(x2, t, 2)
    dL_dX2dot = sp.diff(L, x2dot)
    d_dt_dLd_dX2dot = sp.diff(dL_dX2dot, t)
    Lagrangian_x2 = d_dt_dLd_dX2dot - dL_dx2
    print("x2 ")
    print(Lagrangian_x2.simplify().expand())

    dL_de2 = sp.diff(L, e2)
    e2dot = sp.diff(e2, t)
    e2dotdot = sp.diff(e2, t, 2)
    dL_de2dot = sp.diff(L, e2dot)
    d_dt_dLd_de2dot = sp.diff(dL_de2dot, t)
    Lagrangian_e2 = d_dt_dLd_de2dot - dL_de2
    print("e2 ")
    print(Lagrangian_e2.simplify().expand())

    leq1 = Lagrangian_x1 - u1
    leq2 = Lagrangian_x2 - u2
    leq3 = Lagrangian_e1
    leq4 = Lagrangian_e2

    print("leq1")
    sp.print_latex(leq1.simplify())
    print("leq2")
    sp.print_latex(leq2.simplify())
    print("leq3")
    sp.print_latex(leq3.simplify())
    print("leq4")
    sp.print_latex(leq4.simplify())

    MatrixA = sp.Matrix([[sp.diff(leq1, x1dotdot), sp.diff(leq1, x2dotdot), sp.diff(leq1, e1dotdot), sp.diff(leq1, e2dotdot)],
                         [sp.diff(leq2, x1dotdot), sp.diff(leq2, x2dotdot), sp.diff(leq2, e1dotdot), sp.diff(leq2, e2dotdot)],
                         [sp.diff(leq3, x1dotdot), sp.diff(leq3, x2dotdot), sp.diff(leq3, e1dotdot), sp.diff(leq3, e2dotdot)],
                         [sp.diff(leq4, x1dotdot), sp.diff(leq4, x2dotdot), sp.diff(leq4, e1dotdot), sp.diff(leq4, e2dotdot)]])

    print("matrixA")
    sp.print_latex(MatrixA)
    MatrixB = sp.Matrix([[-(leq1 - sp.diff(leq1, x1dotdot)*x1dotdot - sp.diff(leq1, x2dotdot)*x2dotdot - sp.diff(leq1, e1dotdot)*e1dotdot - sp.diff(leq1, e2dotdot)*e2dotdot)],
                         [-(leq2 - sp.diff(leq2, x1dotdot)*x1dotdot - sp.diff(leq2, x2dotdot)*x2dotdot - sp.diff(leq2, e1dotdot)*e1dotdot - sp.diff(leq2, e2dotdot)*e2dotdot)],
                         [-(leq3 - sp.diff(leq3, x1dotdot)*x1dotdot - sp.diff(leq3, x2dotdot)*x2dotdot - sp.diff(leq3, e1dotdot)*e1dotdot - sp.diff(leq3, e2dotdot)*e2dotdot)],
                         [-(leq4 - sp.diff(leq4, x1dotdot)*x1dotdot - sp.diff(leq4, x2dotdot)*x2dotdot - sp.diff(leq4, e1dotdot)*e1dotdot - sp.diff(leq4, e2dotdot)*e2dotdot)]])
    print("matrixB")
    for item in MatrixB:
        print(" ")
        sp.print_latex(item.simplify())
    MatrixX = MatrixA.inv()@MatrixB
    print("X")
    print(MatrixX)
    x1dotdot = MatrixX[0, 0]
    x2dotdot = MatrixX[1, 0]
    e1dotdot = MatrixX[2, 0]
    e2dotdot = MatrixX[3, 0]

    a51 = sp.diff(x1dotdot, x1)
    a52 = sp.diff(x1dotdot, x2)
    a53 = sp.diff(x1dotdot, e1)
    a54 = sp.diff(x1dotdot, e2)
    a55 = sp.diff(x1dotdot, x1dot)
    a56 = sp.diff(x1dotdot, x2dot)
    a57 = sp.diff(x1dotdot, e1dot)
    a58 = sp.diff(x1dotdot, e2dot)
    u51 = sp.diff(x1dotdot, u1)
    u52 = sp.diff(x1dotdot, u2)

    a61 = sp.diff(x2dotdot, x1)
    a62 = sp.diff(x2dotdot, x2)
    a63 = sp.diff(x2dotdot, e1)
    a64 = sp.diff(x2dotdot, e2)
    a65 = sp.diff(x2dotdot, x1dot)
    a66 = sp.diff(x2dotdot, x2dot)
    a67 = sp.diff(x2dotdot, e1dot)
    a68 = sp.diff(x2dotdot, e2dot)
    u61 = sp.diff(x2dotdot, u1)
    u62 = sp.diff(x2dotdot, u2)

    a71 = sp.diff(e1dotdot, x1)
    a72 = sp.diff(e1dotdot, x2)
    a73 = sp.diff(e1dotdot, e1)
    a74 = sp.diff(e1dotdot, e2)
    a75 = sp.diff(e1dotdot, x1dot)
    a76 = sp.diff(e1dotdot, x2dot)
    a77 = sp.diff(e1dotdot, e1dot)
    a78 = sp.diff(e1dotdot, e2dot)
    u71 = sp.diff(e1dotdot, u1)
    u72 = sp.diff(e1dotdot, u2)

    a81 = sp.diff(e2dotdot, x1)
    a82 = sp.diff(e2dotdot, x2)
    a83 = sp.diff(e2dotdot, e1)
    a84 = sp.diff(e2dotdot, e2)
    a85 = sp.diff(e2dotdot, x1dot)
    a86 = sp.diff(e2dotdot, x2dot)
    a87 = sp.diff(e2dotdot, e1dot)
    a88 = sp.diff(e2dotdot, e2dot)
    u81 = sp.diff(e2dotdot, u1)
    u82 = sp.diff(e2dotdot, u2)

    A = sp.Matrix([[0,0,0,0,1,0,0,0], [0,0,0,0,0,1,0,0], [0,0,0,0,0,0,1,0], [0,0,0,0,0,0,0,1],
                   [a51, a52, a53, a54, a55, a56, a57, a58], [a61, a62, a63, a64, a65, a66, a67, a68],
                   [a71, a72, a73, a74, a75, a76, a77, a78], [a81, a82, a83, a84, a85, a86, a87, a88]]).subs({x1: 0, x2: 0, e1: sp.pi/2, e2: 0,
                                                                                                              x1dot: 0, x2dot: 0, e1dot: 0, e2dot: 0})
    B = sp.Matrix([[0, 0], [0, 0], [0, 0], [0, 0], [u51, u52], [u61, u62], [u71, u72], [u81, u82]]).subs({x1: 0, x2: 0, e1: sp.pi/2, e2: 0,
                                                                                                              x1dot: 0, x2dot: 0, e1dot: 0, e2dot: 0})
    print("A")
    sp.print_latex(A)
    print("B")
    sp.print_latex(B)

    control_mat = B.row_join(A @ B).row_join(A @ A @ B).row_join(A @ A @ A @ B).row_join(A @ A @ A @ A @ B).row_join(A @ A @ A @ A @ A @ B).row_join(A @ A @ A @ A @ A @ A @ B).row_join(A @ A @ A @ A @ A @ A @ A @ B)
    print("control ")
    sp.print_latex(control_mat)
    print(control_mat.rank())
