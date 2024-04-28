from sympy import *
# t = Symbol('t')
# q = symbols('q', cls=Function)
# theta = symbols('theta', cls=Function)
# phi = symbols('phi', cls=Function)
# q = q(t)
# theta = theta(t)
# phi = phi(t)
# x = q * sin(theta) * cos(phi)
# y = q * sin(theta)*sin(phi)
# z = q * cos(theta)
# drdq = Matrix([[x.diff(q)], [y.diff(q)], [z.diff(q)]])
# drdphi = Matrix([[x.diff(phi)], [y.diff(phi)], [z.diff(phi)]])
# drtheta = Matrix([[x.diff(theta)], [y.diff(theta)], [z.diff(theta)]])
#
# print("e1:")
# pprint(drdq)
#
# print("e2:")
# pprint(drtheta)
#
# print("e3:")
# pprint(drdphi)
#
# K = 1/2 * (drdq.dot(drdq)*(q.diff(t))**2+drdphi.dot(drdphi)*(phi.diff(t))**2+drtheta.dot(drtheta)*(theta.diff(t))**2)
#
# dKdq = K.diff(q)
# dKdq_doted = K.diff(q.diff(t))
#
# dKdtheta = K.diff(theta)
# dKdtheta_doted = K.diff(theta.diff(t))
#
# dKdphi = K.diff(phi)
# dKdphi_doted = K.diff(phi.diff(t))
#
#
# ddr_e1 = dKdq_doted.diff(t) - dKdq
# print("ddr * e1:")
# print(expand(simplify(ddr_e1)))
#
# ddr_e2 = dKdtheta_doted.diff(t) - dKdtheta
# print("ddr * e2:")
# print(expand(simplify(ddr_e2)))
#
# ddr_e3 = dKdphi_doted.diff(t) - dKdphi
# print("ddr * e3:")
# print(expand(simplify(ddr_e3)))



# part 3

t = Symbol('t')
tau = symbols('tau', cls=Function)
tau = tau(t)
tau_0 = 1
tau_equation = Eq(tau.diff(t), cos(t)*tau)
x = tau * sin(2*tau)*cos(3*tau)
y = tau * sin(2*tau)*sin(3*tau)
z = tau * cos(2 * tau)

solution = dsolve(tau_equation, tau, ics={tau.subs(t, 0): tau_0})

print("Solution for tau(t):")
print(solution)

x = x.subs(tau, solution.rhs)
y = y.subs(tau, solution.rhs)
z = z.subs(tau, solution.rhs)

drdt = sqrt(x.diff(t)**2 + y.diff(t)**2 + z.diff(t)**2)

integral_result_0_1 = integrate(drdt, (t, 0, 1))
integral_result_simplified_0_1 = nsimplify(integral_result_0_1)
integral_result_numeric_0_1 = integral_result_simplified_0_1.evalf()
# Вывод числового значения интеграла
print("Numeric value of the integral from t=0 to t=1:")
print(integral_result_numeric_0_1)

integral_result_2_3 = integrate(drdt, (t, 2, 3))
integral_result_simplified_2_3 = nsimplify(integral_result_2_3)
integral_result_numeric_2_3 = integral_result_simplified_2_3.evalf()
# Вывод числового значения интеграла
print("Numeric value of the integral from t=2 to t=3:")
print(integral_result_numeric_2_3)