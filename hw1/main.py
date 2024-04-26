from sympy import *
# print(2)
# t = Symbol('t')
# q = symbols('q', cls=Function)
# theta = symbols('theta', cls=Function)
# phi = symbols('phi', cls=Function)
# q = q(t)
# theta = theta (t)
# phi = phi(t)
# x = q * sin(theta) * cos(phi)
# y = q * sin(theta)*sin(phi)
# z = q * cos(theta)
# drdq = Matrix([[x.diff(q)], [y.diff(q)], [z.diff(q)]])
# drdphi = Matrix([[x.diff(phi)], [y.diff(phi)], [z.diff(phi)]])
# drtheta = Matrix([[x.diff(theta)], [y.diff(theta)], [z.diff(theta)]])
# print(drdq)
# print(drdphi)
# print(drtheta)
#
# K = 1/2 * (drdq.dot(drdq)*(q.diff(t))**2+drdphi.dot(drdphi)*phi.diff(t)**2+drtheta.dot(drtheta)*theta.diff(t)**2)
# K = simplify(K)
# print(K)
#
#
# dKdq = K.diff(q)
# dKdq_doted = K.diff(q.diff(t))
# print(dKdq)
# print(dKdq_doted)
#
# dKdphi = K.diff(phi)
# dKdphi_doted = K.diff(phi.diff(t))
# print(dKdphi)
# print(dKdphi_doted)
#
# dKdtheta = K.diff(theta)
# dKdtheta_doted = K.diff(theta.diff(t))
# print(dKdtheta)
# print(dKdtheta_doted)
#
#
# ddr_e1 = dKdq_doted.diff(t) - dKdq
# print(ddr_e1)
#
# ddr_e2 = dKdphi_doted.diff(t) - dKdphi
# print(ddr_e2)
#
# ddr_e3 = dKdtheta_doted.diff(t) - dKdtheta
# print(ddr_e3)
#

#part 3

t = Symbol('t')
tau = symbols('tau', cls=Function)
tau = tau(t)
tau_0 = 1
tau_equation = Eq(tau.diff(t), cos(t)*tau)
x = tau * sin(2*tau)*cos(3*tau)
y = tau * sin(2*tau)*sin(3*tau)
z = tau * cos(2 * tau)

print(x)
print(y)
print(z)

solution = dsolve(tau_equation, tau, ics={tau.subs(t, 0): tau_0})

print("Solution for tau(t):")
print(solution)


