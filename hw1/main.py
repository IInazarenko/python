import sympy as sp
print(2)
q = sp.Symbol('q')
theta = sp.Symbol('theta')
phi = sp.Symbol('phi')
x = q * sp.sin(theta) * sp.cos (phi)
print(x.diff(phi))