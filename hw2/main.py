from sympy import *
from sympy.vector import cross
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Set the backend before importing pyplot
import matplotlib.pyplot as plt


if __name__ == '__main__':
    alpha = np.pi/600  # значение alpha
    t = symbols('t')
    RA = Matrix([[50 * (2 * cos(alpha * t) - 3 * cos(2 * alpha * t))],
                [50 * (2 * sin(alpha * t) - 3 * sin(2 * alpha * t))],
                [5 * sin(2 * alpha * t)]])

    # print(RA)
    # dRA = RA.diff(t)
    # ddRA = dRA.diff(t)
    # e1 = dRA/sqrt(10000 + 150*150 - 15000 * cos (alpha * t)+ 25 * sin(2*alpha*t) * sin(2*alpha*t))
    # print("e1:")
    # print(e1)
    # e2 = cross(dRA, ddRA)
    # print(e2)
    #
    # e2 = e2/e2.norm(2)
    # e3 = cross(e1, e2)
    #
    # print("e2:")
    # print(e2)
    # print("e3:")
    # print(e3)
    print("s:")
    s = integrate(RA.diff(t).norm(2), (t, 0, 1200))
    s_numeric = s.evalf()
    print(s_numeric)
    # w = 1/2 * (cross(e1, e1.diff(t)) + cross (e2, e2.diff(t)) + cross(e3, e3.diff(t)))
    print(cross(e1, e1.diff(t)))
    print(cross(e2, e2.diff(t)))
    print(cross(e3, e3.diff(t)))

# # Определяем параметры
# t_values = np.linspace(0, 1200, 12000)  # значения t от 0 до 10
# alpha = np.pi/600  # значение alpha
#
# # Определяем функции для каждой компоненты вектора RA
# def component1(t):
#     return 50 * (2 * np.cos(alpha * t) - 3 * np.cos(2 * alpha * t))
#
# def component2(t):
#     return 50 * (2 * np.sin(alpha * t) - 3 * np.sin(2 * alpha * t))
#
# def component3(t):
#     return 5 * np.sin(2 * alpha * t)
#
# # Вычисляем значения компонент вектора для каждого t
# RA_component1 = component1(t_values)
# RA_component2 = component2(t_values)
# RA_component3 = component3(t_values)
#
# # Визуализируем вектор RA
# plt.figure(figsize=(8, 6))
# plt.plot(t_values, RA_component1, label='$50(2\\cos(\\alpha t) - 3\\cos(2\\alpha t))i_1$')
# plt.plot(t_values, RA_component2, label='$50(2\\sin(\\alpha t) - 3\\sin(2\\alpha t))i_2$')
# plt.plot(t_values, RA_component3, label='$5\\sin(2\\alpha t)i_3$')
# plt.xlabel('$t$')
# plt.ylabel('$RA$')
# plt.title('Vector RA')
# plt.legend()
# plt.grid(True)
# plt.show()