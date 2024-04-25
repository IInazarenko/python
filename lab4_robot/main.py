import sympy as sp

def RPR_FK(theta1, d2, theta3):
    a =10
    alpha = -135
    # Define transformation matrix A01
    A01 = sp.Matrix([
        [sp.cos(theta1), -sp.sin(theta1), 0, 0],
        [sp.sin(theta1), sp.cos(theta1), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ]) @ sp.Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ]) @ sp.Matrix([
        [1, 0, 0, a],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ]) @ sp.Matrix([
        [1, 0, 0, 0],
        [0, sp.cos(alpha), -sp.sin(alpha), 0],
        [0, sp.sin(alpha), sp.cos(alpha), 0],
        [0, 0, 0, 1]
    ])

    # Print the result
    print(A01)
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    RPR_FK(5,6,7)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
