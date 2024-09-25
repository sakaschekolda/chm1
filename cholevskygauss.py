import cmath
import os

os.system('cls')

def gauss_solve(A, b):
    n = len(A)

    for i in range(n):
        max_el = abs(A[i][i])
        max_row = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > max_el:
                max_el = abs(A[k][i])
                max_row = k

        A[i], A[max_row] = A[max_row], A[i]
        b[i], b[max_row] = b[max_row], b[i]

        for k in range(i + 1, n):
            c = -A[k][i] / A[i][i]
            for j in range(i, n):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]
            b[k] += c * b[i]

    x = [0 for _ in range(n)]
    for i in range(n - 1, -1, -1):
        x[i] = b[i] / A[i][i]
        for k in range(i - 1, -1, -1):
            b[k] -= A[k][i] * x[i]

    return x


def cholesky_solve(A, b):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            sum = 0
            for k in range(j):
                sum += L[i][k] * L[j][k]
            if i == j:
                L[i][i] = cmath.sqrt(A[i][i] - sum)
            else:
                L[i][j] = (A[i][j] - sum) / L[j][j]

    y = [0.0] * n
    for i in range(n):
        sum = 0
        for k in range(i):
            sum += L[i][k] * y[k]
        y[i] = (b[i] - sum) / L[i][i]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum = 0
        for k in range(i + 1, n):
            sum += L[k][i] * x[k]
        x[i] = (y[i] - sum) / L[i][i]

    return x

matrix = [
    [2.93, 2.55, 2.14],
    [3.47, 2.98, 2.50],
    [4.78, 4.22, 3.7]
]

vector = [46.41, 54.78, 75.81]

solutiongauss = gauss_solve(matrix, vector)

print('gauss:\n')
for i in range(len(solutiongauss)):
    print(f"x{i + 1} = {solutiongauss[i]}")
    
solutioncholevsky = cholesky_solve(matrix, vector)

print('\ncholevsky:\n')
for i in range(len(solutioncholevsky)):
    print(f"x{i + 1} = {solutioncholevsky[i]}")