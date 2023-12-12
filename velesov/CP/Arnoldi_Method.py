import numpy as np


def qr_alg(matrix, tol=1e-12, max_iter=1000):
    n = matrix.shape[0]
    eigenvalues = np.zeros(n, dtype=complex)

    for i in range(max_iter):

        shift = matrix[-1, -1]


        Q, R = custom_qr(matrix - shift * np.eye(n))
        matrix = R @ Q + shift * np.eye(n)

        # Check for convergence
        if frobenius_norm_upper_triangle(matrix) < tol:
            break

        eigenvalues = np.diag(matrix)

    return eigenvalues


def custom_qr(matrix):
    n = matrix.shape[0]
    Q = np.eye(n)
    R = matrix.copy()

    for i in range(n - 1):
        x = R[i:, i].copy()
        norm_x = np.linalg.norm(x)
        v = x - norm_x * np.eye(len(x))[:, 0]
        v = v / np.linalg.norm(v)

        # Construct the Householder matrix
        H = np.eye(len(x)) - 2 * np.outer(v, v)

        # Update R and Q
        R[i:, i:] = H @ R[i:, i:]
        Q[:, i:] = Q[:, i:] @ H.T

    return Q, R


def frobenius_norm_upper_triangle(matrix):
    return np.linalg.norm(np.triu(matrix, k=1), ord='fro')

#Метод Арнольди
def arnoldi_method(matrix, k):
    """Computes a basis of the (k + 1)-Krylov subspace of matrix

    Arguments
      matrix: n × n array
      k: dimension of Krylov subspace

    Returns
      Q: n x (k + 1) array, the columns are an orthonormal basis of the
        Krylov subspace.
      h: (n + 1) x n array, A on basis Q. It is upper Hessenberg.
    """
    n = len(matrix)

    # Инициализация начального вектора
    q = np.random.rand(n)
    q = q / np.linalg.norm(q)

    # Массивы для хранения ортогональных векторов и верхнетреугольной матрицы H
    Q = np.zeros((n, k + 1))
    H = np.zeros((k + 1, k))

    # Первый столбец матрицы Q - начальный вектор
    Q[:, 0] = q

    for j in range(k):
        # Вычисление нового вектора
        v = np.dot(matrix, Q[:, j])

        # Процесс ортогонализации по методу Грама-Шмидта
        for i in range(j + 1):
            H[i, j] = np.dot(Q[:, i], v)
            v = v - H[i, j] * Q[:, i]

        # Нормализация нового вектора
        H[j + 1, j] = np.linalg.norm(v)
        Q[:, j + 1] = v / H[j + 1, j]

    # Compute eigenvalues without using np.linalg.eigvals()
    eigenvalues_H = qr_alg(H[:k, :k])

    # Вычисление приближенных собственных векторов матрицы A
    eigenvectors_A = np.dot(Q[:, :k], np.linalg.inv(Q[:, :k].T @ Q[:, :k]) @ Q[:, :k].T)

    return eigenvalues_H, eigenvectors_A


# Получаем матрицу от пользователя
n = int(input("Введите размер матрицы: "))
matrix = np.zeros((n, n))
print("Введите элементы матрицы построчно:")

for i in range(n):
    matrix[i, :] = list(map(float, input().split()))


k = n

# Вызываем функцию метода Арнольди
eigenvalues, eigenvectors = arnoldi_method(matrix, k)

# Выводим результаты
print("\nМатрица:")
print(matrix)
print("\nВсе собственные значения:")
print(eigenvalues)
print("\nПриближенные собственные векторы:")
print(eigenvectors)
