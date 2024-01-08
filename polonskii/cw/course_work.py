import matplotlib.pyplot as plt
import numpy as np


def f1(x):
    return np.cos(2 * x)


def f2(x):
    if x <= 10:
        return np.sin(x) + np.cos(2 * x)
    else:
        return np.sin(x * 0.3) + np.cos(6 * x)


def f3(x):
    return np.cos(x ** 2)


def f4():
    return np.cos(2 * np.pi * (2 ** np.linspace(2, 10, N)) * np.arange(N) / 48000) + np.random.normal(0, 1, N) * 0.15


def discreteWaveletTransform(tabData):
    size = len(tabData) // 2
    cA = np.zeros(size)
    cD = np.zeros(size)

    for i, j in zip(range(0, len(tabData), 2), range(size)):
        c = 2 * (tabData[i] + tabData[i + 1]) / np.sqrt(N)
        cA[j] = c

    for i, j in zip(range(0, len(tabData), 2), range(size)):
        c = 2 * (tabData[i] - tabData[i + 1]) / np.sqrt(N)
        cD[j] = c

    return cA, cD


def waveletDeconstr(tabData, level=1):  # Returns [cA_n, cD_n, cD_n-1, ..., cD2, cD1]
    coeffs = []
    a = tabData
    for i in range(level):
        a, d = discreteWaveletTransform(a)
        coeffs.append(d)
    coeffs.append(a)
    coeffs.reverse()

    return coeffs


def inverseDWT(a, d):
    res = []
    for i in range(len(a)):
        x = (a[i] + d[i]) * np.sqrt(N) / 4
        y = (a[i] - d[i]) * np.sqrt(N) / 4
        res.extend([x, y])
    return np.array(res)


def waveletReconstr(coeffs):
    a, ds = coeffs[0], coeffs[1:]
    for d in ds:
        a = inverseDWT(a, d)

    return a


scale = int(input())
N = int(2**scale)
level = 3
start = 0
stop = N
X = np.linspace(0, stop, N)
# Y = [f3(x) for x in X]
Y = f4()
c = waveletDeconstr(Y, level)
X1 = np.linspace(0, int(stop / (2**level)), int(N / (2**level)))

figure, axis = plt.subplots(2, 2)
axis[0, 0].plot(X, Y)
axis[0, 0].set_title('Original signal')
axis[1, 0].plot(X1, c[0])
axis[1, 0].set_title(f'Approximation Coefficients, level={level}')
axis[1, 1].plot(X1, c[1], c='orange')
axis[1, 1].set_title('Detail Coefficients')

inv = waveletReconstr(c)
axis[0, 1].plot(X, inv)
axis[0, 1].set_title('Reconstructed signal')
plt.show()

