import argparse
from random import randint

import numpy as np
from scipy.sparse import csc_matrix, rand


def main():
    parser = argparse.ArgumentParser()

    output = "matrix2.txt"
    shape = int(input())
    if shape < 3:
        exit()

    matrix = rand(shape, shape, density=0.4, random_state=randint(112, 154))
    matrix = matrix.toarray()
    for i in range(shape):
        for j in range(shape):
            matrix[j][i] = matrix[i][j]
    matrix = csc_matrix(matrix)
    with open(output, "w") as f:
        f.write(f"{shape}\n")
        for i in matrix.toarray().round(3):
            for j in i:
                f.write(f"{j} ")
            f.write("\n")
        d = np.random.randint(5, 53, shape)
        for i in d:
            f.write(f"{i} ")
        f.write("\n")


if __name__ == "__main__":
    main()
