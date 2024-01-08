import argparse
import numpy as np
from scipy.sparse import csc_matrix, rand
from random import randint


class SparseMatrix:
    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self.elements = []

    def add_element(self, row, col, value):
        self.elements.append((row, col, value))

    def set_b(self, b_values):
        self.b = b_values


class MatrixGenerator:
    def generate(self, shape):
        matrix = rand(shape, shape, density=0.4, random_state=randint(112, 154))
        matrix = matrix.toarray()
        for i in range(shape):
            for j in range(shape):
                matrix[j][i] = matrix[i][j]
        matrix = csc_matrix(matrix)
        return matrix


def save_matrix_to_file(matrix, b, output_file):
    with open(output_file, "w") as f:
        f.write(f"{matrix.shape[0]}\n")
        for row in matrix.toarray().round(3):
            f.write(" ".join(map(str, row)) + "\n")
        f.write(" ".join(map(str, b)) + "\n")


def main():
    parser = argparse.ArgumentParser()

    output_file = "matrix.txt"
    shape = int(input())
    if shape < 3:
        exit()

    matrix_generator = MatrixGenerator()
    matrix = matrix_generator.generate(shape)
    
    sparse_matrix = SparseMatrix(matrix.shape[0], matrix.shape[1])
    for i in range(shape):
        for j in range(shape):
            sparse_matrix.add_element(i, j, matrix[i, j])

    d = np.random.randint(5, 53, shape)
    sparse_matrix.set_b(d)

    save_matrix_to_file(matrix, d, output_file)


if __name__ == "__main__":
    main()
