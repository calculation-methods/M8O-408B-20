import numpy as np

t_max = 1

count_x = 60
count_t = 1000

r_coord = np.pi / 2

a = 0.1
c = -1

h = r_coord / count_x
tau = t_max / count_t

sigma = a * tau / (h**2)
