import numpy as np

t_max = 1
count_x = 50
count_t = 1000
r_coord = np.pi
a = 1
b = 2
c = -3
e = 2
h = r_coord / count_x
tau = t_max / count_t
# число Куранта
sigma = a * tau / (h**2)
