import numpy as np

data = {'equation_type': 'fractionalSteps', 'nx': 40, 'ny': 40, 'T': 5, 'K': 200}
equation_type = data['equation_type']
nx, ny, T, K = int(data['nx']), int(data['ny']), int(data['T']), int(data['K'])
params = {
    'a': 1,
    'b': 1,
    'c': 0,
    'd': 0,
    'lx': 1,
    'ly': 1,
    'f': lambda x, y, t: - x * y * np.sin(t),
    'alpha1': 0,
    'alpha2': 1,
    'beta1': 0,
    'beta2': 1,
    'gamma1': 0,
    'gamma2': 1,
    'delta1': 0,
    'delta2': 1,
    'phi11': lambda y, t: 0,
    'phi12': lambda y, t: y * np.cos(t),
    'phi21': lambda x, t: 0,
    'phi22': lambda x, t: x * np.cos(t),
    'psi': lambda x, y: x * y,
    'solution': lambda x, y, t: x * y * np.cos(t),
}