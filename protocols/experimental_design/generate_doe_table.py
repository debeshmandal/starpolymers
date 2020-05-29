import pandas as pd
import numpy as np
import random

ranges = {
    'kap' : np.arange(2, 9, 1),
    'lam' : np.arange(5, 36),
    'charge_max' : np.arange(0.0, 4.0, 1.0),
    'charge_ratio' : np.arange(0.0, 1.1, 0.1),
    'charge_position' : ['core', 'tail'],
    'K': [5., 10., 20., 50.],
    'theta0': [90., 120., 135., 180.0],
    'k1': [0., 1., 5., 10., 100.],
    'k2': [0., 1., 5., 10., 100.]
}

def generate_table(n, seed=1994):
    data = pd.DataFrame(columns=ranges.keys())
    random.seed(seed)
    for i in range(n):
        settings = dict(ranges)
        for key, value in settings.items():
            settings[key] = random.sample(value, 1)
        data = pd.concat(
                [
                    data,
                    pd.DataFrame(settings)
                ],
                sort=False
            ).reset_index(drop=True)
    print data
    return data

if __name__ == '__main__':
    generate_table(10)