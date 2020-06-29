import pandas as pd
import numpy as np
import random
from argparse import ArgumentParser

ranges = {
    'kap' : np.arange(2, 9, 1),
    'lam' : np.arange(5, 26),
    'charge_max' : np.arange(0.0, 4.0, 1.0),
    'charge_ratio' : np.arange(0.0, 1.1, 0.1),
    'charge_position' : ['core', 'tail'],
    'K': [5., 10., 20., 50.],
    'theta0': [90., 120., 135., 180.0],
    'k1': [0., 1., 5., 10., 100.],
    'k2': [0., 1., 5., 10., 100.]
}

def generate_table(n, seed=1994, write=False, fout='table.csv', show=False):
    data = pd.DataFrame(columns=list(ranges.keys()))
    random.seed(seed)
    for i in range(n):
        settings = ranges.copy()
        for key, value in settings.items():
            settings[key] = random.sample(list(value), 1)
        data = pd.concat(
                [
                    data,
                    pd.DataFrame(settings)
                ],
                sort=False
            ).reset_index(drop=True)
    if show:
        print(data)
    if write:
        data.to_csv(fout, index=False)
    return data

if __name__ == '__main__':
    parser = ArgumentParser(
        description="Generate a DoE table based on randomly sampled values found in the file"
    )
    parser.add_argument('-n', '--number', default=100, type=int)
    parser.add_argument('-w', '--write', action='store_true', default=False)
    parser.add_argument('--file', default='./table.csv')
    parser.add_argument('-s', '--seed', default=1994, type=int)
    parser.add_argument('--show', action='store_true', default=False)
    args = parser.parse_args()

    generate_table(
        args.number,
        write=args.write,
        fout=args.file,
        seed=args.seed,
        show=args.show
    )
