import star_gen2 as sg
import subprocess

x = sg.FileGenerator()

star = {'molecule': 'star',
        'kap': 3,
        'lam': 10,
        'charge_style': 'all',
        'charge_max': 1}

DNA = {'molecule': 'DNA',
       'kap': 1,
       'lam': 21,
       'charge_style': 'all',
       'charge_max': -1,
       'counterions': True}

salt = {'molecule': 'salt',
        'kap': 0,
        'lam': 0,
        'charge_max': 1,
        'concentration': 10}

neutralising_ions = {'molecule': 'salt',
        'kap': 0,
        'lam': 0,
        'charge_max': 1,
        'concentration': 'neutralise'}

x.write_system_to_file()
