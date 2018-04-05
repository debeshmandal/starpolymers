import starpolymers.star_gen2 as sg
import subprocess

x = sg.FileGenerator(2000, None)

star = {'molecule': 'star',
        'kap': 10,
        'lam': 3,
        'charge_style': 'all',
        'charge_max': 1,
        'central' : 'all',
        'counterions':False}

DNA = {'molecule': 'DNA',
       'kap': 1,
       'lam': 21,
       'charge_style': 'all',
       'charge_max': -1,
       'counterions': False}

salt = {'molecule': 'salt',
        'kap': 0,
        'lam': 0,
        'charge_max': 1,
        'charge_style': 'all',
        'concentration': 200,
        'neutralise' : True}

system = [star, salt]
x.write_system_to_file(system)
