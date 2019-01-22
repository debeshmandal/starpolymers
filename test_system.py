import star_gen2 as sg
import base

x = sg.FileGenerator(40)

star = {'molecule': 'star',
        'kap': 10,
        'lam': 3,
        'charge_style': 'all',
        'charge_max': 1,
        'central' : 'all',
        'counterions':False,
        'angle_type': 1}

DNA = {'molecule': 'DNA',
       'kap': 1,
       'lam': 21,
       'charge_style': 'all',
       'charge_max': -1,
       'counterions': False,
       'angle_type': 2}

salt = {'molecule': 'salt',
        'kap': 0,
        'lam': 0,
        'charge_max': 1,
        'charge_style': 'all',
        'concentration': 1,
        'neutralise' : True}

brush = {'molecule':'brush',
         'trunk': {'lam':4},
         'branches' : [{'lam':2,
                        'site':2}],
         'start':[0,0,0.5],
         'direction':'up',
         'charge_style':'none',
         'counterions':False,
         'base_id': 12}

brush2 = {'molecule':'brush',
         'trunk': {'lam':4},
         'branches' : [{'lam':2,
                        'site':2}],
         'start':[2,0,0.5],
         'direction':'up',
         'charge_style':'none',
         'counterions':False,
         'base_id': 13}

surface_bottom = {'molecule' : 'base',
                  'dims': [5, 5],
                  'plane': 0,
                  'spacing':0.5,
                  'charge_style':'none'}

surface_top = {'molecule' : 'base',
               'dims': [5, 5],
               'plane': 10,
               'spacing':0.5,
               'charge_style':'none'}

system = [surface_top, surface_bottom, brush, brush2]
x.write_system_to_file(system, angles=False)
