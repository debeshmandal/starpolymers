from starpolymers.generators.input_file import templates, Variable
import pandas as pd

def test_01():
    input_file = templates['01_master']
    variables = input_file.variables
    input_file.update(Variable('k1', 10.0))
    input_file.string

if __name__ == '__main__':
    test_01()