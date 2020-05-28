from starpolymers.generators.input_file import templates, Variable
import pandas as pd

def test_01():
    input_file = templates['00_test']
    variables = input_file.variables
    old_string = str(input_file._string)
    print old_string
    input_file.update(Variable('K', 15.0))
    print input_file._string
    assert input_file._string != old_string

if __name__ == '__main__':
    test_01()