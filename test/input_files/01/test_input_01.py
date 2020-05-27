from starpolymers.generators.input_file import templates, Variable
import pandas as pd

def test_01():
    input_file = templates['01_master']
    variables = input_file.variables
    print pd.Series(variables).apply(lambda x: x.string)

    input_file.update(Variable('k1', 10.0))
    print pd.Series(variables).apply(lambda x: x.string)

    return

if __name__ == '__main__':
    test_01()