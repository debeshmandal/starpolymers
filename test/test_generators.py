from starpolymers.io.configuration import *
from starpolymers.io.input_file import *
import starpolymers.generators.system as system

def test_input_file():
    assert True

def test_system():
    system.System(100)

if __name__ == '__main__':
    test_system()