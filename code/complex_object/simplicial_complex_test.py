import unittest
from simplicial_complex import ChainComplex

# to test, run `$ python3 -m unittest simplicial_complex_test.py` in terminal

class TestComplexMethods(unittest.TestCase):

    def test_simplex_addition(self, test_object=ChainComplex()):
        test_object.add_simplex([0,1])
        self.assertTrue({0,1} in test_object[2].keys())

    def test_simplex_deletion(self, test_object=ChainComplex([0,1])):
        test_object.delete_simplex([0,1])
        self.assertTrue({0,1} not in test_object[2].keys())


if __name__ == '__main__':
    unittest.main()
