import unittest
import hypedsearch

class TestSum(unittest.TestCase):

    def test_add_pass():
        self.assertEqual(1+1, 2, "Should be 2")

    def test_add_fail():
        self.assertEqual(1+1, 1, "Should be 2")
    
if __name__ == '__main__':
    unittest.main()
