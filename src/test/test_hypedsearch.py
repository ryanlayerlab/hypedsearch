import unittest

class TestSum(unittest.TestCase):

    def test_add_pass(self):
        self.assertEqual(2, 2, "Should pass")
    
if __name__ == '__main__':
    unittest.main()
