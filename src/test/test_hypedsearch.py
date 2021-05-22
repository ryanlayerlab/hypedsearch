import unittest

class TestSum(unittest.TestCase):

    def test_add_pass(self):
        self.assertEqual(2, 2, "Should pass")

    def test_add_fail(self):
        self.assertEqual(1, 2, "Should fail")
    
if __name__ == '__main__':
    unittest.main()
