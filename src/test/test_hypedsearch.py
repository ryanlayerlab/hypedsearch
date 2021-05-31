#import unittest

#Fails
#import hypedsearch
#hypedsearch = __import__("../hypedsearch/hypedsearch.py")
#hypedsearch = __import__("hypedsearch.py")


import unittest
import sys    

class Test_Main(unittest.TestCase):

    def setUp(self):
        sys.path.insert(0, "../hypedsearch")
        import hypedsearch
        #from project import item
        # further setup using this import

    def test_add_pass(self):
        self.assertEqual(2, 2, "Should pass")

if __name__ == "__main__":
    unittest.main()

