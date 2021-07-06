import sys  
import os.path  
import unittest
import hypedsearch
from hypedsearch import runner

class Test_Main(unittest.TestCase):

    main_arguments = None
    
    def setUp(self):
        sys.path.insert(0, "../hypedsearch")
        
    def test_main_arguments_specta_folder(self):
        self.assertEqual(1,1,'equality test')
        
    def test_runner(self):
        self.assertEqual(1,1,'hello world')

if __name__ == "__main__":
    unittest.main()
