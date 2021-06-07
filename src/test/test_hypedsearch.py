import unittest
import sys    

class Test_Main(unittest.TestCase):

    main_arguments = None
    
    def setUp(self):
        sys.path.insert(0, "../hypedsearch")
        import hypedsearch
        args = [use_config_file=False,another_arg=True]
        hypedsearch.main(args)
        
    def test_main_arguments_specta_folder(self):
        self.assertEqual(1,1,'equality test')
        
if __name__ == "__main__":
    unittest.main()

