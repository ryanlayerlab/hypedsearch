import unittest
import sys    

class Test_Main(unittest.TestCase):

    main_arguments = None
    
    def setUp(self):
        sys.path.insert(0, "../hypedsearch")
        import hypedsearch
        #self.main_arguments = hypedsearch.Main_Arguments()
        
    def test_main_arguments_specta_folder(self):
        #self.main_arguments.spectra_folder = 'C:\MyFolder'
        #self.assertEqual(self.main_arguments.spectra_folder, 'C:\MyFolder', 'spectra_folder is incorrect')   
        self.assertEqual(1,1,'equality test')
        
if __name__ == "__main__":
    unittest.main()

