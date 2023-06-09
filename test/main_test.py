import sys 
sys.path.append("..")
import unittest
import os
import os.path as osp
from dl.main import main
from dl import args


class MyTestCase(unittest.TestCase):
    def test_main(self):
        main_script_dir = "../dl/"
        os.chdir(osp.dirname(main_script_dir))
        args.max_epochs = 1
        #main(args)

    def test_eval(self):
        main_script_dir = "../dl/"
        os.chdir(osp.dirname(main_script_dir))
        args.model = "pdglstm_bn"    # test model
        args.num_gpus = 1           
        args.model_num = 2       # test model number
        args.mode = "eval"        
        main(args)


if __name__ == '__main__':
    unittest.main()
