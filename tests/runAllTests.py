# -*- coding: utf-8 -*-
import imp
import sys
import os.path
import unittest


prPath = os.path.dirname(sys.modules['__main__'].__file__)

alltests = unittest.TestSuite()
alltests.addTest(unittest.defaultTestLoader.discover(prPath))

runner = unittest.TextTestRunner()
runner.run(alltests)
