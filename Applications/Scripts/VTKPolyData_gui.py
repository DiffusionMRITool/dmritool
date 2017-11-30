#!/usr/bin/env python

'''
Description: GUI for VTKPolyData.py.
'''

import sys

import utlDMRITool as utl
import VTKPolyData

from gooey import Gooey

@Gooey(program_name='VTKPolyData @ dmritool',default_size=(1000, 800))
def main():
    VTKPolyData.main()

if __name__ == '__main__':
    if '-h' in sys.argv or '--help' in sys.argv:
        print(utl.app_doc(__doc__))
    else:
        main()
