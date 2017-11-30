#!/usr/bin/env python
"""
Description: Combine several vtk files into one vtk file.

Usage:
  CombineVTKPolyData.py  <vtkfiles>... -o <combined> [-v]
  CombineVTKPolyData.py (-h | --help)

Options:
  -h --help                Show this screen.
  -v --verbose             Verbose
  -o --output <combined>   Output file

Examples:
CombineVTKPolyData.py file1.vtk file2.vtk -o file_combine.vtk
CombineVTKPolyData.py file1.vtk -o file_combine.vtp

Author(s): Jian Cheng (jian.cheng.1983@gmail.com)
"""

import vtk
from docopt import docopt

import utlVTK
import utlDMRITool as utl


def main():

    args = docopt(utl.app_doc(__doc__), version='1.0')

    if (args['--verbose']):
        print(args)

    append = vtk.vtkAppendPolyData()

    for inputFile in args['<vtkfiles>']:
        polyData = utlVTK.readPolydata(inputFile)
        append.AddInputData(polyData)

    append.Update()
    utlVTK.savePolydata(append.GetOutput(), args['--output'])


if __name__ == '__main__':
    main()
