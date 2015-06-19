#!/usr/bin/env python
"""
Generate outputs for tutorils.
It may take a long time
"""
tutorial_list = ['tutorial_qspacesampling']

import os, subprocess, argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('sphinxDir', nargs=1, help='dmritool sphinx dir')
args = parser.parse_args()

sphinxDir = os.path.abspath(args.sphinxDir[0])
origWD = os.getcwd()


os.chdir(os.path.join(sphinxDir, '_build/html'))
for tutorial in tutorial_list:
    rstfile = os.path.join(sphinxDir, tutorial + '.rst')
    shfile = os.path.join(sphinxDir, '_build/html', tutorial + '.sh')
    subprocess.call(('python', os.path.join(sphinxDir, 'tools/extract_bash_from_rst.py'), rstfile,  shfile ))
    subprocess.call(['sh', shfile])

os.chdir(origWD)
