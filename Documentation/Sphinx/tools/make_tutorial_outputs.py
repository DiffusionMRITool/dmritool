#!/usr/bin/env python
"""
Generate outputs for tutorils.
It may take a long time
"""

import os, subprocess, argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('sphinxDir', nargs=1, help='dmritool sphinx dir')
args = parser.parse_args()

sphinxDir = os.path.abspath(args.sphinxDir[0])
origWD = os.getcwd()

f = open(os.path.join(sphinxDir, 'tutoriallist.txt'),"r")
tutorial_list = f.readlines()
f.close()
tutorial_list = [t.strip() for t in tutorial_list if t.strip() ]

for tutorial in tutorial_list:

    rstfile = os.path.join(sphinxDir, tutorial + '.rst')
    if not os.path.isfile(rstfile):
        raise Exception('no file: ' + rstfile)

    runpath = os.path.join(sphinxDir, '.' + tutorial)
    if os.path.exists(runpath):
        subprocess.call(['rm', '-rf', runpath])
    os.makedirs(runpath)
    os.chdir(runpath)

    shfile = os.path.join(runpath, tutorial + '.sh')
    subprocess.call(('python', os.path.join(sphinxDir, 'tools/extract_bash_from_rst.py'), rstfile,  shfile ))
    subprocess.call(['sh', shfile])

    os.chdir(origWD)

