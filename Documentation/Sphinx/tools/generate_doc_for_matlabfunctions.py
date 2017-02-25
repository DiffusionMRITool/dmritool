#!/usr/bin/env python

"""
Generate documents for matlab functions,
then generate rst file for sphinx.
"""

import argparse
import os, subprocess, re, glob

def get_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def get_filename(file):
    return os.path.splitext(os.path.basename(file))[0]

def get_helpdoc_from_matlabfile(file):
    """Get help doc and brief descriptions from a matlab file"""
    f = open(file,"r")
    lines = f.readlines()
    f.close()

    doc = ''
    description = ''
    firstline = False
    descriptioncheck = re.compile(r'^[% ]*\n$')
    descriptionEnd = False
    for line in lines:
        if line[0]=='%':
            firstline = True
            doc += line;
            if not descriptionEnd and descriptioncheck.match(line):
                descriptionEnd = True
            else:
                if not descriptionEnd:
                    description += line
        else:
            if firstline:
                break
            else:
                continue

    description = description.replace('%', '')
    description = description.replace('\n', ' ')
    description = description.replace(get_filename(file) +':', '')
    description = description.strip()
    return doc, description

def write_helpdoc(doc, file, outputDir):
    """write doc into a file"""

    outfile = get_filename(file)

    f = open(os.path.join(outputDir, outfile+'.rst'), 'w')
    f.write(''.join([outfile, "\n", len(outfile)*"=", "\n\n"]))
    f.write('.. code-block:: none\n\n')
    cmd_help_lines = doc.split('\n')
    f.write('\n'.join(['  ' + line_i for line_i in cmd_help_lines]))
    f.close()

def generate_matlabfile_dictionary(matlabFolder):
    """Generate matlabfile lists from DMRITOOL/Matlab folder."""

    # get subfolders in appfolder
    file_categories = get_subdirectories(matlabFolder)
    file_dict = {folder:[] for folder in file_categories if folder!='Demos'}

    # get file_list for each subfolder
    for file_category in file_dict:
        file_list = glob.glob(os.path.join(matlabFolder, file_category, '*.m'))
        sorted(file_list, key=get_filename)
        file_dict[file_category]= file_list

    # remove Tests folder and folders without m files
    file_dict = dict((k, v) for k, v in file_dict.iteritems() if v and k!='Tests')
    return file_dict

def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('srcDir', nargs=1, help='dmritool source dir')
    parser.add_argument('buildDir', nargs=1, help='dmritool build dir')
    parser.add_argument('outputDir', nargs=1, help='output cmd folder')


    args = parser.parse_args()
    srcDir = args.srcDir[0]
    buildDir = args.buildDir[0]
    outputDir = args.outputDir[0]

    subprocess.call(('mkdir', '-p', outputDir))

    # dictionary of matlab files
    file_dict = generate_matlabfile_dictionary(srcDir + '/Matlab')
    # dictionary of matlab mex files
    filemex_dict = generate_matlabfile_dictionary(srcDir + '/Wrapping/Matlab')

    # replace original mex files using the built files
    for file_category, file_list in filemex_dict.iteritems():
        filebuild_list = []
        for file in file_list:
            filename = get_filename(file)
            filebuild = ''.join([buildDir, '/Wrapping/Matlab/bin/', filename, '.m'])
            filebuild_list.append(filebuild)
        filemex_dict[file_category] = filebuild_list

    # merge two dictionaries
    for file_category, file_list in filemex_dict.iteritems():
        if file_category in file_dict:
            for file in file_list:
                file_dict[file_category].append(file)
            file_dict[file_category] = sorted(file_dict[file_category], key=get_filename)
        else:
            file_dict[file_category] = file_list

    sorted(file_dict)


    rstfile_name = 'matlabfunctions.rst'
    rstfile = open(os.path.join(outputDir, rstfile_name), 'w')
    rstfile.write(r"""
====================
Matlab Function List
====================

.. contents:: Table of Contents
   :depth: 1
   :local:

""")


    rstfile.write('\n\n')

    file_categories = file_dict.keys()
    file_categories = sorted(file_categories)

    for file_category in file_categories:
        file_list = file_dict[file_category]

        rstfile.write(''.join(['.. _matlabfunctions_', file_category, ':' ]))
        rstfile.write('\n\n')
        rstfile.write(''.join([file_category, "\n", len(file_category)*"=", "\n\n"]))
        rstfile.write(r"""
.. toctree::
   :maxdepth: 1
   :hidden:

""")

        for file in file_list:

            rstfile.write(''.join(['   ', get_filename(file), '\n']))



        rstfile.write('\n\n')

        rstfile.write(r"""
.. csv-table::
   :header: Command, Description
   :widths: 20, 20

""")
        for file in file_list:
            doc, description = get_helpdoc_from_matlabfile(file)
            write_helpdoc(doc, file, outputDir)

            filename = get_filename(file)
            rstfile.write(''.join(['   :doc:`', filename, ' <',filename, '>`, "', description,'"\n']))

        rstfile.write('\n\n')

    rstfile.close()


if __name__ == '__main__':
    main()
