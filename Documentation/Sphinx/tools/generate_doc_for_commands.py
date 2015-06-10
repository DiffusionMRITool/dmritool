#!/usr/bin/env python

import argparse
import os, subprocess, re, glob


def which(exe):
    for dir in os.getenv("PATH").split(':'):
        if (os.path.exists(os.path.join(dir, exe))):
            return os.path.join(dir, exe)
    return None

def is_exe(exe):
    return which(exe)!=None

def get_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def generate_cmd_lists(appfolder):
    """Generate cmd lists from application folder."""

    # get subfolders in appfolder
    app_categories = get_subdirectories(appfolder)
    app_categories = sorted(app_categories)

    # get app_list for each subfolder
    apps = []
    for app_category in app_categories:
        app_list = glob.glob(os.path.join(appfolder, app_category, '*.cxx'))
        app_list = [os.path.splitext(os.path.basename(name))[0] for name in app_list]

        # add vtkviewer into Visualization
        if app_category=='Visualization':
            app_list.append('vtkviewer')

        app_list = [name for name in app_list if is_exe(name)]
        app_list = sorted(app_list)
        apps.append(app_list)

    return app_categories, apps


def generate_cmd_helpfiles(cmds, outfolder):
    """Generate docs from commands, then store them to the output foder.
    Output the description of the help doc.
    """
    cmd_descriptions = []
    for cmd in cmds:
        cmd_help = subprocess.check_output((cmd, '-h'))
        # cmd_help = '\nDescription: adadad\n'

        outfile = open(os.path.join(outfolder, cmd+'.txt'), 'w')
        outfile.write(cmd + ' -h\n\n')
        outfile.write(cmd_help)
        outfile.close()

        cmd_help_lines = cmd_help.split('\n');
        description_found = False
        for line in cmd_help_lines:
            if not line and not description_found:
                continue
            if re.search(r"Description:", line):
                description_found = True
                cmd_help_desciption = line[line.find('Description:')+13:]
            else:
                if description_found:
                    if line:
                        cmd_help_desciption = cmd_help_desciption + line
                    else:
                        cmd_descriptions.append(cmd_help_desciption)
                        break

    return cmd_descriptions



def main():

    parser = argparse.ArgumentParser(description='''\
                                     Generate documents for executables in application folder,
                                     then generate rst file for sphinx.
                                     ''')
    parser.add_argument('output', nargs=1, help='output folder')

    args = parser.parse_args()
    outfolder = args.output[0]

    # mkdir output folder
    subprocess.call(('mkdir', '-p', outfolder))

    (app_categories, apps) = generate_cmd_lists('../../Applications')

    rstfile_name = 'commandlist.rst'
    rstfile = open(os.path.join(outfolder, rstfile_name), 'w')
    rstfile.write(r"""
============
Command List
============

.. contents:: Table of Contents
   :depth: 1
   :local:

""")


    rstfile.write('\n\n')

    for app_category, cmds in zip(app_categories, apps):

        rstfile.write(''.join([app_category, "\n", len(app_category)*"=", "\n\n"]))
        rstfile.write(r"""
.. toctree::
   :maxdepth: 1
   :hidden:

""")

        for cmd in cmds:
            rstfile.write(''.join(['   ', cmd, '\n']))

        rstfile.write('\n\n')

        rstfile.write(r"""
.. csv-table::
   :header: Command, Description
   :widths: 20, 20

""")

        cmd_descriptions = generate_cmd_helpfiles(cmds, outfolder)
        # print cmds
        # print cmd_descriptions

        for cmd in cmds:
            rstfile_cmd_name = cmd + '.rst'
            rstfile_cmd = open(os.path.join(outfolder, rstfile_cmd_name), 'w')
            rstfile_cmd.write(''.join([cmd, "\n", len(cmd)*"=", "\n\n",
                          ".. literalinclude:: ", cmd, ".txt"]))
            rstfile_cmd.close()

        for cmd, description in zip(cmds, cmd_descriptions):
            rstfile.write(''.join(['   `', cmd, ' <',cmd, '.html>`__, "', description,'"\n']))

        rstfile.write('\n\n')


    rstfile.close()


if __name__ == '__main__':
    main()


