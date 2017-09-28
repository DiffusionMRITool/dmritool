#!/usr/bin/env python

"""
Generate documents for executables in application folder,
then generate rst file for sphinx.
"""

import argparse
import os, subprocess, re, glob

apps_added={
    'Visualization' :
        ['vtkviewer', 'VTKPolyData.py', 'VTKPolyData_gui.py', 'CombineVTKPolyData.py']
}

def which(exe):
    for dir in os.getenv("PATH").split(':'):
        if (os.path.exists(os.path.join(dir, exe))):
            return os.path.join(dir, exe)
    return None

def is_exe(exe):
    exe_real = which(exe)
    return (exe_real is not None) and (os.access(exe_real, os.X_OK))

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
        if app_category=='Scripts':
            app_list = glob.glob(os.path.join(appfolder, app_category, '*.py'))
            app_list = [os.path.basename(name) for name in app_list]
        else:
            app_list = [os.path.splitext(os.path.basename(name))[0] for name in app_list]

        # add vtkviewer into Visualization
        if app_category in apps_added:
            for app_tmp in apps_added[app_category]:
                app_list.append(app_tmp)

        app_list = [name for name in app_list if is_exe(name)]
        app_list = sorted(app_list)
        apps.append(app_list)

    return app_categories, apps


def generate_cmd_helpfiles(cmds, outfolder):
    """Generate docs from commands, then store them to the output foder.
    Output the description of the help doc.
    """
    cmd_descriptions = []
    recheck = re.compile(r"Description:")
    for cmd in cmds:
        cmd_help = subprocess.check_output((cmd, '-h'))

        rstfile_cmd_name = cmd + '.rst'
        rstfile_cmd = open(os.path.join(outfolder, rstfile_cmd_name), 'w')
        rstfile_cmd.write(''.join([cmd, "\n", len(cmd)*"=", "\n\n"]))
        rstfile_cmd.write('.. code-block:: none\n\n')

        cmd_help_lines = cmd_help.split('\n')
        rstfile_cmd.write('\n'.join(['  ' + line_i for line_i in cmd_help_lines]))

        rstfile_cmd.close()

        description_found = False
        for line in cmd_help_lines:
            if not line and not description_found:
                continue
            if recheck.search(line):
                description_found = True
                cmd_help_desciption = line[line.find('Description:')+13:]
            else:
                if description_found:
                    if line:
                        cmd_help_desciption = ' '.join([cmd_help_desciption, line])
                    else:
                        cmd_descriptions.append(cmd_help_desciption)
                        break

    return cmd_descriptions



def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs=1, help='input cmd source folder')
    parser.add_argument('output', nargs=1, help='output cmd folder')

    args = parser.parse_args()
    infolder = args.input[0]
    outfolder = args.output[0]

    # mkdir output folder
    subprocess.call(('mkdir', '-p', outfolder))

    (app_categories, apps) = generate_cmd_lists(infolder)

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


        rstfile.write(''.join(['.. _commandlist_', app_category, ':' ]))
        rstfile.write('\n\n')
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

        for cmd, description in zip(cmds, cmd_descriptions):
            rstfile.write(''.join(['   :doc:`', cmd, ' <',cmd, '>`, "', description,'"\n']))

        rstfile.write('\n\n')


    rstfile.close()


if __name__ == '__main__':
    main()


