#!/usr/bin/env python

import argparse
import os, subprocess, re


def generate_cmdhelp(cmds, outfolder):
    """Generate docs from commands, then store them to the output foder.
    Output the description of the help doc.
    """
    cmd_descriptions = []
    for cmd in cmds:
        if cmd:
            cmd_help = subprocess.check_output((cmd, '-h'))

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


def generate_rst_for_cmds(cmds, cmd_descriptions, outfolder):
    """Generate rst files.
    """
    for cmd in cmds:
        rstfilename = cmd + '.rst'
        rstfile = open(os.path.join(outfolder, rstfilename), 'w')
        rstfile.write(cmd + "\n"+len(cmd)*"="+"\n\n" +
            ".. literalinclude:: " + cmd + ".txt")
        rstfile.close()

    totalrstfilename = 'commandlist.rst'
    rstfile = open(os.path.join(outfolder, totalrstfilename), 'w')
    rstfile.write(r"""
============
Command List
============

.. csv-table::
   :header: Command, Description
   :widths: 20, 20

""")

    for cmd, description in zip(cmds, cmd_descriptions):
        rstfile.write('   `' + cmd + ' <' + cmd + '.html>`__, "' + description + '"\n')

    rstfile.write(r"""
.. toctree::
   :maxdepth: 1
   :hidden:

""")

    for cmd in cmds:
        rstfile.write('   '+cmd+'\n')

    rstfile.close()


def main():

    parser = argparse.ArgumentParser(description='''\
                                     Generate documents for a list of commands,
                                     then generate rst file for sphinx.
                                     ''')
    parser.add_argument('list', nargs=1, help='a file with a list of commands')
    parser.add_argument('output', nargs=1, help='output folder')

    args = parser.parse_args()
    cmdlistfile = args.list[0]
    outfolder = args.output[0]

    # mkdir output folder
    subprocess.call(('mkdir', '-p', outfolder))

    # read commands
    lines = [line.rstrip('\n') for line in open(cmdlistfile)]
    cmds = sorted(lines)
    cmds = [x for x in cmds if x]

    # generate help doc files
    cmd_descriptions = generate_cmdhelp(cmds, outfolder)

    # generate rst files
    generate_rst_for_cmds(cmds, cmd_descriptions, outfolder)


if __name__ == '__main__':
    main()


