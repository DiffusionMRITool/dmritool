
import os
import subprocess


def app_doc(strIn):
    '''add copyright information into doc'''

    str_copy='''Copyright (c) 2015-2017 the dmritool contributors.
For more information, see https://diffusionmritool.github.io'''

    return '\n'.join([strIn, str_copy])


def get_file_extension(file):
    """Get file path and extension.
    GetFileExtension("aa/bb.nii.gz") ->  ('aa/bb', '.nii.gz')
    GetFileExtension("aa/bb.nii") ->  ('aa/bb', '.nii')

    Parameters:
        file : filename string

    Returns:
        fileExt:   extension string
        filePath:  file path string
    """

    [filePath, fileExt] = os.path.splitext(file)
    (filePath2, fileExt2) = os.path.splitext(filePath)

    if fileExt2 == "":
        return filePath, fileExt
    else:
        return filePath2, fileExt2 + fileExt


def process(cmd, print_cmd=True, print_sto=True, print_ste=True):
    '''Wrap subprocess.Popen to process a shell command

    cmd :   string,
              command line to be run

    print_sto : boolean
         Print standard output (stdout) or not (default: True)

    print_ste : boolean
         Print standard error (stderr) or not (default: True)
    '''

    if print_cmd:
        print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sto = p.stdout.readlines()
    ste = p.stderr.readlines()
    outCode = p.wait()
    if print_sto and sto:
        print(''.join(sto))
    if print_ste and ste:
        print(''.join(ste))
    if outCode!=0:
        raise Exception("Fail: " + cmd + "!\n")


def remove_empty_lines(strIn):
    '''Remove empty lines from a string'''
    return os.linesep.join([s for s in strIn.splitlines() if s.strip()])


