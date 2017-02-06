#!/usr/bin/env python
"""
Get a list of files. The script makes assumptions about the path and file
extension that are probably not what you want.

Usage:
    python walk_for_file_list.py /path/to/files filename-extension out-file

One may optionally supply the args, but note the expected position, etc.
"""
from __future__ import print_function
import os
import sys


def enumerate_files(path, keystring):
    '''
    Get a list of all files under `path` containing the string `keystring`.
    '''
    import re

    file_collection = []
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if re.match(keystring, filename):
                fullpath = os.path.join(dirpath, filename)
                file_collection.append(fullpath)

    return file_collection


def write_list_of_files_to_file(list_of_files, filenam):
    with open(filenam, 'w') as f:
        for file in list_of_files:
            print(file, file=f)


if __name__ == '__main__':

    if '-h' in sys.argv or '--help' in sys.argv:
        print(__doc__)
        sys.exit(1)

    data_path  = '/pnfs/minerva/scratch/users/minervapro'
    data_path += '/mc_production_genie_DFR_v10r8p4/grid/central_value/minerva'
    data_path += '/genie/v10r8p4/00/01/00/'

    if len(sys.argv) > 1:
        data_path = sys.argv[1]

    search_extension = r'(.*)v10r8p4_DFR_ghep.root'

    if len(sys.argv) > 2:
        search_extension = sys.argv[2]

    out_file = 'events_file_list.txt'

    if len(sys.argv) > 3:
        out_file = sys.argv[3]

    print('data path:', data_path)
    print('search extension:', search_extension)
    print('out file:', out_file)

    files = enumerate_files(data_path, search_extension)
    write_list_of_files_to_file(files, out_file)
