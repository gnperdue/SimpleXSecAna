#!/usr/bin/env python
from __future__ import print_function
import os


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

    data_path = '/genie/app/perdue/GENIE/lamp_svn_286/lamp/genie_runs'

    # try to match: any name + <4 digit number> + '.ghep.root'
    search_string = r'(.*)[0-9]{4}\.ghep\.root'

    files = enumerate_files(data_path, search_string)
    write_list_of_files_to_file(files, 'ghep_file_list.txt')
