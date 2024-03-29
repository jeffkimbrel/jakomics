import os
import sys
import shutil 

from natsort import natsorted
import uuid
import hashlib
from jakomics.file import FILE
import operator
import subprocess
from jakomics import colors


def test():
    print("utilities module loaded correctly")


def get_unique_ID():
    return uuid.uuid4().hex

def string_to_hash(s):
    return hashlib.md5(s.encode('UTF-8')).hexdigest()

def get_file_list(input_file_list, ending, file_list=[]):
    '''
    input_file_list: list of files
    file_list: an existing list to append to (empty by default)
    '''
    for file_name in input_file_list:
        if file_name.endswith(tuple(ending)):
            abs_path = os.path.abspath(file_name)
            file_list.append(abs_path)

    return natsorted(list(set(file_list)))


def get_directory_file_list(directory, ending, file_list=[]):
    '''
    ending: a string, list, or tuple. All will get converted to a tuple
    directory: path to search
    existing_list: existing list to append to (empty by default)
    '''
    if directory != "":
        directory = os.path.abspath(directory) + '/'
        directory_list = os.listdir(directory)
        for file_name in directory_list:
            if file_name.endswith(tuple(ending)):
                file_list.append(os.path.join(directory, file_name))

    return natsorted(list(set(file_list)))


def get_files(file_names, directory, file_type=""):
    files = {}

    for file in file_names:
        file_obj = FILE(os.path.abspath(file))
        files[file_obj.file_path] = file_obj

    if directory != "":
        directory_list = os.listdir(os.path.abspath(directory) + '/')
        for file in directory_list:
            if file.endswith(tuple(file_type)):
                file_obj = FILE(os.path.abspath(directory) + '/' + file)
                files[file_obj.file_path] = file_obj

    # return list of values, should be unique by file path
    sorted_by_path = natsorted(list(files.values()), key=operator.attrgetter('file_path'))

    if len(sorted_by_path) == 0:
        sys.exit(f"{colors.bcolors.RED}Error: No valid files were found!{colors.bcolors.END}")

    return sorted_by_path

def check_executable(app):
    return shutil.which(app) != None


def system_call(call, echo=False, run=True, return_type='err'):
    if echo:
        print(call, file=sys.stderr)

    if run:
        p1 = subprocess.Popen(call,
                              shell=True,
                              stdin=None,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        out, err = p1.communicate()
        if return_type == 'err':
            err = err.decode()
            lines = err.split('\n')
            return lines
        elif return_type == 'out':
            out = out.decode()
            lines = out.split('\n')
            return lines
        else:
            return [out, err]

    else:
        return []

if __name__ == "__main__":
    print(check_executable(sys.argv[1]))