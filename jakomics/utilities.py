import os
from natsort import natsorted
import uuid


def test():
    print("utilities module loaded correctly")


def get_unique_ID():
    return uuid.uuid4().hex


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

    directory = os.path.abspath(directory) + '/'
    directory_list = os.listdir(directory)
    for file_name in directory_list:
        if file_name.endswith(tuple(ending)):
            file_list.append(os.path.join(directory, file_name))

    return natsorted(list(set(file_list)))
