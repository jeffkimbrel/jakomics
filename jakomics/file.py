import os
import uuid
import subprocess
import sys
from jakomics import colors


class FILE:

    def __str__(self):
        return "<JAKomics FILE class>"

    def __init__(self, file_path):
        validate_path(file_path)
        self.file_path = file_path
        self.name = os.path.basename(self.file_path)
        self.dir = os.path.abspath(os.path.dirname(file_path))
        self.short_name = os.path.splitext(self.name)[0]
        self.suffix = os.path.splitext(self.name)[1]
        self.id = uuid.uuid4().hex
        self.stats = os.stat(self.file_path)
        self.temp_files = {}

    def check_files_exist(self, exit_if_false = True):
        if not os.path.exists(self.file_path):
            if exit_if_false:
                raise FileNotFoundError
            else:
                return(False)
        else:
            return(True)

    def remove_temp(self):
        for temp_file in self.temp_files:
            os.remove(self.temp_files[temp_file])

    def view(self):
        print(self.short_name, self.name, self.file_path, self.id, sep="\t")

    def get_md5(self):
        call = "md5 " + self.file_path
        p1 = subprocess.Popen(call, shell=True,
                              stdin=None,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        out, err = p1.communicate()
        out = out.decode()
        junk, md5 = out.split(" = ")
        self.md5 = md5.strip()


def validate_path(path):
    try:
        f = open(path)
        f.close()
    except IOError:
        sys.exit(f'{colors.bcolors.RED}ERROR - File not found: {path}{colors.bcolors.END}')
