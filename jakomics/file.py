import os
import uuid
import subprocess


class FILE:

    def __init__(self, file_path):
        self.file_path = file_path
        self.name = os.path.basename(self.file_path)
        self.dir = os.path.dirname(file_path)
        self.short_name = os.path.splitext(self.name)[0]
        self.suffix = os.path.splitext(self.name)[1]
        self.id = uuid.uuid4().hex
        self.stats = os.stat(self.file_path)

    def __str__(self):
        return "<JAKomics FILE class>"

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
