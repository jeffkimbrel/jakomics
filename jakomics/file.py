import os
import uuid


class FILE:

    def __init__(self, file_path):
        self.file_path = file_path
        self.name = os.path.basename(self.file_path)
        self.short_name = os.path.splitext(self.name)[0]
        self.suffix = os.path.splitext(self.name)[1]
        self.id = uuid.uuid4().hex
        self.stats = os.stat(self.file_path)

    def __str__(self):
        return "<JAKomics GENE class>"

    def view(self):
        print(self.short_name, self.name, self.file_path, self.id, sep="\t")
