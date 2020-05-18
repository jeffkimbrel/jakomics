import os
import uuid

class FILE:

    def __init__(self, file):
        self.file_path = os.path.abspath(file)
        self.name = os.path.basename(self.file_path)
        self.id = uuid.uuid4().hex
        self.stats = os.stat(self.file_path)
