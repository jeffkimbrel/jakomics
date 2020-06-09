import os
import uuid


class FILE:

    def __init__(self, file_path):
        self.file_path = file_path
        self.name = os.path.basename(self.file_path)
        self.short_name = os.path.splitext(self.name)[0]
        self.id = uuid.uuid4().hex
        self.stats = os.stat(self.file_path)
