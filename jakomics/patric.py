import json
import os
from natsort import natsorted
import operator


class Patric():

    def __str__(self):
        return "<JAKomics Patric Class>"

    def __init__(self, main_folder, full_name):
        self.full_name = full_name
        self.features_path = os.path.join(main_folder, full_name, full_name + ".txt")
        self.metadata_path = os.path.join(main_folder, full_name, "load_files", "genome.json")
        self.subsystems_path = os.path.join(main_folder, full_name, "load_files", "subsystem.json")
        self.apply_metadata()

    def apply_metadata(self):
        metadata = json.loads(open(self.metadata_path, "r").read())
        self.genome_name = metadata[0]['genome_name']
        self.genome_id = metadata[0]['genome_id']


def get_files(main_folder):
    genome_list = []
    folder_list = [d for d in os.listdir(
        main_folder) if os.path.isdir(os.path.join(main_folder, d))]

    for folder in folder_list:
        genome_list.append(Patric(main_folder, folder))

    sorted_by_name = natsorted(genome_list, key=operator.attrgetter('genome_name'))
    return sorted_by_name


class Subsystem():

    def __str__(self):
        return "<JAKomics Subsystem Class>"

    def __init__(self, json_object):
        for key in json_object:
            if key == 'class':
                setattr(self, "subsystem_class", json_object["class"])
            else:
                setattr(self, key, json_object[key])

        self.locus_tag = self.patric_id.replace("fig|" + self.genome_id, self.genome_name)

    def print_hierarchy(self, sep="\t"):
        print(self.genome_name, self.superclass, self.subsystem_class,
              self.subclass, self.subsystem_name, self.active, self.product, self.locus_tag, sep=sep)


class Patric_Gene():

    def __str__(self):
        return "<JAKomics Patric Gene Class>"

    def __init__(self, id):
        self.id = id

    def result(self):
        '''
        method to get results suitable for gator
        '''

        return {'gene': self.id,
                'annotation': self.annotation,
                'score': getattr(self, "score", None),
                'evalue': getattr(self, "evalue", None)}


if __name__ == '__main__':

    subsystem = json.loads(open(
        "/Users/kimbrel1/Dropbox/LLNL/Projects/Biofuels_SFA/CSP1939/isolates/shared/patric/20201223/PT6CLA/load_files/subsystem.json", "r").read())

    for a in subsystem:
        s = Subsystem(a)
        # print(s.patric_id)
        s.print_hierarchy()
