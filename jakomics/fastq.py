import subprocess
from jakomics.file import FILE


class FASTQ(FILE):

    def __str__(self):
        return "<JAKomics FASTQ class>"


if __name__ == "__main__":
    import argparse
    from jakomics import utilities

    parser = argparse.ArgumentParser(
        description='Check the Illumina run information on a directory (--in_dir) or list (-f) of fastq.gz files')

    parser.add_argument('--in_dir',
                        help="Directory with fastq.gz files",
                        required=False,
                        default=None)

    args = parser.parse_args()

    file_list = utilities.get_files([], args.in_dir, ["fastq.gz", "fastq"], file_class="FASTQ")

    for file in file_list:
        print(file)
