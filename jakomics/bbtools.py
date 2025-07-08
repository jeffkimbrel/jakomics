import sys

import pandas as pd
from jakomics.utilities import system_call, check_executable
from jakomics import colors

def bbtools_version(echo=False, run=True, return_type='err'):
    if check_executable("bbmap.sh") == False:
        sys.exit(f"{colors.bcolors.RED}Error: BBTools not found!{colors.bcolors.END}")
    else :
        lines = system_call("bbmap.sh --version", echo=echo, run=run, return_type = return_type)[4]
        print(f"{colors.bcolors.GREEN}{lines}{colors.bcolors.END}")

def bbduk_stats_parser(stats_file_path):
    stats = pd.read_csv(stats_file_path,
                        sep=None, names=['TYPE', 'READS', 'PERCENT'],
                        header=None,
                        engine='python')
    stats = stats.set_index('TYPE')

    total = stats.loc['#Total', 'READS']
    matched = stats.loc['#Matched', 'READS']

    return(total, matched)


def extract_stats(lines):
    stats = {}

    for line in lines:
        if line.startswith("Input:"):
            stats["IN"] = {"READS": line.split()[1],  "BP": line.split()[3]}
        elif line.startswith("FTrimmed:"):
            stats["FTrimmed"] = {"READS": line.split()[1],  "BP": line.split()[4]}
        elif line.startswith("KTrimmed:"):
            stats["KTrimmed"] = {"READS": line.split()[1],  "BP": line.split()[4]}
        elif line.startswith("Trimmed by overlap::"):
            stats["Overlap_Trimmed"] = {"READS": line.split()[3],  "BP": line.split()[6]}
        elif line.startswith("Contaminants:"):
            stats["CF"] = {"READS": line.split()[1],  "BP": line.split()[4]}
        elif line.startswith("QTrimmed:"):
            stats["QT"] = {"READS": line.split()[1],  "BP": line.split()[4]}
        elif line.startswith("Total Removed:"):
            stats["REMOVED"] = {"READS": line.split()[2],  "BP": line.split()[5]}
        elif line.startswith("Result:"):
            stats["OUT"] = {"READS": line.split()[1],  "BP": line.split()[4]}

    return stats

if __name__ == "__main__":

    print(bbtools_version())
