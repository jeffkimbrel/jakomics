import pandas as pd


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
