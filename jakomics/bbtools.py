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
    # stats = {"input": 0, "lowQC": 0, "contamination": 0, "removed": 0, "remain": 0}
    stats = {}

    for line in lines:
        if line.startswith("Input:"):
            stats["IN"] = {"READS": line.split()[1],  "BP": line.split()[3]}
        elif line.startswith("FTrimmed:"):
            stats["FTrimmed"] = {"READS": line.split()[1],  "BP": line.split()[4]}
            # stats["READS_FTrimmed"] = line.split()[1]
            # stats["BP_FTrimmed"] = line.split()[4]
        elif line.startswith("KTrimmed:"):
            stats["KTrimmed"] = {"READS": line.split()[1],  "BP": line.split()[4]}
            # stats["READS_KTrimmed"] = line.split()[1]
            # stats["BP_KTrimmed"] = line.split()[4]
        elif line.startswith("Trimmed by overlap::"):
            stats["Overlap_Trimmed"] = {"READS": line.split()[3],  "BP": line.split()[6]}
            # stats["READS_Overlap_Trimmed"] = line.split()[3]
            # stats["BP_Overlap_Trimmed"] = line.split()[6]
        elif line.startswith("Contaminants:"):
            stats["CF"] = {"READS": line.split()[1],  "BP": line.split()[4]}
            # stats["READS_CF"] = line.split()[1]
            # stats["BP_CF"] = line.split()[4]
        elif line.startswith("QTrimmed:"):
            stats["QT"] = {"READS": line.split()[1],  "BP": line.split()[4]}
            # stats["READS_QT"] = line.split()[1]
            # stats["BP_QT"] = line.split()[4]
        elif line.startswith("Total Removed:"):
            stats["REMOVED"] = {"READS": line.split()[2],  "BP": line.split()[5]}
            # stats["READS_REMOVED"] = line.split()[2]
            # stats["BP_REMOVED"] = line.split()[5]
        elif line.startswith("Result:"):
            stats["OUT"] = {"READS": line.split()[1],  "BP": line.split()[4]}
            # stats["READS_OUT"] = line.split()[1]
            # stats["BP_OUT"] = line.split()[4]

    return(stats)
