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
    stats = {"input": 0, "lowQC": 0, "contamination": 0, "removed": 0, "remain": 0}

    for line in lines:
        if line.startswith("Input:"):
            stats["input"] = line.split()[1]
        elif line.startswith("Contaminants:"):
            stats["contamination"] = line.split()[1]
        elif line.startswith("Low quality discards:"):
            stats["lowQC"] = line.split()[3]
        elif line.startswith("Total Removed:"):
            stats["removed"] = line.split()[2]
        elif line.startswith("Result:"):
            stats["remain"] = line.split()[1]

    return(stats)
