from jakomics import utilities, blast, hmm, gene, colors
from bin.genomes import gbk_to_fasta
import argparse
import os
import re
import sys
import pandas as pd
from natsort import natsorted
from Bio import SearchIO
from collections import Counter
import math

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa genomes",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--out_dir',
                    help="Directory to write results to",
                    required=True)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

# TASmania_HMMs_path = "/Users/kimbrel1/Science/repos/jakomics/db/TASmania_HMMs/"
# TASmania_metadata_A = pd.read_excel(os.path.join(
#     TASmania_HMMs_path, "journal.pcbi.1006946.s008.xlsx"), sheet_name='antitoxin_profile_annotation')
# TASmania_metadata_T = pd.read_excel(os.path.join(
#     TASmania_HMMs_path, "journal.pcbi.1006946.s008.xlsx"), sheet_name='toxin_profile_annotation')

tadb_data = pd.read_excel("/Users/kimbrel1/Science/repos/jakomics/db/TADB_2.0/TADB2.xlsx", sheet_name = "merged")
tadb_type = pd.read_excel("/Users/kimbrel1/Science/repos/jakomics/db/TADB_2.0/TADB2.xlsx", sheet_name = "type")

tadb_type = tadb_type.set_index('TA_ID')

# CLASSES #####################################################################

class PAIR():

    def __init__(self, gene1, gene2, id):
        self.id = id
        self.gene1 = gene1
        self.gene2 = gene2

    def get_replicon(self):
        return self.gene1.replicon

    def get_range(self):
        smallest = min(self.gene1.start, self.gene1.stop, self.gene2.start, self.gene2.stop)
        largest = max(self.gene1.start, self.gene1.stop, self.gene2.start, self.gene2.stop)

        return(str(smallest) + "-" + str(largest))

    def get_type(self):
        gene1_type = max(set(self.gene1.tadb_type), key = self.gene1.tadb_type.count)
        gene2_type = max(set(self.gene2.tadb_type), key = self.gene2.tadb_type.count)

        if gene1_type == gene2_type:
            return gene1_type
        else:
            return [gene1_type, gene2_type]

    def all_toxin_tabd_ids(self):
        if self.toxin == self.gene1.id:
            return self.gene1.tadb_ids
        elif self.toxin == self.gene2.id:
            return self.gene2.tadb_ids

    def all_antitoxin_tabd_ids(self):
        if self.antitoxin == self.gene1.id:
            return self.gene1.tadb_ids
        elif self.antitoxin == self.gene2.id:
            return self.gene2.tadb_ids

    def get_all_toxin_names(self):
        if self.toxin == self.gene1.id:
            return dict(Counter(self.gene1.tadb_gene_name))
        elif self.toxin == self.gene2.id:
            return dict(Counter(self.gene2.tadb_gene_name))

    def get_all_antitoxin_names(self):
        if self.antitoxin == self.gene1.id:
            return dict(Counter(self.gene1.tadb_gene_name))
        elif self.antitoxin == self.gene2.id:
            return dict(Counter(self.gene2.tadb_gene_name))

    def view(self, genome_name):
        print(genome_name,
              self.id,
              self.get_replicon(),
              self.get_range(),
              self.get_type(),
              self.toxin,
              self.toxin_name,
              self.antitoxin,
              self.antitoxin_name,
              self.shared_tadb_families,
              self.shared_tadb_ids,
              self.score,
              scores[self.score],
              self.all_toxin_tabd_ids(),
              self.all_antitoxin_tabd_ids(),
              self.get_all_toxin_names(),
              self.get_all_antitoxin_names(),
              sep = "\t"
        )

# FUNCTIONS ###################################################################

def gene_distance(gene1, gene2):

    '''
    Removed same strand requirements... some type I systems that are protein and RNA can be on different strands
    '''

    if gene1.id >= gene2.id: # only compare once
        return math.inf
    elif gene1.replicon != gene2.replicon:
        return math.inf
    # elif gene1.strand != gene2.strand:
    #     return math.inf
    else:
        shortest_distance = min(
            abs(gene1.start - gene2.start),
            abs(gene1.stop - gene2.stop),
            abs(gene1.stop - gene2.start),
            abs(gene1.start - gene2.stop)
        )
        return shortest_distance

def get_tadb_family(tadb):
    tadb = tadb.replace("TADB|", "")
    role, ta_id, _ = re.split('(\d+)', tadb)
    df = tadb_data.loc[tadb_data['TA_ID'] == int(ta_id)]
    family = df['TA_FAMILY'].tolist()
    prediction = df[role].tolist()

    return family, prediction, role, ta_id

def score_tadb_pairs(pairs):
    print(f'{colors.bcolors.PURPLE}Scoring Potential Pairs{colors.bcolors.END}')
    processed_pairs = []

    for pair in pairs:
        pair.gene1_type = max(set(pair.gene1.tadb_type), key = pair.gene1.tadb_type.count)
        pair.gene1_role = max(set(pair.gene1.tadb_roles), key = pair.gene1.tadb_roles.count)
        pair.gene1_family = max(set(pair.gene1.tadb_families), key = pair.gene1.tadb_families.count)
        pair.gene1_name = max(set(pair.gene1.tadb_gene_name), key = pair.gene1.tadb_gene_name.count)
        pair.shared_tadb_ids = {}
        pair.shared_tadb_families = {}

        if pair.gene2 != None:
            pair.gene2_type = max(set(pair.gene2.tadb_type), key = pair.gene2.tadb_type.count)
            pair.gene2_role = max(set(pair.gene2.tadb_roles), key = pair.gene2.tadb_roles.count)
            pair.gene2_family = max(set(pair.gene2.tadb_families), key = pair.gene2.tadb_families.count)
            pair.gene2_name = max(set(pair.gene2.tadb_gene_name), key = pair.gene2.tadb_gene_name.count)
            pair.shared_tadb_ids = list(set(pair.gene1.tadb_ids).intersection(set(pair.gene2.tadb_ids)))
            pair.shared_tadb_families = list(set(pair.gene1.tadb_families).intersection(set(pair.gene2.tadb_families)))

            if pair.gene1_role == "T":
                pair.toxin = pair.gene1.id
                pair.toxin_name = pair.gene1_name
                pair.antitoxin = pair.gene2.id
                pair.antitoxin_name = pair.gene2_name
            else:
                pair.toxin = pair.gene2.id
                pair.toxin_name = pair.gene2_name
                pair.antitoxin = pair.gene1.id
                pair.antitoxin_name = pair.gene1_name

            if pair.gene1_role == pair.gene2_role:
                if len(set(pair.gene1.tadb_roles)) > 1 or len(set(pair.gene2.tadb_roles)) > 1:
                    pair.gene1_role = dict(Counter(pair.gene1.tadb_roles))
                    pair.gene2_role = dict(Counter(pair.gene2.tadb_roles))
                    pair.score = 8
                else:
                    pair.score = 1

                #! need to put logic for when multiple roles are given (family 6353)
            else:
                # start more strict
                if len(pair.shared_tadb_ids) > 0:
                    pair.score = 20
                elif pair.gene1_family == pair.gene2_family:
                    pair.score = 15
                elif len(pair.shared_tadb_families) > 0:
                    pair.score = 10
                else:
                    pair.score = 5

        # if pair.score > 5:
        #     print(f'--- PAIR #{pair.id} ---')
        #     print(f'Genes: {pair.gene1.id}\t{pair.gene2.id}')
        #     print(f'TA Type: {pair.gene1_type}\t{pair.gene2_type}')
        #     print(f'TA Role: {pair.gene1_role}\t{pair.gene2_role}')
        #     print(pair.gene1_family, pair.gene2_family, sep = "\t")
        #     print(pair.gene1_name, pair.gene2_name, sep = "\t")
        #     print(f'Shared TADB IDs: {pair.shared_tadb_ids}')
        #     print(f'SCORE: {pair.score} - {scores[pair.score]}')
        #     print("***")
        #     print(pair.__dict__)

        processed_pairs.append(pair)

    return(pairs)

scores = {1:  "Not candidates - have same role",
          5:  "Same TYPE and opposite ROLEs",
          8:  "Some ambiguity in ROLEs",
          10: "Have at least one TADB FAMILY prediction in common",
          15: "Best TADB FAMILY predictions are the same",
          20: "TA predictions have TADB IDs in common"}


def parse_tadb_results(gene):
    gene.tadb_families = []
    gene.tadb_gene_name = []
    gene.tadb_roles = []
    gene.tadb_ids = []
    gene.tadb_type = []

    if hasattr(gene, 'tadb_blast'):
        for blast_result in gene.tadb_blast:
            # blast_result.print_full_result()
            if blast_result.filter(e = 1e-15, p = 25):
                family, prediction, role, ta_id = get_tadb_family(blast_result.subject)
                gene.tadb_families += family
                gene.tadb_gene_name += prediction
                gene.tadb_roles.append(role)
                gene.tadb_ids.append(ta_id)
                gene.tadb_type.append(tadb_type.loc[int(ta_id), 'TA_TYPE'])

    return gene

def get_potential_tadb_pairs(genome, overlapping_ids = 1, max_distance = 500):
    print(f'{colors.bcolors.PURPLE}Finding Potential Pairs{colors.bcolors.END}')

    pairs = []

    genome.potential_TA_list = list(set(genome.potential_TA_list))
    for gene1 in natsorted(genome.potential_TA_list):
        for gene2 in natsorted(genome.potential_TA_list):
            distance = gene_distance(genome.genes[gene1], genome.genes[gene2])
            if distance <= max_distance:
                pairs.append(PAIR(genome.genes[gene1], genome.genes[gene2], len(pairs) + 1))

                # working to remove the code below

                # print(f'---\n{genome.genes[gene1].id} ({genome.genes[gene1].start}-{genome.genes[gene1].stop}) is {distance} bp from {genome.genes[gene2].id} ({genome.genes[gene2].start}-{genome.genes[gene2].stop})')
                # print(f'{genome.genes[gene1].id} is likely a type {max(set(genome.genes[gene1].tadb_type), key = genome.genes[gene1].tadb_type.count)} {max(set(genome.genes[gene1].tadb_roles), key = genome.genes[gene1].tadb_roles.count)}, and {genome.genes[gene2].id} is likely a type {max(set(genome.genes[gene2].tadb_type), key = genome.genes[gene2].tadb_type.count)} {max(set(genome.genes[gene2].tadb_roles), key = genome.genes[gene2].tadb_roles.count)}')
                # print(f'They have {len(set(genome.genes[gene1].tadb_ids).intersection(set(genome.genes[gene2].tadb_ids)))} TADB IDs in common.')
                # print(f'They are both predicted to be in the {set(genome.genes[gene1].tadb_families).intersection(set(genome.genes[gene2].tadb_families))} TADB family')
                # if len(genome.genes[gene1].tadb_gene_name) > 0:
                #     print(f'{genome.genes[gene1].id} is likely {max(set(genome.genes[gene1].tadb_gene_name), key = genome.genes[gene1].tadb_gene_name.count)}')
                # if len(genome.genes[gene2].tadb_gene_name) > 0:
                #     print(f'{genome.genes[gene2].id} is likely {max(set(genome.genes[gene2].tadb_gene_name), key = genome.genes[gene2].tadb_gene_name.count)}')
                #
                #
                # if max(set(genome.genes[gene1].tadb_roles), key = genome.genes[gene1].tadb_roles.count) == max(set(genome.genes[gene2].tadb_roles), key = genome.genes[gene2].tadb_roles.count):
                #     print(f'But they are the same role, so not counting as pairs!')
                # else:
                #     print(f'So let\'s make them pairs')
                #     genome.genes[gene1].tadb_pair.append(genome.genes[gene2].id)
                #     genome.genes[gene2].tadb_pair.append(genome.genes[gene1].id)

    # for gene in natsorted(genome.potential_TA_list):
    #     if len(genome.genes[gene].tadb_pair) == 0:
    #
    #         # add genes without pairs here
    #
    #         print(f'---\n{genome.genes[gene].id} ({genome.genes[gene].start}-{genome.genes[gene].stop}) isn\'t close to others, so is maybe an orphan gene')
    #         print(f'{genome.genes[gene].id} is likely a type {max(set(genome.genes[gene].tadb_type), key = genome.genes[gene].tadb_type.count)} {max(set(genome.genes[gene].tadb_roles), key = genome.genes[gene].tadb_roles.count)}')
    #         print(f'It has {len(set(genome.genes[gene].tadb_ids))} TADB IDs.')
    #         print(f'It is predicted to be in the {set(genome.genes[gene].tadb_families)} TADB family ')
    #         if len(genome.genes[gene].tadb_gene_name) > 0:
    #             print(f'{genome.genes[gene].id} is likely {max(set(genome.genes[gene].tadb_gene_name), key = genome.genes[gene].tadb_gene_name.count)}')

    return pairs

def find_TAs(genome):
    print(f'{colors.bcolors.PURPLE}Finding TAs{colors.bcolors.END}')

    # write genes to genomes and gene class dictionary
    genome.faa_path = os.path.join(args.out_dir, genome.name + ".faa")
    genome.nt_path  = os.path.join(args.out_dir, genome.name + ".ffn")
    genome.contig_path = os.path.join(args.out_dir, genome.name + ".fa")
    genome.genes = gbk_to_fasta.main(genome.file_path, write_faa = genome.faa_path, write_nt = genome.nt_path, write_contig = genome.contig_path, return_gene_dict = True)
    genome.potential_TA_list = []

    blast_tadb(genome, aa = True, nt = True)

    for gene in natsorted(genome.genes):
        genome.genes[gene] = parse_tadb_results(genome.genes[gene])

        if len(genome.genes[gene].tadb_roles) > 0:
            genome.potential_TA_list.append(genome.genes[gene].id)
            # print(genome.genes[gene].id, genome.genes[gene].replicon, genome.genes[gene].start, genome.genes[gene].stop, genome.genes[gene].strand, dict(Counter(genome.genes[gene].tadb_roles)), sorted(genome.genes[gene].tadb_ids), dict(Counter(genome.genes[gene].tadb_gene_name)), dict(Counter(genome.genes[gene].tadb_families)), sep = "\t")
        # else:
        #     print('***** NO TYPES *****', genome.genes[gene].id, genome.genes[gene].replicon, genome.genes[gene].start, genome.genes[gene].stop, genome.genes[gene].strand, dict(Counter(genome.genes[gene].tadb_roles)), sorted(genome.genes[gene].tadb_ids), dict(Counter(genome.genes[gene].tadb_gene_name)), dict(Counter(genome.genes[gene].tadb_families)), sep = "\t")

    pairs = get_potential_tadb_pairs(genome)

    processed_pairs = score_tadb_pairs(pairs)

    # print results
    print('GENOME', 'GENOME_PAIR', 'REPLICON', 'RANGE', 'TYPE', 'TOXIN_LOCUS', 'TOXIN_NAME', 'ANTITOXIN_LOCUS', 'ANTITOXIN_NAME', 'SHARED_FAMILIES', 'SHARED_TADB_IDS', 'SCORE', 'SCORE_NOTE', 'TOXIN_TADB_IDS', 'ANTITOXIN_TADB_IDS', 'TOXIN_NAMES', 'ANTITOXIN_NAMES', sep = "\t")
    for pair in processed_pairs:
        if pair.score > 5:
            pair.view(genome.name)

    # for pair in processed_pairs:
    #     print(pair.gene1.id, pair.gene2.id, pair.score, sep = "\t")


    # annotate_tasmania(genome.faa_path, TASmania_metadata_A)
    # annotate_tasmania(genome.faa_path, TASmania_metadata_T)

def blast_tadb(genome, aa = True, nt = False):
    print(f'{colors.bcolors.PURPLE}Starting BLASTs{colors.bcolors.END}')

    if aa:
        tadb_aa_results = blast.run_blast(type="prot",
                                  q=genome.faa_path,
                                  db="/Users/kimbrel1/Science/repos/jakomics/db/TADB_2.0/TADB2.faa",
                                  e=1e-7,
                                  make = True)

        for query in natsorted(tadb_aa_results.keys()):
            genome.genes[query].tadb_blast = tadb_aa_results[query]
            # for hit in tadb_aa_results[query]:
            #     hit.print_rough_result()

    if nt:
        tadb_nt_results = blast.run_blast(type="nucl",
                                  q=genome.nt_path,
                                  db="/Users/kimbrel1/Science/repos/jakomics/db/TADB_2.0/TADB2.ffn",
                                  e=1e-7,
                                  make = True)

        for query in natsorted(tadb_nt_results.keys()):
            genome.genes[query].tadb_blast = tadb_nt_results[query]
            # for hit in tadb_nt_results[query]:
            #     hit.print_rough_result()
    #
    #
    # tadb_contig_results = blast.run_blast(type="nucl",
    #                           q=genome.contig_path,
    #                           db="/Users/kimbrel1/Science/repos/jakomics/db/TADB_2.0/TADB2.ffn",
    #                           e=1e-7,
    #                           make = True)
    #
    # for query in natsorted(tadb_contig_results.keys()):
    #     genome.genes[query].tadb.append(tadb_contig_results[query])

def annotate_tasmania(file_path, df):
    print(file_path)

    for toxin, row in df.iterrows():
        hmm.run_hmmsearch(file_path, "temp.log", "temp_output.txt",
                          TASmania_HMMs_path + row['hmm_profile_id'] + '.hmm3')
        with open('temp_output.txt', 'rU') as input:
            for qresult in SearchIO.parse(input, 'hmmsearch3-domtab'):

                query_name = qresult.id  # sequence ID from fasta
                qlen = qresult.seq_len
                hits = qresult.hits
                # ?? full_score = qresult.score
                # ?? full_eval = qresult.evalue
                num_hits = len(hits)
                if num_hits > 0:
                    for i in range(0, num_hits):
                        target_name = hits[i].id
                        hit_evalue = hits[i].evalue  # evalue
                        hit_score = hits[i].bitscore
                        hmm_name = hits[i].accession
                        tlen = hits[i].seq_len

                        for hsp in hits[i]:
                            query_range = str(hsp.query_start + 1) + "-" + str(hsp.query_end)
                            hit_range = str(hsp.hit_start + 1) + "-" + str(hsp.hit_end)
                            print(row['hmm_cluster_group_id'], row['hmm_profile_id'], row['nearest_pfam_identifier'], row['nearest_pfam_desc'], query_name, qlen, num_hits, query_range,
                                  target_name, tlen, hit_evalue, hit_score, hmm_name, hit_range, sep="\t")

# MAIN LOOP ###################################################################


if __name__ == "__main__":

    genome_list = utilities.get_files(args.files, args.in_dir, ["gbk", "gbff", "gb"])

    for genome in genome_list:
        print(f'Processing {genome.name}')

        find_TAs(genome)
