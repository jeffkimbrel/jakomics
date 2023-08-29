from Bio import SeqIO
import sys

from jakomics import colors
from jakomics.gene import GENE
from jakomics.utilities import string_to_hash


class GENOME():

    def __str__(self):
        return "<JAKomics GENOME class>"

    def __init__(self, file):
        self.file = file

    def genbank_to_fasta(self, write_contig=None, write_faa=None, write_nt=None, feature_identifier='locus_tag', return_gene_dict=False):
        '''
        self.gbk is a JAKomics FILE
        '''

        # RAST/Patric output by default has an incorrect header which leads to a BioPython warning. This will suppress all warnings.
        import warnings
        from Bio import BiopythonWarning
        warnings.simplefilter('ignore', BiopythonWarning)

        if write_contig != None:
            ct_file = open(write_contig, 'w')

        if write_faa != None:
            aa_file = open(write_faa, 'w')

        if write_nt != None:
            nt_file = open(write_nt, 'w')

        if return_gene_dict:
            gene_dict = {}

        for seq_record in SeqIO.parse(self.file.file_path, "genbank"):
            if write_contig != None:
                ct = str(seq_record.seq)
                ct_file.write(">" + seq_record.name + "\n" + ct + "\n")

            for feature in seq_record.features:
                if feature.type in ['CDS', 'ncRNA']:
                    """
                    Initially sets the gene id to the hash of the nucleotide sequence. 
                    This gets overwritten if a feature identifier is found
                    """

                    nt = str(feature.location.extract(seq_record).seq)
                    id = string_to_hash(nt) 

                    if feature_identifier in feature.qualifiers:
                        id = feature.qualifiers[feature_identifier][0]
                    elif 'db_xref' in feature.qualifiers:
                        for f in feature.qualifiers['db_xref']:
                            if f.startswith('RAST2:fig|'):
                                id = f.replace('RAST2:fig|', '')
                    else:
                        print(f"{colors.bcolors.YELLOW}WARNING: Gene without a {feature_identifier} was given a locus_tag of {id}{colors.bcolors.END}", file=sys.stderr)

                    if return_gene_dict:
                        gene = GENE(id)
                        gene.replicon = seq_record.name
                        gene.parse_gbk_location(feature.location)
                        if 'product' in feature.qualifiers:
                            gene.product = feature.qualifiers['product']
                        if 'EC_number' in feature.qualifiers:
                            gene.EC_number = feature.qualifiers['EC_number']

                    if write_faa != None:
                        if 'translation' in feature.qualifiers:
                            aa = feature.qualifiers['translation'][0]
                            aa_file.write(">" + id + "\n" + aa + "\n")

                            if return_gene_dict:
                                gene.aa = aa

                    if write_nt != None:
                        nt_file.write(">" + id + "\n" + nt + "\n")

                        if return_gene_dict:
                            gene.nt = nt

                    if return_gene_dict:
                        gene_dict[id] = gene

        if write_contig != None:
            ct_file.close()

        if write_faa != None:
            aa_file.close()

        if write_nt != None:
            nt_file.close()

        if return_gene_dict:
            return gene_dict


def increment_gbk_version(seqRecord):
    if "sequence_version" in seqRecord.annotations:
        version = seqRecord.annotations["sequence_version"]
        version += 1
        seqRecord.id = str(seqRecord.name) + "." + str(version)

    else:
        seqRecord.id = str(seqRecord.name) + "." + str(1)

    return(seqRecord)

def add_gbk_comment(seqRecord, newComment):

    # is there a comments field?
    if "comment" in seqRecord.annotations:
        comments = seqRecord.annotations["comment"]
    else:
        comments = []

    # check to see if it has already been formatted to a list
    if not isinstance(comments, list):
        comments = comments.split("\n")

    comments.append(newComment)
    seqRecord.annotations["comment"] = comments

    return(seqRecord)