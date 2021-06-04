import pysam
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', dest='targets', type=str, default="/mnt/nas/users/julie/baitinginfo/targets.fa", help="""fasta file to which reads were aligned""")
parser.add_argument('-b', dest='bam', type=str, help="""bam alignment file""")
args = parser.parse_args()

def get_contig_names(fasta_file):
    """
    Gets contig names from a fasta file using SeqIO.
    :param fasta_file: Full path to uncompressed, fasta-formatted file
    :return: List of contig names.
    """
    contig_names = list()
    for contig in SeqIO.parse(fasta_file, 'fasta'):
        contig_names.append(contig.id)
    return contig_names

print('Gene\tHits\tGene Fraction')
contig_names = get_contig_names(args.targets)
bamfile = pysam.AlignmentFile(args.bam, 'rb')
# print(bamfile.reference_name())
for contig_name in contig_names:
    read_count = bamfile.count(contig=contig_name, read_callback='all')
    cov_counts = bamfile.count_coverage(contig=contig_name)
    cov_number = 0
    total_number = 0
    for a, c, g, t in zip(*cov_counts):   # doesn't matter whether the order is right
        total_number = total_number + 1
        if (((a + g) + c) + t) != 0:
            cov_number = cov_number + 1
    prop_covered = float(cov_number) / float(total_number)
    print('{}\t{}\t{}'.format(contig_name, read_count, prop_covered))
