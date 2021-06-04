# Julie Shay
# January 18, 2019
# This script will take a cluster file and a folder with target alignment summaries as input,
# and will output a table with counts for the most abundant gene per cluster per sample
import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='infolder', nargs='+', required=True, help="""Input Folder(s)""")
parser.add_argument('-s', dest='suffix', nargs='+', help="""Suffix for samples coming from each folder""")
parser.add_argument('-n', dest='fname', type=str, default="target_new_counts.txt", help="""Name of target alignment summary file within each subfolder""")
parser.add_argument('-flagstat', type=str, help="""Optional flagstat file name within each subfolder. If you set this, the script will output a column for proportion of reads within each sample.""")
parser.add_argument('-c', dest='clusters', nargs='+', default=["/mnt/nas/users/julie/baitinginfo/210426_calvinrequest/allthethings_amr.tsv", "/mnt/nas/users/julie/baitinginfo/210426_calvinrequest/allthethings_plasmid.tsv"], help="""Tab-delimited files with genes and clusters""")
parser.add_argument('-p', dest='minprop', type=float, default=0, help="""minimum proportion of gene covered""")
parser.add_argument('-o', dest='outprefix', type=str, default="temp_", help="""output file prefix""")
parser.add_argument('-drop', type=bool, default=False, help="""drop Gene Fraction column from output""")
args = parser.parse_args()

annodfs = []
grouped = []
# all of the annotation files must be tab-separated.
# They must also all have columns called "ID" and "Cluster".
# There can also be other columns with pretty much anything else in them.
for a in args.clusters:
    annodfs.append(pd.read_csv(a, sep="\t"))
    grouped.append(pd.DataFrame([], columns=["Cluster"]))

suf = ""
for d in range(0, len(args.infolder)):
    if args.suffix is not None:
        if args.suffix[d]:
            suf = "_" + args.suffix[d]
        else:
            suf = ""
    else:
        suf = ""
    for sample in os.listdir(args.infolder[d]):
        sd = args.infolder[d] + "/" + sample + "/"
        if os.path.isfile(sd + args.fname):
            tempdf = pd.read_csv(sd + args.fname, header=0, names=["ID", "Hits_" + sample + suf, "Gene_Fraction_" + sample + suf], sep="\t")
            tempdf = tempdf[tempdf["Hits_" + sample + suf] != 0]
            if args.flagstat is not None:
                if os.path.isfile(sd + args.flagstat):
                    # get total number of reads
                    with open(sd + args.flagstat) as f:
                        nreads = float(f.readline().split()[0])
                    # divide hit counts by total number of reads
                    tempdf["Prop_reads_" + sample + suf] = tempdf["Hits_" + sample + suf].apply(lambda x: float(x) / nreads)
            sampledf = tempdf[tempdf["Gene_Fraction_" + sample + suf] >= args.minprop]
            for i in range(0, len(annodfs)):
                tempdf = pd.merge(annodfs[i], sampledf, on = "ID", how="inner")
                tempdf.drop(annodfs[i].columns.drop("Cluster").tolist(), axis=1, inplace=True)
                tempdf2 = tempdf.groupby("Cluster").aggregate(np.max)
                if args.drop:
                    tempdf2 = tempdf2.drop(["Gene_Fraction_" + sample + suf], axis=1).rename(columns={"Hits_" + sample + suf: sample + suf})
                grouped[i] = pd.merge(grouped[i], tempdf2, on="Cluster", how="outer")
            del tempdf

for i in range(0, len(annodfs)):
	g = pd.merge(annodfs[i].drop(["ID"], axis=1).drop_duplicates(), grouped[i], on="Cluster", how="right")
	g.set_index("Cluster").to_csv((args.outprefix + os.path.splitext(os.path.basename(args.clusters[i]))[0] + ".tsv"), sep="\t")

