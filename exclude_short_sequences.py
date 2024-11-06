import argparse

bases = ["A", "G", "C", "T", "R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N"]

parser = argparse.ArgumentParser()
parser.add_argument("-len", "--length", type=int, dest="length", default=5)
parser.add_argument("-g", dest="genelist", default=None)
args = parser.parse_args()

genelist = open(args.genelist).readlines()

for genename in genelist:
    alignment = open(genename[:-1]).readlines()
    retained_sequences = []
    excluded_sequences = []
    for x in range(len(alignment)):
        if alignment[x][0] == ">":
            count = 0
            for base in alignment[x+1]:
                if base in bases:
                    count += 1
            if count >= args.length:
                retained_sequences.append(alignment[x])
                retained_sequences.append(alignment[x+1])
            else:
                excluded_sequences.append(alignment[x])
                excluded_sequences.append(alignment[x+1])
    if retained_sequences:
        with open(genename[:-1] + ".retained.fasta", "w") as f:
            for ff in retained_sequences:
                f.write(ff)
        f.close()
    if excluded_sequences:
        with open(genename[:-1] + ".excluded.fasta", "w") as f:
            for ff in excluded_sequences:
                f.write(ff)
        f.close()