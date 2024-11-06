import argparse

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-g", dest="gene_list", help="Input individual gene trees to be investigated.", default=None)
parser.add_argument("-s", dest="sample_list", help="A list of samples to be investigated", default=None)
args = parser.parse_args()

for gene in open(args.gene_list).readlines():
    alignment = open(gene[:-1]).readlines()
    newalignment = []
    for line in alignment:
        if line[0] == ">":
            #index = line.index(" ")
            newalignment.append(line[:])
        else:
            string = ""
            for base in line:
                if base in ["n", "N"]:
                    string += "-"
                else:
                    string += base
            newalignment.append(string)
    with open(gene[:-1] + ".new", "w") as f:
        for x in newalignment:
            f.write (x)
    f.close()
