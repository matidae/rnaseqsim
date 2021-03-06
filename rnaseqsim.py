import gff_parser, splicer, sys
from Bio import SeqIO
from decimal import *
import distrib
from subprocess import call

gffvar = "brucei.gff"
readsvar = sys.argv[1]
gene_expression = sys.argv[2]
utr5 = 50
utr3 = 50
read_length = int(sys.argv[3])

#Parsing GFF file
gff_parser.parse(gffvar, gene_expression, utr5, utr3)

expresion_list = map(lambda x: int(x.rsplit()[1]), open(gene_expression, "r"))
#total_reads = sum(expresion_list)

def __printReads(gffvar, coef_var):
    n=0
#    exp_file = open(readsvar+".exp_normal","w")
    reads_file = open(readsvar+".fasta","w")
#    for i in SeqIO.parse(gffvar + ".utr.fasta", "fasta"):
    b = open("brucei.gff.utr.fasta","r")
    for i in SeqIO.parse(b, "fasta"):
        read_var = round(distrib.valueNormal(coef_var, int(expresion_list[n])))
        reads_list = splicer.generate_reads(i.seq, int(read_var), read_length)
        # exp_normal en RPKM
#        exp_normal = (Decimal(str(read_var)) / len(i.seq)) / total_reads
#        exp_file.write(i.name + "\t" + str(round(exp_normal*10000, 4)) + "\n")
        j = 0
        for read in reads_list:
            reads_file.write(">" + i.name + "_read_"+ str(j) + "\n")
            reads_file.write(str(read) + "\n")
            j+=1
        n+=1

__printReads(gffvar, 0.0025)
def test(coef_var, exp, out):
    exp_file = open(out,"w")
    expresion_list = map(lambda x: int(x.rsplit()[1]), open(exp, "r"))
    gene_list = map(lambda x: x.rsplit()[0], open(exp, "r"))
    n = 0
    for i in expresion_list:
        read_var = round(distrib.valueNormal(coef_var, int(expresion_list[n])))
        exp_file.write(str(read_var)+"\n")
        n+=1
