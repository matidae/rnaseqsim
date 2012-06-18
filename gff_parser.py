import sys, os, re
from Bio import SeqIO

def __getNumChromosomes(gffvar):
   chromo_list=[]
   gff_file = open (gffvar, "r")
   for line in gff_file:
      if "##sequence-region" in line:
         chromo_list.append(line.split("\t")[1])
   return chromo_list

def __fastaChromosomes(gffvar):
   flag = False
   chromo_list = __getNumChromosomes(gffvar)
   gff_file = open(gffvar, "r")
   chromo_file = open (gffvar+".chr.fasta", "w")
   for line in gff_file:
      if flag:
         chromo_file.write(line)
      if ">"+chromo_list[0] in line:
         chromo_file.write(line)
         flag = True

def __printGene(gffvar, gff_file, gene_name, utr5, utr3):
    chrom = start = end = strand = ""
    fasta_file = open( gffvar + ".utr.fasta", "a")
    for line in gff_file:
       if line[0]!="#":
          if (re.findall("\\b" + gene_name + "\\b",line)) and line.split("\t")[2] == "gene":
             vector = line.split("\t")
             chrom = vector[0]
             start = vector[3]
             end = vector[4]
             strand = vector[6]
             break
    for i in SeqIO.parse(gffvar + ".chr.fasta", "fasta"):
       if i.id == chrom:
           if strand == "+":
              fasta_file.write(">" + gene_name + "\n")
              seq_line = i.seq[int(start) - 1 - utr5 : int(end) + utr3] + "\n"
              fasta_file.write(str(seq_line))
           else:
              fasta_file.write(">" + gene_name + "\n")
              seq = "\n" + i.seq[int(start) - 1 - utr3 : int(end) + utr5]
              fasta_file.write(str(seq.reverse_complement()))
    
def parse(gffvar, listvar, u5, u3):
    try:
        gff_file = open(gffvar, "r").readlines()
        list_file = open(listvar, "r").readlines()
        utr5 = int(u5)
        utr3 = int(u3)
    except IOError:
        print "File not found"
        return 0
    except ValueError:
        print "NaN !!"
        return 0

    if not os.path.isfile(gffvar + ".chr.fasta") or os.path.getsize(gffvar + ".chr.fasta") == 0:
         __fastaChromosomes(gffvar)
    genes_list = []
    if not os.path.isfile(gffvar + ".utr.fasta"):
        for i in list_file:
            gene_name = i.rsplit()[0].split("\t")[0]
            __printGene(gffvar, gff_file, gene_name, utr5, utr3)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Missing parameters (gff_file list_file utr5 utr3)"
    else:
        gffvar = sys.argv[1]
        listvar = sys.argv[2]
        u5 = sys.argv[3]
        u3 = sys.argv[4]
        parse(gffvar, listvar, u5, u3)
