#!/usr/bin/python2
import random

def __weighted_index(weights):
    rnd = random.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

#Wang et al. (2009) RNA-Seq: a revolutionary tool for transcriptomics. Nature reviews. Genetics, 10(1), pp.57-63.
#RNA fragmentation
def __coverage_bias (gene_len, read_length):
    weigth = []
    for i in xrange(0, gene_len-read_length):
        if i >= 0 and i < gene_len/5:
            weigth.append(0.05/(gene_len/5))
        elif i >= gene_len/5 and i < (gene_len - (gene_len/5)):
            weigth.append(0.9/(gene_len - (gene_len/5)))
        else:
            weigth.append(0.05/(gene_len/5))
    return weigth

#Aird, D. et al. (2011) Analyzing and minimizing PCR amplification bias in Illumina sequencing libraries. Genome Biology, 12(2) 
def __gc_bias(sequence, weigth, read_length):
    aux_read = ""
    aux_gc = 0
    for i in xrange(0,10):
        pos = __weighted_index(weigth)
        read = sequence[pos:pos+read_length]
        gc = (read.count("G") + read.count("C")) / float(len(read))
        if gc >= 0.12 and gc <= 0.65:
            return read
        else:
            aux_read = read
        return aux_read

#Dohm et al.(2008) Substantial biases in ultra-short read data sets from high-throughput DNA sequencing. Nucleic Acids Research, 36(16), p.e105.
def __weigth_errors(read_length):
    weigth_pos = list()
    n = 0
    while n <= read_length * .5:
        weigth_pos.append(0.01/(read_length*.5))
        n+=1
    while n > read_length * .5 and n <= read_length * .6:
        weigth_pos.append(0.015/(read_length*.6 - read_length*.5))
        n+=1
    while n > read_length * .6  and n <= read_length * .7:
        weigth_pos.append(0.02/(read_length*.7 - read_length*.6))
        n+=1
    while n > read_length * .7 and n <= read_length * .8:
        weigth_pos.append(0.03/(read_length*.8 - read_length*.7))
        n+=1
    while n > read_length * .8 and n <= read_length * .9:
        weigth_pos.append(0.04/(read_length*.9 - read_length*.8))
        n+=1
    while n > read_length * .9 and n <= read_length * .95:
        weigth_pos.append(0.065/(read_length*.95 - read_length*.9))
        n+=1
    while n > read_length * .95 and n < read_length:
        weigth_pos.append(0.085/(read_length - read_length*.95))
        n+=1
    weigth_pos.append(50)
    return weigth_pos

def __seq_error(read):
    weigth_pos = __weigth_errors(len(read))
    pos = __weighted_index(weigth_pos)
    if pos != len(read):
        return __wrong_base_call(read, pos)
    else:
        return read    

def __wrong_base_call(read, pos):
    read_aux = list(read)
    #print str(len(read)) + "-----"+ str(pos)
    if read_aux[pos] == "A":
        err = __weighted_index([0.05, 0.1, 0.25])#A>G, A>T, A>C
        if err == 0:
            read_aux[pos] =  "G"
        elif err == 1:
            read_aux[pos] = "T"
        else:
            read_aux[pos] = "C"

    elif read_aux[pos] == "T":
        err = __weighted_index([0.05, 0.06, 0.1])#T>A, T>G, T>C
        if err == 0:
            read_aux[pos] =  "A"
        elif err == 1:
            read_aux[pos] = "G"
        else:
            read_aux[pos] = "C"

    elif read_aux[pos] == "G":
        err = __weighted_index([0.02, 0.05, 0.17])#G>A, G>C, G>T
        if err == 0:
            read_aux[pos] =  "A"
        elif err == 1:
            read_aux[pos] = "C"
        else:
            read_aux[pos] = "T"

    elif read_aux[pos] == "C":
        err = __weighted_index([0.02, 0.05, 0.07])#C>G, C>T, C>A
        if err == 0:
            read_aux[pos] =  "G"
        elif err == 1:
            read_aux[pos] = "T"
        else:
            read_aux[pos] = "A"

    read_err = "".join(read_aux)
    return read_err

def generate_reads(sequence, numreads, read_length):
    gene_len = len(sequence)
    weigth = __coverage_bias(gene_len, read_length)
    read_list = []
    for i in xrange(0, numreads):
        read = __gc_bias(sequence, weigth, read_length)
        read_seq =  __seq_error(read)
        read_list.append(read_seq)
    return read_list

