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
def __coverage_bias (gene_len):
    weigth = []
    for i in xrange(0, gene_len-30):
        if i >= 0 and i < gene_len/10:
            weigth.append(0.1/(gene_len/10))
        elif i >= gene_len/10 and i < (gene_len - (gene_len/10)):
            weigth.append(0.5/(gene_len - (gene_len/10)))
        else:
            weigth.append(0.1/(gene_len/10))
    return weigth
 
#Dohm et al.(2008) Substantial biases in ultra-short read data sets from high-throughput DNA sequencing. Nucleic Acids Research, 36(16), p.e105.
def __gc_bias(sequence, weigth):
    aux_read = ""
    aux_gc = 0
    for i in xrange(0,3):
        pos = __weighted_index(weigth)
        read = sequence[pos:pos+30]
        gc = (read.count("G") + read.count("C")) / float(len(read))
        if aux_gc <= gc:
            aux_read = read
            aux_gc = gc
    return aux_read

def __seq_error(read):
    weigth_pos = [0.028, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.012, 0.015, 0.017, 0.02, 0.019, 0.021, 0.022, 0.025, 0.03, 0.035, 0.04, 0.042, 0.061, 0.065, 0.082, 0.087, 70]    
    pos = __weighted_index(weigth_pos)
    if pos != 30:
        return __wrong_base_call(read, pos)
    else:
        return read    

def __wrong_base_call(read, pos):
    read_aux = list(read)
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

def generate_reads(sequence, numreads):
    gene_len = len(sequence)
    weigth = __coverage_bias(gene_len)
    read_list = []
    for i in xrange(0, numreads):
        read = __gc_bias(sequence, weigth)
        read_seq =  __seq_error(read)
        read_list.append(read_seq)
    return read_list

