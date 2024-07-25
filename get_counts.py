#!/usr/bin/python3
import re
import os
import sys


def process_star_alignfile(f, strand='reverse'):
    GENE = 0
    STRAND_ = {'unstranded':1,
               'forward':2,
               'reverse':3}
    gd = {}
    with open(f,'r') as aligned:
        for line in aligned:
            if line.startswith('N_'):
                continue
            cur = line.strip().split('\t')
            gd[cur[GENE]] = cur[STRAND_[strand]]
            
    return gd

def main():

    # go through alignments tab and collect all files which will be read from
    alignfolder = sys.argv[1]
    if not os.path.exists(alignfolder):
        raise OSError("Align folder doesn't exist.")
        
    mappedreadfiles = [x for x in os.listdir(alignfolder) if re.search('ReadsPerGene',x)]
    filecounts = {}
    for x in mappedreadfiles:
        filecounts[x] = process_star_alignfile(f=os.path.join(alignfolder, x))
        
    allgenes = list(filecounts[mappedreadfiles[0]].keys())
    allgenes.sort()
    with open("counts.txt",'w') as counts:
        # write header
        header = ['Ensembl_ID'] + mappedreadfiles
        counts.write('\t'.join(header)+'\n')
        
        # write genes
        for gene in allgenes:
            genecounts = []
            genecounts.append(gene)
            for sample in mappedreadfiles:
                genecounts.append(str(filecounts[sample][gene]))
            counts.write('\t'.join(genecounts)+'\n')

if __name__=='__main__':
    main()