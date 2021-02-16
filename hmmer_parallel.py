import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures

argv = OptionParser()

argv.add_option("-I", "--Input", action = "store", dest = "file", type ="string",
                   help = "list of all fasta files do align")
(argumentos, palha_que_nao_interessa) = argv.parse_args()

clusterfile = open (argumentos.file, 'r')
allClusterLines = clusterfile.readlines()

#for clusterline in allClusterLines:
def process_rnacode(clusterline):
    clusterline = clusterline.rstrip('\n')
    print (clusterline, "\n")

    cmd1 = "hmmscan -o %s\_hmmer3.out --cpu 30 -E 0.01 --tblout %s\_hmmer3.tblout PlantSSPv1.hmm %s"%(clusterline, clusterline, clusterline)
    print ("\t", cmd1, "\n")
    os.system(cmd1)



with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_rnacode, allClusterLines)
