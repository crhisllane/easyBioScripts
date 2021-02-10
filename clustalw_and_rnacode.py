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
    
    print (clusterline, "\n")

    cmd1 = "clustalw2 %s"%(clusterline)
    print ("\t", cmd1, "\n")
    os.system(cmd1)
    clusterline_Aln = re.sub(".fasta", ".aln" ,clusterline)
    clusterline_Pure = re.sub(".fasta", "" ,clusterline)
    cmd2 = "RNAcode %s -p 0.05 >> %s.tab"%(clusterline_Aln, clusterline_Pure)
    print ("\t", cmd2, "\n")
    os.system(cmd2)


with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_rnacode, allClusterLines)
