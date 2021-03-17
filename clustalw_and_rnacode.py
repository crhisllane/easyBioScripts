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

    cmd1 = "/home/crhisllane.vasconcelos/clustalw-2.1/src/clustalw2 %s"%(clusterline)
    print ("\t", cmd1, "\n")
    os.system(cmd1)
    clusterline_Aln = re.sub(".fasta", ".aln" ,clusterline)
    clusterline_Pure = re.sub(".fasta", "" ,clusterline)
    cmd2 = "/home/crhisllane.vasconcelos/RNAcode-0.3/src/RNAcode %s -p 0.05 >> %s.tab"%(clusterline_Aln, clusterline_Pure)
    print ("\t", cmd2, "\n")
    os.system(cmd2)


with concurrent.futures.ProcessPoolExecutor(max_workers=50) as executor:
    executor.map(process_rnacode, allClusterLines)
