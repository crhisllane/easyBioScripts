import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 

argv = OptionParser()

argv.add_option("-I", "--Input", action = "store", dest = "file", type ="string",
                   help = "cluster cdhit output file")
argv.add_option("-l","--list", action="store", dest="list", type="string",
                    help = "cdhit output fasta file" )
(argumentos, palha_que_nao_interessa) = argv.parse_args()

KOs_file = open (argumentos.list, 'r')
KOs_line = KOs_file.readlines()


def process_count(rpsfile):
    print("file", rpsfile, "\n")
  
    rps_open = open (rpsfile, 'r')
    rps_all_lines = rps_open.readlines()
    Query = 0
    for rps_lines in rps_all_lines:

        rps_elements = rps_lines.split(" ")
        #print("teste", rps_elements[0])
    
        if re.match(r"Query=", rps_elements[0]):
            Query = rps_elements[1]
            Query = Query.rstrip('\n')
            print ("\n", Query)
        if re.match(r"CDD:", rps_elements[0]):
            KO = rps_elements[2]
            KO = KO.rstrip('\n')
            print ("\t", KO)
            

path = os.getcwd()
path = str(path) + '/*.rpsblast' 
files_rpsblast =  glob.glob(path) 
    
with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_count, files_rpsblast)