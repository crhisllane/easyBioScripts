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
    
    for rps_lines in rps_all_lines:
        #print("LINHA", rps_lines)

        #rps_one = rps_lines.rstrip('\n')
        #rps_elements, = rps_lines.split('\t')
        print("teste", rps_lines[0])
    
        #if re.match(r"Query=", rps_elements[0]):
        #    print ("bateu aqui", rps_elements[0], "\n")
        
        #if re.match(r"CDD:", rps_elements[0]):
        #    print ("-------- aqui", rps_elements[0], "\n")
            

path = os.getcwd()
path = str(path) + '/*.rpsblast' 
files_rpsblast =  glob.glob(path) 
    
with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_count, files_rpsblast)