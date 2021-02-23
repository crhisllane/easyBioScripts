import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 
import subprocess

argv = OptionParser()

argv.add_option("-K", "--KOcountfile", action = "store", dest = "file", type ="string",
                   help = "KO count tab")
argv.add_option("-f","--first", action="store", dest="firstid", type="string",
                    help = "first id in eggnog out" )

(argumentos, palha_que_nao_interessa) = argv.parse_args()

eggnog_file = open (argumentos.file, 'r')
eggnog_lines = eggnog_file.readlines()


for eggnog_line in eggnog_lines:
    if (eggnog_line.startswith("Genome")):
        print (eggnog_line)

    else:    
        rps_elements = eggnog_line.split("\t")
        name =  rps_elements[0]
        print (name, "\n")
         
        comando = "grep %s *.faa | sed 's/_[0-9]* .*//' | sort | uniq"%(name)
        resp = subprocess.check_output(comando, shell=True)
        resp = str(resp)
        resp = re.sub(".fna.out.faa..*", "" ,resp)
        resp = re.sub("^b\'", "" ,resp)
        outputLine = re.sub(name, resp ,eggnog_line)
        print (outputLine, "\n")