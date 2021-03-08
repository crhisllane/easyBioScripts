import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 
import subprocess

argv = OptionParser()

argv.add_option("-I", "--Input", action = "store", dest = "file", type ="string",
                   help = "fasta to change name")
argv.add_option("-f","--first", action="store", dest="out", type="string",
                    help = "name of output fasta" )

(argumentos, palha_que_nao_interessa) = argv.parse_args()

namefileout = argumentos.out
logout = argumentos.out
logout = logout + '.log' 
outfasta = open (namefileout, 'w+')
logfile = open (logout, 'w+')

count = 1

def take_sequence(sequenceFasta):
      
    name_seq, seqs = sequenceFasta.id, str(sequenceFasta.seq)
    
    new_sequence = '>' + count + '__' + '\n' + seqs + '\n'
    new_logline = count + '\t' + name_seq + '\n'
    outfasta.write(new_sequence)
    logfile.write(new_logline)

            

sequences = SeqIO.parse(open(argumentos.file),'fasta')
with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
    executor.map(take_sequence, sequences)