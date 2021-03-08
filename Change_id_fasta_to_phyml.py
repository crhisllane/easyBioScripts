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

oldfasta = argumentos.file
namefileout = argumentos.out

logout = argumentos.out
logout = logout + '.log' 
outfasta = open (namefileout, 'w+')
logfile = open (logout, 'w+')

count = 1

for sequences in SeqIO.parse(open(oldfasta),'fasta'):
      
    name_seq, seqs = sequences.id, str(sequences.seq)
    
    new_sequence = '>' + str(count) + '__' + '\n' + seqs + '\n'
    new_logline = str(count) + '__' + '\t' + name_seq + '\n'
    outfasta.write(new_sequence)
    logfile.write(new_logline)
    count = count + 1
