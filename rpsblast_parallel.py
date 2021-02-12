import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 

#before using this script it is necessary to break the fasta file into smaller numbers
#awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%60000==0){file=sprintf("group_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < completeSmall_cdhit.out 


argv = OptionParser()

argv.add_option("-d", "--database", action = "store", dest = "database", type ="string",
                   help = "data base to rpsblast")
(argumentos, palha_que_nao_interessa) = argv.parse_args()

db = argumentos.database
db = db.rstrip('\n')



#for clusterline in allClusterLines:
def process_rpsblast(clusterline):
    clusterline = clusterline.rstrip('\n')
    print (clusterline, "\n")
    clusterline_Pure = re.sub(".fasta", "" ,clusterline)
    cmd1 = "rpsblast+ -db %s -query %s -out %s_CDD.rpsblast -num_threads 40 -evalue 0.01  -outfmt '7 qseqid sseqid pident nident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxis frames qseq sseq'"%(db, clusterline, clusterline_Pure)
    os.system(cmd1)


path = os.getcwd()
path = str(path) + '/group*' 
files_fasta =  glob.glob(path) 
    
with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_rpsblast, files_fasta)
