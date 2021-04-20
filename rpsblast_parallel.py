import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 

#before using this script it is necessary to break the fasta file into smaller numbers
#using splitfasta?

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
    
    #o rpsblast trunca o id quando ocorre um - por isso ta sendo substituido
    cmd0 = "sed -i '/^>/s/-/____/' %s"%(clusterline)
    os.system(cmd0)

    clusterline_Pure = re.sub(".faa", "" ,clusterline)
    db_name = re.sub(".*/", "" ,db)
    cmd1 = "/strg1/home/crhisllane.vasconcelos/programas/ncbi-blast-2.11.0+/bin/rpsblast -db %s -query %s -out %s_%s.rpsblast -num_threads 40 -evalue 0.01"%(db, clusterline, clusterline_Pure, db_name)
    #cmd1 = "blastp -db %s  -query %s -out %s_%s.blast -num_threads 4 -evalue 0.01  -outfmt '7 qseqid sseqid pident nident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxis frames qseq sseq'"%(db, clusterline, clusterline_Pure, db_name)
    os.system(cmd1)

    #o rpsblast trunca o id quando ocorre um - foi substituido no inicio do script e agora estou colocando de volta
    cmd2 = "sed -i '/^>/s/____/-/' %s"%(clusterline)
    os.system(cmd2)
    cmd3 = "sed -i '/^Query=/s/____/-/' %s_%s.rpsblast"%(clusterline_Pure, db_name)


path = os.getcwd()
path = str(path) + '/*.faa' 
files_fasta =  glob.glob(path) 
    
with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_rpsblast, files_fasta)
