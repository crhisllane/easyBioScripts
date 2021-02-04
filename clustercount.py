import sys
import os
import re
from optparse import OptionParser

argv = OptionParser()

argv.add_option("-I", "--Input", action = "store", dest = "file", type ="string",
                   help = "cluster cdhit output file")
argv.add_option("-s", "--size", action="store", dest="size", type="int",
                    help = "minimum cluster size")
argv.add_option("-f","--fasta", action="store", dest="fasta", type="string",
                    help = "cdhit output fasta file" )
(argumentos, palha_que_nao_interessa) = argv.parse_args()

clusterfile = open (argumentos.file, 'r')
allClusterLines = clusterfile.readlines()
count = 0
lastcount = 0
idclust = 0
limit = argumentos.size
clusterok = []

for clusterline in allClusterLines:
    elements = clusterline.rstrip('\n')
    element = elements.split()

    if re.match(r"^>", element[0]):
        if (lastcount >= limit):
            #print (idclust, "=", lastcount)
            clusterok.append(idclust)

        count = 0
        idclust = element[0] + " " + str(element[1])
        #print ("idclus", element[0], element[1])
    else:
        #print ("  idseq", element[0], element[2])
        lastcount=int(element[0]) + 1

if (lastcount >= limit):
    #print (idclust, "=", lastcount)
    clusterok.append(idclust)

print(clusterok)
log = open ('ERROR.log', 'w+')

teste = 0
for idcluok in clusterok: 
    for clusterline in allClusterLines:
        elements = clusterline.rstrip('\n')
        element = elements.split()
        fileout2 = 0
        if re.match(r"^>", element[0]):
            teste = 0
            idclust = element[0] + " " + str(element[1])
            
            if idcluok == idclust:
                teste = 1
                fileout = element[0] + "_" + str(element[1]) + ".fasta"
                fileout2 = re.sub(">", "" ,fileout)

                outfasta = open (fileout2, 'w+')  
        else:
            if teste == 1:
                fastaname = re.sub("\.\.\.", "" ,element[2])
                fastaname1 = re.sub("^>[A-Z][A-Z]_", "" ,fastaname)
                fastaname1 = re.sub(".proteins.ffa", "" ,fastaname1)
                cabid = re.search(r"_.*", fastaname1)
                cabid = re.sub("^_", "" ,cabid[0])
                fastaFileName = re.search(r".*\.fna", fastaname1)

                if re.match(r"^>PD", fastaname):
                    fastaFileName = fastaFileName[0] + ".genes.fna" 
                elif re.match(r"^>SM", fastaname):
                    fastaFileName = fastaFileName[0] + ".ffa" 
                else:
                    log.write(fastaname + ' no ^>PD or ^>SM')


                print (fastaFileName, "\t", cabid)

                fasta_sequences = SeqIO.parse(open(fastaFileName),'fasta')
                for fasta in fasta_sequences:
                    name, sequence = fasta.id, str(fasta.seq)
                    new_sequence = fastaname + '\n' + sequence + '\n'
                    outfasta.write(new_sequence)

