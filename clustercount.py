import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO

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

cmd = "sed 's/>Cluster />Cluster_/' %s > FILECLUSTER.temp"%(argumentos.file)
os.system(cmd)

#taking the clusters with number of members according to the size informed.
for clusterline in allClusterLines:
    elements = clusterline.rstrip('\n')
    element = elements.split()

    if re.match(r"^>", element[0]):
        if (lastcount >= limit):
            #print (idclust, "=", lastcount)
            clusterok.append(idclust)

        count = 0
        idclust = element[0] + "_" + str(element[1])
        #print ("idclus", element[0], element[1])
    else:
        #print ("  idseq", element[0], element[2])
        lastcount=int(element[0]) + 1

if (lastcount >= limit):
    #print (idclust, "=", lastcount)
    clusterok.append(idclust)

print("clusters ok")
log = open ('ERROR.log', 'w+')

teste = 0

cdhit_sequences = SeqIO.parse(open("FILECLUSTER.temp"),'fasta')
for cdhitseq in cdhit_sequences:
    
    name_clu, seqs_clus = cdhitseq.id, str(cdhitseq.seq)
    name_clu = re.sub("^", ">" ,name_clu)
    print ("analise", name_clu, "\n")
    if name_clu in clusterok:
        print ("name_clu", name_clu, "\n")
        seqs_clus = re.sub("[0-9]*aa,", "----" ,seqs_clus)
        each_id_seq = seqs_clus.split("----")
        for id_seq in each_id_seq:            
            print ("\t---", id_seq, "\n")
    else:
        print ("OUT cluster: less than size", name_clu, "\n")    


"""
for clusterline in allClusterLines:   
    elements = clusterline.rstrip('\n')
    element = elements.split()
    fileout2 = 0
    if re.match(r"^>", element[0]):
        print ("take sequences from",  clusterline, "\n") 
        idclust = element[0] + " " + str(element[1])
        fileout = element[0] + "_" + str(element[1]) + ".fasta"
        fileout2 = re.sub(">", "" ,fileout)
        equalcluster = "nao"

        if idclust in clusterok:
            equalcluster = "tem"
            outfasta = open (fileout2, 'w+')
            print ("\t making a file",  idclust, " ", fileout2, "\n")
        else:
            print ("\t OUT cluster: less than size",  idclust, "\n") 
  
    else:
        if equalcluster == "tem":
            fastaname = re.sub("\.\.\.", "" ,element[2])
            fastaname1 = re.sub("^>[A-Z][A-Z]_", "" ,fastaname)
            fastaname1 = re.sub(".proteins.ffa", "" ,fastaname1)
            cabid = re.search(r"_.*", fastaname1)
            cabid = re.sub("^_", "" ,cabid[0])
            fastaFileName = re.search(r".*\.fna", fastaname1)
            if re.match(r"^>PD", fastaname):
                fastaFileName = fastaFileName[0] + ".genes.fna"

            elif re.match(r"^>SM", fastaname):
                fastaFileName = fastaFileName[0] + ".ffn"

            else:
                log.write(fastaname + ' no ^>PD or ^>SM')


            print ("\t", fastaFileName, "\t", cabid)

            fasta_sequences = SeqIO.parse(open(fastaFileName),'fasta')
            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                if re.search(r'.ffn', fastaFileName): 
                    name = re.sub(">[a-zA-Z0-9]*_[a-zA-Z0-9]* ", ">" ,name)
                if cabid == name:
                    new_sequence = fastaname + '\n' + sequence + '\n'
                    outfasta.write(new_sequence)
                    print ("\t\t inserting ",  cabid, " ", name, " ", fastaname, "\n") 

"""
