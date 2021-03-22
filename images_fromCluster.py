import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import numpy as np
import pandas as pd

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
prodigal = open ('prodigal.clusters', 'w+')
smallOrf = open ('smallOrf.clusters', 'w+')


#for cdhitseq in cdhit_sequences:
def process_cdhitcluster(cdhitseq):
    ini = ">"    
    name_clu, seqs_clus = cdhitseq.id, str(cdhitseq.seq)
    name_clu = re.sub("^", ">" ,name_clu)
    print ("analise", name_clu, "\n")
    if name_clu in clusterok:
        print ("\tname_cluster ok", name_clu, "\n")

        #fileout = name_clu + ".fasta"
        #fileout = re.sub(">", "" ,fileout)
        #outfasta = open (fileout, 'w+')
        seqs_clus = re.sub("aa,", "\t" ,seqs_clus)
        each_id_seq = seqs_clus.split("\t")
        ClusterName = []
        GeneName = []
        Tool = []
        OriginalFile = []
        Origem = []
        GeneLenght = []
        for id_seq in each_id_seq:
            #print (".........", name_clu, "\t", id_seq)
            if re.search('[a-zA-Z]', id_seq):
                completeN = id_seq.split("_")
                Tool.append(re.sub('>', '', completeN[0]))
                OrigF = re.sub("\..*", "", completeN[1])
                OriginalFile.append(OrigF)
                Origem.append(completeN[2])
                GeneName.append(id_seq)
                ClusterName.append(name_clu)

                #print (".........", name_clu, "\t", "namegene", nameGene)
            elif re.search('[1-9]', id_seq):
                GeneLenght.append(id_seq)
                #print (".........", name_clu, "\t", "lenghGene", lenghGene)
            
            '''only200seqs = int(only200seqs) + 1
            print ("--------squantidade - ", only200seqs)
            if (int(only200seqs) <= 200): 
                if ini in id_seq:           
                    #print ("\t---", id_seq, "\n")
                    fastaname = re.sub("\.\.\..*", "" ,id_seq)
                    fastaname1 = re.sub("^>[A-Z][A-Z]_", "" ,fastaname)
                    fastaname1 = re.sub(".proteins.ffa", "" ,fastaname1)
                    cab_id = re.search(r"_.*", fastaname1)
                    realseqid = re.sub("^_", "" ,cab_id[0])
                    fastaFileName = re.search(r".*\.fna", fastaname1)
                    #print (fastaname, " ", realseqid, " ",  fastaFileName[0], "\n")
                    if re.match(r"^>PD", fastaname):
                        fastaFileName = fastaFileName[0] + ".genes.fna"

                    elif re.match(r"^>SM", fastaname):
                        fastaFileName = fastaFileName[0] + ".ffn"

                    else:
                        log.write(fastaname + ' no ^>PD or ^>SM')
                    
                    fasta_sequences = SeqIO.parse(open(fastaFileName),'fasta')
                    for fasta in fasta_sequences:
                        name, sequence = fasta.id, str(fasta.seq)
                        if re.search(r'.ffn', fastaFileName): 
                            name = re.sub(">[a-zA-Z0-9]*_[a-zA-Z0-9]* ", ">" ,name)
                        if realseqid == name:
                            new_sequence = fastaname + '\n' + sequence + '\n'
                            outfasta.write(new_sequence)
                            #print ("\t\t inserting ",  realseqid, " ", name, " ", fastaname, "\n")'''
        
        #matriz = np.concatenate((ClusterName, GeneName, GeneLenght))

        for i in range(len(ClusterName)):
            print('----', ClusterName[i], OriginalFile[i], Tool[i], Origem[i], GeneLenght[i],  GeneName[i],'\n') 
            if re.search("PD", Tool[i]):
                #prodigal.write(ClusterName[i],'\n')
                print('PD', ClusterName[i], OriginalFile[i], Tool[i], Origem[i], GeneLenght[i],  GeneName[i],'\n')
                prodigal.write(ClusterName[i] + '\n') 
            else:
                #smallOrf.write(ClusterName[i],'\n')
                print('SM', ClusterName[i], OriginalFile[i], Tool[i], Origem[i], GeneLenght[i],  GeneName[i],'\n')
                smallOrf.write(ClusterName[i] + '\n') 


    else:
        print ("OUT cluster: less than size", name_clu, "\n")    


cdhit_sequences = SeqIO.parse(open("FILECLUSTER.temp"),'fasta')
with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
    executor.map(process_cdhitcluster, cdhit_sequences)


'''
#for cdhitseq in cdhit_sequences:
def process_cdhitcluster(cdhitseq):
    ini = ">"    
    name_clu, seqs_clus = cdhitseq.id, str(cdhitseq.seq)
    name_clu = re.sub("^", ">" ,name_clu)
    print ("analise", name_clu, "\n")
    if name_clu in clusterok:
        print ("\tname_cluster ok", name_clu, "\n")

        #fileout = name_clu + ".fasta"
        #fileout = re.sub(">", "" ,fileout)
        #outfasta = open (fileout, 'w+')
        seqs_clus = re.sub("aa,", "\t" ,seqs_clus)
        each_id_seq = seqs_clus.split("\t")
        ClusterName = []
        GeneName = []
        GeneLenght = []
        GeneLenght = []
        for id_seq in each_id_seq:
            #print (".........", name_clu, "\t", id_seq)
            if re.search('[a-zA-Z]', id_seq):
                completeN = id_seq.split("_")
               
                GeneName.append(id_seq)
                ClusterName.append(name_clu)

                #print (".........", name_clu, "\t", "namegene", nameGene)
            elif re.search('[1-9]', id_seq):
                GeneLenght.append(id_seq)
                #print (".........", name_clu, "\t", "lenghGene", lenghGene)
            
            only200seqs = int(only200seqs) + 1
            print ("--------squantidade - ", only200seqs)
            if (int(only200seqs) <= 200): 
                if ini in id_seq:           
                    #print ("\t---", id_seq, "\n")
                    fastaname = re.sub("\.\.\..*", "" ,id_seq)
                    fastaname1 = re.sub("^>[A-Z][A-Z]_", "" ,fastaname)
                    fastaname1 = re.sub(".proteins.ffa", "" ,fastaname1)
                    cab_id = re.search(r"_.*", fastaname1)
                    realseqid = re.sub("^_", "" ,cab_id[0])
                    fastaFileName = re.search(r".*\.fna", fastaname1)
                    #print (fastaname, " ", realseqid, " ",  fastaFileName[0], "\n")
                    if re.match(r"^>PD", fastaname):
                        fastaFileName = fastaFileName[0] + ".genes.fna"

                    elif re.match(r"^>SM", fastaname):
                        fastaFileName = fastaFileName[0] + ".ffn"

                    else:
                        log.write(fastaname + ' no ^>PD or ^>SM')
                    
                    fasta_sequences = SeqIO.parse(open(fastaFileName),'fasta')
                    for fasta in fasta_sequences:
                        name, sequence = fasta.id, str(fasta.seq)
                        if re.search(r'.ffn', fastaFileName): 
                            name = re.sub(">[a-zA-Z0-9]*_[a-zA-Z0-9]* ", ">" ,name)
                        if realseqid == name:
                            new_sequence = fastaname + '\n' + sequence + '\n'
                            outfasta.write(new_sequence)
                            #print ("\t\t inserting ",  realseqid, " ", name, " ", fastaname, "\n")
        
        #matriz = np.concatenate((ClusterName, GeneName, GeneLenght))

        for i in range(len(ClusterName)):
            print(ClusterName[i], GeneName[i], GeneLenght[i],'\n') 
    else:
        print ("OUT cluster: less than size", name_clu, "\n")    


cdhit_sequences = SeqIO.parse(open("FILECLUSTER.temp"),'fasta')
with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
    executor.map(process_cdhitcluster, cdhit_sequences)'''