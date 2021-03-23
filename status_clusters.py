import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures

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
cmd2 = "sed -i '/^[0-9]/ s/ /----/' FILECLUSTER.temp"
os.system(cmd2)

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



#for cdhitseq in cdhit_sequences:
def process_cdhitcluster(cdhitseq):
    ini = ">"    
    name_clu, seqs_clus = cdhitseq.id, str(cdhitseq.seq)
    name_clu = re.sub("^", ">" ,name_clu)
    print ("analise", name_clu, "\n")
    if name_clu in clusterok:
        print ("\tname_cluster ok", name_clu, "\n seqsclu", seqs_clus, "\n\n")
        ClusterName = []
        GeneName = []
        Tool = []
        OriginalFile = []
        Origem = []
        GeneLenght = []
        
        seqs_clus = re.sub("\.\.\.\*", "----" ,seqs_clus)
        seqs_clus = re.sub("%", "----" ,seqs_clus)

        each_id_seq = seqs_clus.split("----")

        for id_seq in each_id_seq:
            print ("idseq-", id_seq)
            if ini in id_seq:           
                print ("---- idser com", id_seq, "\n")
                linecomplt = seqs_clus.split(",")
                GeneLenght.append(re.sub("aa", "" ,linecomplt[0]))

                '''fastaname = re.sub("\.\.\..*", "" ,id_seq)
                fastaname1 = re.sub("^>[A-Z][A-Z]_", "" ,fastaname)
                fastaname1 = re.sub(".proteins.ffa", "" ,fastaname1)
                cab_id = re.search(r"_.*", fastaname1)
                realseqid = re.sub("^_", "" ,cab_id[0])
                fastaFileName = re.search(r".*\.fna", fastaname1)
                #print ("Antes-", fastaname, " ", realseqid, " ",  fastaFileName[0], "\n")
                if re.match(r"^>PD", fastaname):
                    fastaFileName = fastaFileName[0] + ".genes.fna"
                    #print ("Depois-", fastaname, " ", realseqid, " ",  fastaFileName, "\n")

                elif re.match(r"^>SM", fastaname):
                    fastaFileName = fastaFileName[0] + ".ffn"
                    #print ("Depois-", fastaname, " ", realseqid, " ",  fastaFileName, "\n")


                else:
                    log.write(fastaname + ' no ^>PD or ^>SM')'''
                
                #fasta_sequences = SeqIO.parse(open(fastaFileName),'fasta')
                #for fasta in fasta_sequences:
                #    name, sequence = fasta.id, str(fasta.seq)
                #    if re.search(r'.ffn', fastaFileName):
                #        print ("teste do SM0-", name, fastaname, " ", realseqid, " ",  fastaFileName, "\n") 
                #        name = re.sub("^[a-zA-Z0-9]*_[0-9]*_", "" ,name)
                #        print ("teste do SM1-", name, fastaname, " ", realseqid, " ",  fastaFileName, "\n")
                #    if realseqid == name:
                #        print ("teste do SM2-", name, fastaname, " ", realseqid, " ",  fastaFileName, "\n")
                #        new_sequence = fastaname + '\n' + sequence + '\n'
                #        outfasta.write(new_sequence)
                #        #print ("\t\t inserting ",  realseqid, " ", name, " ", fastaname, "\n") 
            else:

    else:
        print ("OUT cluster: less than size", name_clu, "\n")    


cdhit_sequences = SeqIO.parse(open("FILECLUSTER.temp"),'fasta')
with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
    executor.map(process_cdhitcluster, cdhit_sequences)
