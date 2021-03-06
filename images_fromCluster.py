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
argv.add_option("-m","--metafile", action="store", dest="metaf", type="string",
                    help = "metafile" )
(argumentos, palha_que_nao_interessa) = argv.parse_args()

clusterfile = open (argumentos.file, 'r')
allClusterLines = clusterfile.readlines()

metafilel = open (argumentos.metaf, 'r')
allmetaLines = metafilel.readlines()

count = 0
lastcount = 0
idclust = 0
limit = argumentos.size
clusterok = []

cmd = "sed 's/>Cluster />Cluster_/' %s > FILECLUSTER.temp"%(argumentos.file)
os.system(cmd)
cmd2 = "sed -i '/^[0-9]/ s/\t/----/' FILECLUSTER.temp"
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
    print (idclust, "=", lastcount)
    clusterok.append(idclust)

print("clusters ok")

proffile = "prodigal.clusters" 
smalfile = "smallOrf.clusters" 
prodigal = open (proffile, 'w+')
smallOrf = open (smalfile, 'w+')

NPA = "npa.clusters" 
PA = "pa.clusters" 
RA = "ra.clusters"
SOIL = "soil.clusters"
NPA = open (NPA, 'w+')
PA = open (PA, 'w+')
RA = open (RA, 'w+')
SOIL = open (SOIL, 'w+')

#for cdhitseq in cdhit_sequences:
def process_cdhitcluster(cdhitseq):
    ini = ">"    
    name_clu, seqs_clus = cdhitseq.id, str(cdhitseq.seq)
    name_clu = re.sub("^", ">" ,name_clu)
    print ("analise", name_clu, "\n")
    if name_clu in clusterok:
        print ("\tname_cluster ok", name_clu, "\n seqsclu", seqs_clus, "\n\n")
        name_clu2 = re.sub(">", "", name_clu)
     
        seqs_clus = re.sub("\.\.\.\*", "----" ,seqs_clus)
        seqs_clus = re.sub("%", "----" ,seqs_clus)

        each_id_seq = seqs_clus.split("----")

        for id_seq in each_id_seq:
            print ("idseq-", id_seq)
            if ini in id_seq:           
                linecomplt = id_seq.split(",")
                GeneLenght = re.sub("aa", "" ,linecomplt[0])
                completeN = linecomplt[1].split("_")
                Tool = re.sub('>', '', completeN[0])
                OrigF = re.sub("\..*", "", completeN[1])
                Origem = completeN[2]
                GeneName = linecomplt[1]
                print ("---- idser com", id_seq, "\n", "aqui-", Tool, OrigF, Origem, GeneName, '\n')

                with open('prodigal.clusters', 'a+') as p:
                    if re.search("PD",completeN[0]):
                        print (".........", name_clu, "\t", "toolPD", completeN[0])
                        p.write(name_clu + "\n")
                        p.close()


                with open('smallOrf.clusters', 'a+') as s:
                    if re.search("SM",completeN[0]):
                        print (".........", name_clu, "\t", "toolSM", completeN[0])
                        s.write(name_clu + "\n")
                        s.close()
                
                with open('npa.clusters', 'a+') as npa:
                    if (completeN[2]=="NPA"):
                        print (".........", name_clu, "\t", "npateste", completeN[2])
                        npa.write(name_clu + "\n")
                        npa.close()
                
                with open('pa.clusters', 'a+') as pa:
                    if (completeN[2]=="PA"):
                        print (".........", name_clu, "\t", "pateste", completeN[2])
                        pa.write(name_clu + "\n")
                        pa.close()

                with open('ra.clusters', 'a+') as ra:
                    if (completeN[2]=="RA"):
                        print (".........", name_clu, "\t", "rateste", completeN[2])
                        ra.write(name_clu + "\n")
                        ra.close()
                
                with open('soil.clusters', 'a+') as soil:
                    if (completeN[2]=="soil"):
                        print (".........", name_clu, "\t", "soilteste", completeN[2])
                        soil.write(name_clu + "\n")
                        soil.close()
                
                for metaLine in allmetaLines:
                    completeMeta = metaLine.split("\t")
                    if (completeMeta[0]==OrigF):
                        with open(completeMeta[4], 'a+') as Ori:
                                print (".........", name_clu, "\t", "Origemteste", OrigF, completeMeta[0], completeMeta[4])
                                Ori.write(name_clu + "\n")
                                Ori.close()




    else:
        print ("OUT cluster: less than size", name_clu, "\n")    


cdhit_sequences = SeqIO.parse(open("FILECLUSTER.temp"),'fasta')
with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
    executor.map(process_cdhitcluster, cdhit_sequences)


cmd3 = "sort prodigal.clusters | uniq > prodigal_sorted.clusters && mv prodigal_sorted.clusters prodigal.clusters"
os.system(cmd3)
cmd4 = "sort smallOrf.clusters | uniq > smallOrf_sorted.clusters && mv smallOrf_sorted.clusters smallOrf.clusters"
os.system(cmd4)
cmd6 = "sort npa.clusters | uniq > npa_sorted.clusters && mv npa_sorted.clusters npa.clusters"
os.system(cmd6)
cmd7 = "sort pa.clusters | uniq > pa_sorted.clusters && mv pa_sorted.clusters pa.clusters"
os.system(cmd7)
cmd8 = "sort ra.clusters | uniq > ra_sorted.clusters && mv ra_sorted.clusters ra.clusters"
os.system(cmd8)
cmd9 = "sort soil.clusters | uniq > soil_sorted.clusters && mv soil_sorted.clusters soil.clusters"
os.system(cmd9)