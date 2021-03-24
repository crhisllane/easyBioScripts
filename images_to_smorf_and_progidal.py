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
lenghAA = 0
limiteAA = 50

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
        #print ("elemento0", element[0], "elemento1", element[1], "elemento2", element[2])
        if lastcount == 0:
            lenghAA = element[1]
            lenghAA = re.sub("aa,", "" ,lenghAA)
            lenghAA = int(lenghAA)

        lastcount=int(element[0]) + 1

if ((lastcount >= limit) and (lenghAA <= limiteAA)):
    print (idclust, "=", lastcount, "and", lenghAA)

    clusterok.append(idclust)

print("clusters ok")

proffile = "prodigal.clusters_menorque50" 
smalfile = "smallOrf.clusters_menorque50" 
prodigal = open (proffile, 'w+')
smallOrf = open (smalfile, 'w+')

NPA = "npa.clusters_menorque50" 
PA = "pa.clusters_menorque50" 
RA = "ra.clusters_menorque50"
SOIL = "soil.clusters_menorque50"
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
        print ("\tname_cluster ok", name_clu)
        name_clu2 = re.sub(">", "", name_clu)
     
        seqs_clus = re.sub("\.\.\.\*", "----" ,seqs_clus)
        seqs_clus = re.sub("%", "----" ,seqs_clus)

        each_id_seq = seqs_clus.split("----")

        for id_seq in each_id_seq:
            print ("\tidseq-", id_seq)
            if ini in id_seq:  
                linecomplt = id_seq.split(",")
                GeneLenght = re.sub("aa", "" ,linecomplt[0])
                print ("\ttamanho ok", GeneLenght)         
                completeN = linecomplt[1].split("_")
                Tool = re.sub('>', '', completeN[0])
                OrigF = re.sub("\..*", "", completeN[1])
                Origem = completeN[2]
                GeneName = linecomplt[1]
                print ("\t---- idser com", id_seq, "\n", "aqui-", Tool, OrigF, Origem, GeneName, '\n')

                with open('prodigal.clusters_menorque50', 'a+') as p:
                    if re.search("PD",completeN[0]):
                        print ("\t.........", name_clu, "\t", "toolPD", completeN[0])
                        p.write(name_clu + "\n")
                        p.close()


                with open('smallOrf.clusters_menorque50', 'a+') as s:
                    if re.search("SM",completeN[0]):
                        print ("\t.........", name_clu, "\t", "toolSM", completeN[0])
                        s.write(name_clu + "\n")
                        s.close()
                
                with open('npa.clusters_menorque50', 'a+') as npa:
                    if (completeN[2]=="NPA"):
                        print ("\t.........", name_clu, "\t", "npateste", completeN[2])
                        npa.write(name_clu + "\n")
                        npa.close()
                
                with open('pa.clusters_menorque50', 'a+') as pa:
                    if (completeN[2]=="PA"):
                        print ("\t.........", name_clu, "\t", "pateste", completeN[2])
                        pa.write(name_clu + "\n")
                        pa.close()

                with open('ra.clusters_menorque50', 'a+') as ra:
                    if (completeN[2]=="RA"):
                        print ("\t.........", name_clu, "\t", "rateste", completeN[2])
                        ra.write(name_clu + "\n")
                        ra.close()
                
                with open('soil.clusters_menorque50', 'a+') as soil:
                    if (completeN[2]=="soil"):
                        print ("\t.........", name_clu, "\t", "soilteste", completeN[2])
                        soil.write(name_clu + "\n")
                        soil.close()
                
                #for metaLine in allmetaLines:
                #    completeMeta = metaLine.split("\t")
                #    if (completeMeta[0]==OrigF):
                #        with open(completeMeta[4], 'a+') as Ori:
                #                print (".........", name_clu, "\t", "Origemteste", OrigF, completeMeta[0], completeMeta[4])
                #                Ori.write(name_clu + "\n")
                #                Ori.close()




    else:
        print ("OUT cluster: less than size", name_clu, "\n")    


cdhit_sequences = SeqIO.parse(open("FILECLUSTER.temp"),'fasta')
with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
    executor.map(process_cdhitcluster, cdhit_sequences)


cmd3 = "sort prodigal.clusters_menorque50 | uniq > prodigal_sorted.clusters_menorque50 && mv prodigal_sorted.clusters_menorque50 prodigal.clusters_menorque50"
os.system(cmd3)
cmd4 = "sort smallOrf.clusters_menorque50 | uniq > smallOrf_sorted.clusters_menorque50 && mv smallOrf_sorted.clusters_menorque50 smallOrf.clusters_menorque50"
os.system(cmd4)
cmd6 = "sort npa.clusters_menorque50 | uniq > npa_sorted.clusters_menorque50 && mv npa_sorted.clusters_menorque50 npa.clusters_menorque50"
os.system(cmd6)
cmd7 = "sort pa.clusters_menorque50 | uniq > pa_sorted.clusters_menorque50 && mv pa_sorted.clusters_menorque50 pa.clusters_menorque50"
os.system(cmd7)
cmd8 = "sort ra.clusters_menorque50 | uniq > ra_sorted.clusters_menorque50 && mv ra_sorted.clusters_menorque50 ra.clusters_menorque50"
os.system(cmd8)
cmd9 = "sort soil.clusters_menorque50 | uniq > soil_sorted.clusters_menorque50 && mv soil_sorted.clusters_menorque50 soil.clusters_menorque50"
os.system(cmd9)