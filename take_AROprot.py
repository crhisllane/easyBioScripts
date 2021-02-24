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
                   help = "eggnog output file")
argv.add_option("-f","--first", action="store", dest="firstid", type="string",
                    help = "first id in eggnog out" )

(argumentos, palha_que_nao_interessa) = argv.parse_args()

eggnog_file = open (argumentos.file, 'r')
eggnog_lines = eggnog_file.readlines()

Query = 0
Aro =["K00800", "K24018", "K13830"]

#KOs_complete = [Ptransport, OrganicPmine, InorganicPsolub, NitroToAmonia, AmoniaToHydroxy, HydroxyToNitrite, NitriteToNitrate, NitrateToNitrite, NitriteToNitricOxide, NitricOxideToNitrousOxide, NitrousOxideToNitrogen, NitriteToAmonia]

#OUT = open ('countKOeggnog.log', 'w+')
#OUT.write("Genome\tPtransport\tOrganicPmine\tInorganicPsolub\tNitroToAmonia\tAmoniaToHydroxy\tHydroxyToNitrite\tNitriteToNitrate\tNitrateToNitrite\tNitriteToNitricOxide\tNitricOxideToNitrousOxide\tNitrousOxideToNitrogen\tNitriteToAmonia\n")
#print("Genome\tPtransport\tOrganicPmine\tInorganicPsolub\tNitroToAmonia\tAmoniaToHydroxy\tHydroxyToNitrite\tNitriteToNitrate\tNitrateToNitrite\tNitriteToNitricOxide\tNitricOxideToNitrousOxide\tNitrousOxideToNitrogen\tNitriteToAmonia\n")
#print("Genome\t", Ptransport, "\t", OrganicPmine, "\t", InorganicPsolub, "\t", NitroToAmonia, "\t", AmoniaToHydroxy, "\t", HydroxyToNitrite, "\t", NitriteToNitrate, "\t", NitrateToNitrite, "\t", NitriteToNitricOxide, "\t", NitricOxideToNitrousOxide, "\t", NitrousOxideToNitrogen, "\t", NitriteToAmonia)
namefirst = argumentos.firstid


#ARO
Aro_all = len(Aro)
Aro_count = [0] * Aro_all



for eggnog_line in eggnog_lines:
    if not (eggnog_line.startswith("#")):
        rps_elements = eggnog_line.split("\t")
        name =  rps_elements[0]
        name = name.split("_")
        kos_annot = re.sub("ko:", "" ,rps_elements[14])
        kos_annot = kos_annot.split(",")

        combine1 = set(kos_annot) & set(Aro)
        tamanho1 = len(combine1)
             
        if tamanho1 >= 1:
            print(name, "\n")
            #comando = "grep %s *.faa | sed 's/_[0-9]* .*//' | sort | uniq"%(name)
            #resp = subprocess.check_output(comando, shell=True)
            #resp = str(resp)
            #print (resp)

#cdhit_sequences = SeqIO.parse(open("FILECLUSTER.temp"),'fasta')
#with concurrent.futures.ProcessPoolExecutor(max_workers=300) as executor:
#    executor.map(process_cdhitcluster, cdhit_sequences)