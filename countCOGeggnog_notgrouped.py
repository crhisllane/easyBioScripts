import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 
from numpy import empty
from array import array

argv = OptionParser()

argv.add_option("-I", "--Input", action = "store", dest = "file", type ="string",
                   help = "eggnog output file")
argv.add_option("-f","--first", action="store", dest="firstid", type="string",
                    help = "first id in eggnog out" )

(argumentos, palha_que_nao_interessa) = argv.parse_args()

eggnog_file = open (argumentos.file, 'r')
eggnog_lines = eggnog_file.readlines()

Query = 0
teste =["K03469", "K11814", "K03408", "K03415"]
Ptransport sugarTATPase= ["COG1129", "COG3839"]
OrganicPmine sugarTpermease = ["COG1175", "COG4158"]
InorganicPsolub sugarTperiplasmic = ["COG1879"]
NitroToAmonia glycerolTpermease = ["COG0395"]
AmoniaToHydroxy glycerolTperiplasmic = ["COG1653"]
HydroxyToNitrite TripartiteT= ["COG3181"]
NitriteToNitrate = ["K00370", "K00371"]
NitrateToNitrite = ["K00370", "K00371", "K00374", "K02567", "K02568", "K00367", "K10534", "K00372", "K00360"]
NitriteToNitricOxide = ["K00368", "K15864"]
NitricOxideToNitrousOxide = ["K04561", "K02305"]
NitrousOxideToNitrogen = ["K00376"]
NitriteToAmonia = ["K00366", "K00362", "K00363", "K03385", "K15876", "K17877"]

#KOs_complete = [Ptransport, OrganicPmine, InorganicPsolub, NitroToAmonia, AmoniaToHydroxy, HydroxyToNitrite, NitriteToNitrate, NitrateToNitrite, NitriteToNitricOxide, NitricOxideToNitrousOxide, NitrousOxideToNitrogen, NitriteToAmonia]

#OUT = open ('countKOeggnog.log', 'w+')
#OUT.write("Genome\tPtransport\tOrganicPmine\tInorganicPsolub\tNitroToAmonia\tAmoniaToHydroxy\tHydroxyToNitrite\tNitriteToNitrate\tNitrateToNitrite\tNitriteToNitricOxide\tNitricOxideToNitrousOxide\tNitrousOxideToNitrogen\tNitriteToAmonia\n")
print("Genome\tPtransport\tOrganicPmine\tInorganicPsolub\tNitroToAmonia\tAmoniaToHydroxy\tHydroxyToNitrite\tNitriteToNitrate\tNitrateToNitrite\tNitriteToNitricOxide\tNitricOxideToNitrousOxide\tNitrousOxideToNitrogen\tNitriteToAmonia\n")
print("Genome\t", Ptransport, "\t", OrganicPmine, "\t", InorganicPsolub, "\t", NitroToAmonia, "\t", AmoniaToHydroxy, "\t", HydroxyToNitrite, "\t", NitriteToNitrate, "\t", NitrateToNitrite, "\t", NitriteToNitricOxide, "\t", NitricOxideToNitrousOxide, "\t", NitrousOxideToNitrogen, "\t", NitriteToAmonia)
namefirst = argumentos.firstid

tab1 = 0
tab2 = 0
tab3 = 0
tab4 = 0
tab5 = 0
tab6 = 0
tab7 = 0
tab8 = 0
tab9 = 0
tab10 = 0
tab11 = 0
tab12 = 0

#Ptransport
Ptransport_all = len(Ptransport)
Ptransport_count = [0] * Ptransport_all

#OrganicPmine
OrganicPmine_all = len(OrganicPmine)
OrganicPmine_count = [0] * OrganicPmine_all

#InorganicPsolub
InorganicPsolub_all = len(InorganicPsolub)
InorganicPsolub_count = [0] * InorganicPsolub_all

#NitroToAmonia
NitroToAmonia_all = len(NitroToAmonia)
NitroToAmonia_count = [0] * NitroToAmonia_all

#AmoniaToHydroxy
AmoniaToHydroxy_all = len(AmoniaToHydroxy)
AmoniaToHydroxy_count = [0] * AmoniaToHydroxy_all

#HydroxyToNitrite
HydroxyToNitrite_all = len(HydroxyToNitrite)
HydroxyToNitrite_count = [0] * HydroxyToNitrite_all

#NitriteToNitrate
NitriteToNitrate_all = len(NitriteToNitrate)
NitriteToNitrate_count = [0] * NitriteToNitrate_all

#NitrateToNitrite
NitrateToNitrite_all = len(NitrateToNitrite)
NitrateToNitrite_count = [0] * NitrateToNitrite_all

#NitriteToNitricOxide
NitriteToNitricOxide_all = len(NitriteToNitricOxide)
NitriteToNitricOxide_count = [0] * NitriteToNitricOxide_all

#NitricOxideToNitrousOxide
NitricOxideToNitrousOxide_all = len(NitricOxideToNitrousOxide)
NitricOxideToNitrousOxide_count = [0] * NitricOxideToNitrousOxide_all

#NitrousOxideToNitrogen
NitrousOxideToNitrogen_all = len(NitrousOxideToNitrogen)
NitrousOxideToNitrogen_count = [0] * NitrousOxideToNitrogen_all

#NitriteToAmonia
NitriteToAmonia_all = len(NitriteToAmonia)
NitriteToAmonia_count = [0] * NitriteToAmonia_all


for eggnog_line in eggnog_lines:
    if not (eggnog_line.startswith("#")):
        rps_elements = eggnog_line.split("\t")
        name =  rps_elements[0]
        name = name.split("_")
        kos_annot = re.sub("ko:", "" ,rps_elements[14])
        kos_annot = kos_annot.split(",")

        if namefirst == name[0]:
            for ko in kos_annot:
             
                #Ptransport
                for i in range(Ptransport_all):
                    if ko == Ptransport[i]:
                        Ptransport_count[i] = Ptransport_count[i] + 1     
                
                #OrganicPmine
                for i in range(OrganicPmine_all):
                    if ko == OrganicPmine[i]:
                        OrganicPmine_count[i] = OrganicPmine_count[i] + 1 

                #InorganicPsolub
                for i in range(InorganicPsolub_all):
                    if ko == InorganicPsolub[i]:
                        InorganicPsolub_count[i] = InorganicPsolub_count[i] + 1  
                
                #NitroToAmonia
                for i in range(NitroToAmonia_all):
                    if ko == NitroToAmonia[i]:
                        NitroToAmonia_count[i] = NitroToAmonia_count[i] + 1 

                #AmoniaToHydroxy
                for i in range(AmoniaToHydroxy_all):
                    if ko == AmoniaToHydroxy[i]:
                        AmoniaToHydroxy_count[i] = AmoniaToHydroxy_count[i] + 1 
            
                #HydroxyToNitrite
                for i in range(HydroxyToNitrite_all):
                    if ko == HydroxyToNitrite[i]:
                        HydroxyToNitrite_count[i] = HydroxyToNitrite_count[i] + 1

                #NitriteToNitrate
                for i in range(NitriteToNitrate_all):
                    if ko == NitriteToNitrate[i]:
                        NitriteToNitrate_count[i] = NitriteToNitrate_count[i] + 1 

                #NitrateToNitrite
                for i in range(NitrateToNitrite_all):
                    if ko == NitrateToNitrite[i]:
                        NitrateToNitrite_count[i] = NitrateToNitrite_count[i] + 1 

                #NitriteToNitricOxide
                for i in range(NitriteToNitricOxide_all):
                    if ko == NitriteToNitricOxide[i]:
                        NitriteToNitricOxide_count[i] = NitriteToNitricOxide_count[i] + 1

                #NitricOxideToNitrousOxide
                for i in range(NitricOxideToNitrousOxide_all):
                    if ko == NitricOxideToNitrousOxide[i]:
                        NitricOxideToNitrousOxide_count[i] = NitricOxideToNitrousOxide_count[i] + 1 

                #NitrousOxideToNitrogen
                for i in range(NitrousOxideToNitrogen_all):
                    if ko == NitrousOxideToNitrogen[i]:
                        NitrousOxideToNitrogen_count[i] = NitrousOxideToNitrogen_count[i] + 1  

                #NitriteToAmonia
                for i in range(NitriteToAmonia_all):
                    if ko == NitriteToAmonia[i]:
                        NitriteToAmonia_count[i] = NitriteToAmonia_count[i] + 1

        else:
            print (namefirst, "\t", Ptransport_count, "\t", OrganicPmine_count, "\t", InorganicPsolub_count, "\t", NitroToAmonia_count, "\t", AmoniaToHydroxy_count, "\t", HydroxyToNitrite_count, "\t", NitriteToNitrate_count, "\t", NitrateToNitrite_count, "\t", NitriteToNitricOxide_count, "\t", NitricOxideToNitrousOxide_count, "\t", NitrousOxideToNitrogen_count, "\t", NitriteToAmonia_count)
            namefirst = name[0]

            #Ptransport
            Ptransport_count = [0] * Ptransport_all

            #OrganicPmine
            OrganicPmine_count = [0] * OrganicPmine_all

            #InorganicPsolub
            InorganicPsolub_count = [0] * InorganicPsolub_all

            #NitroToAmonia
            NitroToAmonia_count = [0] * NitroToAmonia_all 

            #AmoniaToHydroxy
            AmoniaToHydroxy_count = [0] * AmoniaToHydroxy_all 

            #HydroxyToNitrite
            HydroxyToNitrite_count = [0] * HydroxyToNitrite_all

            #NitriteToNitrate
            NitriteToNitrate_count = [0] * NitriteToNitrate_all

            #NitrateToNitrite
            NitrateToNitrite_count = [0] * NitrateToNitrite_all

            #NitriteToNitricOxide
            NitriteToNitricOxide_count = [0] * NitriteToNitricOxide_all

            #NitricOxideToNitrousOxide
            NitricOxideToNitrousOxide_count = [0] * NitricOxideToNitrousOxide_all

            #NitrousOxideToNitrogen
            NitrousOxideToNitrogen_count = [0] * NitrousOxideToNitrogen_all

            #NitriteToAmonia
            NitriteToAmonia_count = [0] * NitriteToAmonia_all

            for ko in kos_annot:
                                
                #Ptransport
                for i in range(Ptransport_all):
                    if ko == Ptransport[i]:
                        Ptransport_count[i] =  1  

                #OrganicPmine
                for i in range(OrganicPmine_all):
                    if ko == OrganicPmine[i]:
                        OrganicPmine_count[i] = 1 

                #InorganicPsolub
                for i in range(InorganicPsolub_all):
                    if ko == InorganicPsolub[i]:
                        InorganicPsolub_count[i] = 1 

                #NitroToAmonia
                for i in range(NitroToAmonia_all):
                    if ko == NitroToAmonia[i]:
                        NitroToAmonia_count[i] = 1  

                #AmoniaToHydroxy
                for i in range(AmoniaToHydroxy_all):
                    if ko == AmoniaToHydroxy[i]:
                        AmoniaToHydroxy_count[i] = 1 

                #HydroxyToNitrite
                for i in range(HydroxyToNitrite_all):
                    if ko == HydroxyToNitrite[i]:
                        HydroxyToNitrite_count[i] = 1 

                #NitriteToNitrate
                for i in range(NitriteToNitrate_all):
                    if ko == NitriteToNitrate[i]:
                        NitriteToNitrate_count[i] = 1

                #NitrateToNitrite
                for i in range(NitrateToNitrite_all):
                    if ko == NitrateToNitrite[i]:
                        NitrateToNitrite_count[i] = 1

                #NitriteToNitricOxide
                for i in range(NitriteToNitricOxide_all):
                    if ko == NitriteToNitricOxide[i]:
                        NitriteToNitricOxide_count[i] = 1 

                #NitricOxideToNitrousOxide
                for i in range(NitricOxideToNitrousOxide_all):
                    if ko == NitricOxideToNitrousOxide[i]:
                        NitricOxideToNitrousOxide_count[i] = 1 

                #NitrousOxideToNitrogen
                for i in range(NitrousOxideToNitrogen_all):
                    if ko == NitrousOxideToNitrogen[i]:
                        NitrousOxideToNitrogen_count[i] = 1 

                #NitriteToAmonia
                for i in range(NitriteToAmonia_all):
                    if ko == NitriteToAmonia[i]:
                        NitriteToAmonia_count[i] = 1  

print (namefirst, "\t", Ptransport_count, "\t", OrganicPmine_count, "\t", InorganicPsolub_count, "\t", NitroToAmonia_count, "\t", AmoniaToHydroxy_count, "\t", HydroxyToNitrite_count, "\t", NitriteToNitrate_count, "\t", NitrateToNitrite_count, "\t", NitriteToNitricOxide_count, "\t", NitricOxideToNitrousOxide_count, "\t", NitrousOxideToNitrogen_count, "\t", NitriteToAmonia_count)