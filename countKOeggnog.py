import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 

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
Ptransport = ["K08176", "K16322", "K03306", "K03324", "K14683", "K02036", "K02037", "K02038", "K02040", "K02041", "K02044", "K02042", "K05814", "K05813", "K05816", "K05815"]
OrganicPmine = ["K02043", "K06166", "K06165", "K06164", "K06163", "K05781", "K05780", "K06162", "K05774", "K09994", "K06167", "K06193", "K03430", "K05306", "K01083", "K01093", "K01126", "K07048", "K01077", "K01113", "K09474", "K03788", "K01078"]
InorganicPsolub = ["K01507", "K01524", "K00117", "PF03600", "PF16980"]
NitroToAmonia = ["K02588", "K02586", "K02591"]
AmoniaToHydroxy = ["K10944", "K10945", "K10946"]
HydroxyToNitrite = ["K10535", "PF13435"]
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

for eggnog_line in eggnog_lines:
    if not (eggnog_line.startswith("#")):
        rps_elements = eggnog_line.split("\t")
        name =  rps_elements[0]
        name = name.split("_")
        kos_annot = re.sub("ko:", "" ,rps_elements[14])
        kos_annot = kos_annot.split(",")

        combine1 = set(kos_annot) & set(Ptransport)
        tamanho1 = len(combine1)

        combine2 = set(kos_annot) & set(OrganicPmine)
        tamanho2 = len(combine2)

        combine3 = set(kos_annot) & set(InorganicPsolub)
        tamanho3 = len(combine3)

        combine4 = set(kos_annot) & set(NitroToAmonia)
        tamanho4 = len(combine4)

        combine5 = set(kos_annot) & set(AmoniaToHydroxy)
        tamanho5 = len(combine5)

        combine6 = set(kos_annot) & set(HydroxyToNitrite)
        tamanho6 = len(combine6)

        combine7 = set(kos_annot) & set(NitriteToNitrate)
        tamanho7 = len(combine7)

        combine8 = set(kos_annot) & set(NitrateToNitrite)
        tamanho8 = len(combine8)

        combine9 = set(kos_annot) & set(NitriteToNitricOxide)
        tamanho9 = len(combine9)

        combine10 = set(kos_annot) & set(NitricOxideToNitrousOxide)
        tamanho10 = len(combine10)

        combine11 = set(kos_annot) & set(NitrousOxideToNitrogen)
        tamanho11 = len(combine11)

        combine12 = set(kos_annot) & set(NitriteToAmonia)
        tamanho12 = len(combine12)

        #print("-----tem", combine, "\t", tamanho, "\n")
        if namefirst == name[0]:
            #print("-----tem", combine, "\t", tamanho, "\n")
            tab1 = tamanho1 + tab1
            tab2 = tamanho2 + tab2
            tab3 = tamanho3 + tab3
            tab4 = tamanho4 + tab4
            tab5 = tamanho5 + tab5
            tab6 = tamanho6 + tab6
            tab7 = tamanho7 + tab7
            tab8 = tamanho8 + tab8
            tab9 = tamanho9 + tab9
            tab10 = tamanho10 + tab10
            tab11 = tamanho11 + tab11
            tab12 = tamanho12 + tab12

        else:
            print(namefirst,"\t",tab1,"\t",tab2,"\t",tab3,"\t",tab4,"\t",tab5,"\t",tab6,"\t",tab7,"\t",tab8,"\t",tab9,"\t",tab10,"\t",tab11,"\t",tab12,"\n")
            #saida = (" %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n") %(namefirst, tab1, tab2,  tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10, tab11, tab12)
            #OUT.write = saida
            namefirst = name[0]
            tab1 = tamanho1
            tab2 = tamanho2
            tab3 = tamanho3
            tab4 = tamanho4
            tab5 = tamanho5
            tab6 = tamanho6
            tab7 = tamanho7
            tab8 = tamanho8
            tab9 = tamanho9
            tab10 = tamanho10
            tab11 = tamanho11
            tab12 = tamanho12
        #if re.match(r"Query=", rps_elements[0]):
            #print(Query, "\t", KO, "\n")
        #    KO = []
        #    Query = rps_elements[1]
        #    Query = Query.rstrip('\n')
            #print ("\n", Query)
        #if re.match(r"CDD:", rps_elements[0]):
        #    KOone = rps_elements[2]
        #    KOone = KOone.rstrip('\n')
        #    KO.append(KOone)
        #   if KOone in KOs_line:
        #        print (Query, "\t", KOone, "\n")
        
        #print(Query, "\t", KO, "\n")

print(namefirst,"\t",tab1,"\t",tab2,"\t",tab3,"\t",tab4,"\t",tab5,"\t",tab6,"\t",tab7,"\t",tab8,"\t",tab9,"\t",tab10,"\t",tab11,"\t",tab12,"\n")
#OUT.write(namefirst,"\t",tab1,"\t",tab2,"\t",tab3,"\t",tab4,"\t",tab5,"\t",tab6,"\t",tab7,"\t",tab8,"\t",tab9,"\t",tab10,"\t",tab11,"\t",tab12,"\n")
        
