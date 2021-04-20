import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 

argv = OptionParser()

argv.add_option("-I", "--Input", action = "store", dest = "file", type ="string",
                   help = "rpsblast file complete")
argv.add_option("-f","--first", action="store", dest="firstid", type="string",
                    help = "first id in eggnog out" )

(argumentos, palha_que_nao_interessa) = argv.parse_args()


rps_open = open (argumentos.file, 'r')
rps_all_lines = rps_open.readlines()

namefirst = argumentos.firstid

sugar_transport_system=["COG1129", "COG1175", "COG1879", "COG3839", "COG4158"]
glycerol3phosphate_transport=["COG0395", "COG1653"]
tricarboxylate_transporter=["COG3181"]

aminoacid_transport=["COG0410", "COG0411", "COG0683", "COG4177"] 
putrescine_transport=["COG1176", "COG1177"]
glycine_betaine_transport=["COG1125", "COG1174", "COG2113", "COG4175", "COG4176"]

bicarbonate_transport=["COG0600", "COG0715", "COG1116"]
Fe3_transport=["COG0609", "COG0614", "COG1178", "COG1840", "COG3842", "COG4594"]

uncharacterized_transport=["COG1101", "COG2984", "COG3218", "COG3225", "COG3683", "COG3694", "COG4120", "COG4132", "COG4134", "COG4135", "COG4136", "COG4137", "COG4152", "COG4178", "COG4586", "COG4587", "COG4590", "COG4674"]
transport=["COG0226", "COG0310", "COG0390", "COG0395", "COG0410", "COG0411", "COG0444", "COG0555", "COG0559", "COG0573", "COG0577", "COG0581", "COG0600", "COG0601", "COG0609", "COG0614", "COG0683", "COG0715", "COG0725", "COG0747", "COG0755", "COG0765", "COG0842", "COG1101", "COG1108", "COG1116", "COG1117", "COG1118", "COG1119", "COG1120", "COG1121", "COG1123", "COG1124", "COG1125", "COG1126", "COG1129", "COG1131", "COG1132", "COG1134", "COG1135", "COG1172", "COG1173", "COG1174", "COG1175", "COG1176", "COG1177", "COG1178", "COG1277", "COG1464", "COG1613", "COG1653", "COG1732", "COG1840", "COG1879", "COG1930", "COG2011", "COG2113", "COG2386", "COG2984", "COG2998", "COG3127", "COG3218", "COG3221", "COG3225", "COG3638", "COG3639", "COG3683", "COG3694", "COG3833", "COG3839", "COG3840", "COG3842", "COG4107", "COG4120", "COG4132", "COG4133", "COG4134", "COG4135" "COG4136", "COG4137", "COG4138", "COG4139", "COG4143", "COG4148", "COG4149", "COG4150", "COG4152", "COG4158", "COG4160", "COG4161", "COG4166", "COG4172", "COG4174", "COG4175", "COG4176", "COG4177", "COG4178", "COG4181", "COG4208", "COG4209", "COG4211", "COG4213", "COG4214", "COG4215", "COG4239", "COG4473", "COG4521", "COG4525", "COG4531", "COG4555", "COG4558", "COG4559", "COG4586", "COG4587", "COG4590", "COG4591", "COG4592", "COG4594", "COG4597", "COG4598", "COG4604", "COG4605", "COG4606", "COG4607", "COG4608", "COG4618", "COG4662", "COG4674", "COG4721", "COG4779", "COG4985", "COG4986", "COG4987", "COG4988", "COG5265"]

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

print("Genome\tsugar_transport_system\tglycerol3phosphate_transport\ttricarboxylate_transporter\taminoacid_transport\tspermidine_putrescine_transport\tproline_glycine_betaine_transport\tnitrate_sulfonate_bicarbonate_transport\tFe3+_transport\tuncharacterized_transport\tABC-type_transport")

Query = 0
KO = []
for rps_lines in rps_all_lines:

    rps_elements = rps_lines.split(" ")
    #print("teste", rps_elements[0])

    if re.match(r"Query=", rps_elements[0]):
        #print(Query, "\t", KO, "\n")
        KO = []
        Query = rps_elements[1]
        Query = Query.rstrip('\n')
        name = Query.split(".")
        #print ("\n", Query)
    if re.match(r"CDD:", rps_elements[0]):
        GOone = rps_elements[2]
        GOone = GOone.rstrip('\n')
        GOone = GOone.split(",")
        #print("------", GOone, "\n")
        
        combine1 = set(GOone) & set(sugar_transport_system)
        tamanho1 = len(combine1)

        combine2 = set(GOone) & set(glycerol3phosphate_transport)
        tamanho2 = len(combine2)

        combine3 = set(GOone) & set(tricarboxylate_transporter)
        tamanho3 = len(combine3)

        combine4 = set(GOone) & set(aminoacid_transport)
        tamanho4 = len(combine4)

        combine5 = set(GOone) & set(putrescine_transport)
        tamanho5 = len(combine5)

        combine6 = set(GOone) & set(glycine_betaine_transport)
        tamanho6 = len(combine6)

        combine7 = set(GOone) & set(bicarbonate_transport)
        tamanho7 = len(combine7)

        combine8 = set(GOone) & set(Fe3_transport)
        tamanho8 = len(combine8)

        combine9 = set(GOone) & set(uncharacterized_transport)
        tamanho9 = len(combine9)

        combine10 = set(GOone) & set(transport)
        tamanho10 = len(combine10)

        
        if namefirst == name[0]:
        print("-----tem", combine, "\t", tamanho, "\n")
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
            
        
        else:
            print(namefirst,"\t",tab1,"\t",tab2,"\t",tab3,"\t",tab4,"\t",tab5,"\t",tab6,"\t",tab7,"\t",tab8,"\t",tab9,"\t",tab10)
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
    

print(namefirst,"\t",tab1,"\t",tab2,"\t",tab3,"\t",tab4,"\t",tab5,"\t",tab6,"\t",tab7,"\t",tab8,"\t",tab9,"\t",tab10)
     

