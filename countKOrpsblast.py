import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 

argv = OptionParser()

argv.add_option("-I", "--Input", action = "store", dest = "file", type ="string",
                   help = "cluster cdhit output file")
argv.add_option("-l","--list", action="store", dest="list", type="string",
                    help = "cdhit output fasta file" )
(argumentos, palha_que_nao_interessa) = argv.parse_args()

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

def process_count(rpsfile):
    print("file", rpsfile, "\n")
  
    rps_open = open (rpsfile, 'r')
    rps_all_lines = rps_open.readlines()
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
            name = Query.split("_")
            #print ("\n", Query)
        if re.match(r"CDD:", rps_elements[0]):
            KOone = rps_elements[2]
            KOone = KOone.rstrip('\n')
            KO.append(KOone)
            print (name[0], "\t", KOone, "\n")
            #if KOone in KOs_line:
                #print (Query, "\t", KOone, "\n")
    
    #print(Query, "\t", KO, "\n")
        

            

path = os.getcwd()
path = str(path) + '/*.rpsblast' 
files_rpsblast =  glob.glob(path) 
    
with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_count, files_rpsblast)