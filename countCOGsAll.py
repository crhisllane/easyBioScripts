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
argv.add_option("-c","--coglst", action="store", dest="cogs", type="string",
                    help = "cog table" )
argv.add_option("-f","--first", action="store", dest="firstid", type="string",
                    help = "first id in eggnog out" )

(argumentos, palha_que_nao_interessa) = argv.parse_args()


rps_open = open (argumentos.file, 'r')
rps_all_lines = rps_open.readlines()

cogs_open = open (argumentos.cogs, 'r')
cogs_all_lines = cogs_open.readlines()

namefirst = argumentos.firstid


Query = 0

tamanhoCOGs = len(cogs_all_lines)
cogs_lst=[]
cogs_name=[]
for cogs_line in cogs_all_lines:
    cogs_elements = cogs_line.split("\t")
    cogs_lst.append(cogs_elements[0])
    cogs_completeName = cogs_elements[0] + " " + cogs_elements[2]
    cogs_name.append(cogs_completeName)

cogs_count = [0] * (len(cogs_lst))


print("Genome\t", cogs_name)


for rps_lines in rps_all_lines:

    rps_elements = rps_lines.split(" ")
    

    if re.match(r"Query=", rps_elements[0]):
        #print(Query, "\t", KO, "\n")
        Query = rps_elements[1]
        Query = Query.rstrip('\n')
        name = Query.split("_")
        #print ("\n", Query, '\t---', name[0])
    if re.match(r"^CDD:", rps_elements[0]):
        COGone = rps_elements[2]
        COGone = COGone.rstrip('\n')
        COGone = COGone.split(",")
        #print("------", COGone, "\n")

        if namefirst == name[0]:   
            for i in range(len(cogs_lst)):
                if COGone[0] == cogs_lst[i]:
                    cogs_count[i] = cogs_count[i] + 1
        else:
            print (namefirst, "\t", cogs_count)
            namefirst = name[0]
            cogs_count = [0] * (len(cogs_lst))          
            for i in range(len(cogs_lst)):
                if COGone[0] == cogs_lst[i]:
                    cogs_count[i] = cogs_count[i] + 1

print (namefirst, "\t", cogs_count)
          