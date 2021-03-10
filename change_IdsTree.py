import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 
import subprocess

argv = OptionParser()

argv.add_option("-t", "--tree", action = "store", dest = "tree", type ="string",
                   help = "tree file from phyml")
argv.add_option("-m","--metadata", action="store", dest="meta", type="string",
                    help = "metadata with old ids and new ids" )

(argumentos, palha_que_nao_interessa) = argv.parse_args()

tree_File = argumentos.tree

namefileout = tree_File + '.renamed' 
outfile = open (namefileout, 'w+')
print(tree_File, namefileout)
comando = ("cp %s %s")%(tree_File, namefileout)
os.system(comando)

metadata_File = open (argumentos.meta, 'r')
metadata_Lines = metadata_File.readlines()

for metadata_Line in metadata_Lines:
      
    meta_Columns = metadata_Line.split("\t")
    print(meta_Columns[0], meta_Columns[1], meta_Columns[2], meta_Columns[3])
    new_name = meta_Columns[2] + meta_Columns[1] + meta_Columns[3]
    #comando2 = ('sed -i -e "s/%s/%s/" %s')%(meta_Columns[0], new_name, namefileout)
    oldID = meta_Columns[0]
    oldID = oldID.rstrip('\n')
    new_name = new_name.rstrip('\n')
    namefileout = namefileout.rstrip('\n')

    #subprocess.call("sed -i 's/^" + meta_Columns[0] + "/" + new_name + "/g' " + namefileout, shell=True)
    subprocess.call("sed -i 's/(" + meta_Columns[0] + "/(" + new_name + "/g' " + namefileout, shell=True)
    subprocess.call("sed -i 's/," + meta_Columns[0] + "/," + new_name + "/g' " + namefileout, shell=True)


    #subprocess.call(comando2, shell=True)

    #os.system(comando2)




