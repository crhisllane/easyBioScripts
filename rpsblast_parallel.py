import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import concurrent.futures
import glob 

argv = OptionParser()

argv.add_option("-I", "--Input", action = "store", dest = "file", type ="string",
                   help = "cdhit fasta")
(argumentos, palha_que_nao_interessa) = argv.parse_args()

#clusterfile = open (argumentos.file, 'r')
#allClusterLines = clusterfile.readlines()
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                next(entry)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

#for clusterline in allClusterLines:
def process_rpsblast(clusterline):
    clusterline = clusterline.rstrip('\n')
    print (clusterline, "\n")
    clusterline_Pure = re.sub(".fasta", "" ,clusterline)
    cmd1 = "rpsblast+ -db db/Cdd -query %s -out %s_spp.rpsblast -num_threads 40 -evalue 0.01  -outfmt '7 qseqid sseqid pident nident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxis frames qseq sseq'"%(clusterline, clusterline_Pure)
    os.system(cmd1)


record_iter = SeqIO.parse(open(argumentos.file), "fasta")
for i, batch in enumerate(batch_iterator(record_iter, 60000)):
    filename = "group_%i.fasta" % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))


path = os.getcwd()
path = str(path) + '/group*' 
files_fasta =  glob.glob(path) 
    
with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_rpsblast, files_fasta)
