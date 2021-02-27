#accessing the HCMV Genbank genome (NCBI accession EF999921)
import logging
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "jessiechen0513@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb", retmode="text") #retrieving in Ge>
record = SeqIO.read(handle, "genbank")
handle.close()

logging.basicConfig(filename="miniProject.log", level=logging.INFO) #opening file for logging info
num_CDS = 0
for feature in record.features: #loop all features and counting ones that are CDS features
        if feature.type  == "CDS":
                num_CDS += 1
logging.info('The HCMV genome (EF999921) has ' + str(num_CDS) + ' CDS.') #write to log file

#writing out the full sequence from FASTA into separate file
outfile = open("hcmv_complete_genome.txt", "w")
outfile.write(str(record.seq))
outfile.close()
