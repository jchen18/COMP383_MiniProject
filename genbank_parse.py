#accessing the HCMV Genbank genome (NCBI accession EF999921)
import logging
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "jessiechen0513@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb", retmode="text") #retrieving in Ge>
record = SeqIO.read(handle, "genbank")
handle.close()

outfile = open("hcmvCompleteGenome.txt", "w") #file to save the CDS regions only

logging.basicConfig(filename="miniProject.log", level=logging.INFO) #opening file for logging info
num_CDS = 0
for feature in record.features: #loop all features and counting ones that are CDS features
        if feature.type  == "CDS":
                num_CDS += 1
                name = str(feature.qualifiers['protein_id']) #saving the protein id as the name for seq
                seq = feature.extract(record.seq) #extract out seq
                outfile.write(">" + name + "\n" + str(seq) + "\n" #fasta output in file without line wrapping
logging.info('The HCMV genome (EF999921) has ' + str(num_CDS) + ' CDS.') #write to log file

outfile.close()
