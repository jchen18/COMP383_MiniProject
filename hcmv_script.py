import os
#creating outfile folder
os.system("mkdir miniProject_Jessie_Chen")
os.chdir("miniProject_Jessie_Chen")

#retrieving HCMV transcriptomes of patients
#Donor 1 (2dpi)
os.system("wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1")
#Donor 1 (6dpi)
os.system("wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1")
#Donor 3 (2dpi)
os.system("wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1")
#Donor 3 (6dpi)
os.system("wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1")

#convert to paired-end fastq files
os.system("fastq-dump -I --split-files SRR5660030.1")
os.system("fastq-dump -I --split-files SRR5660033.1")
os.system("fastq-dump -I --split-files SRR5660044.1")
os.system("fastq-dump -I --split-files SRR5660045.1")

#accessing the HCMV Genbank ref genome (NCBI accession EF999921)
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "jessiechen0513@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb", retmode="text") #retrieving in Genbank format
record = SeqIO.read(handle, "genbank")
handle.close()

logging.basicConfig(filename="miniProject.log", level=logging.INFO) #opening file for logging info
num_CDS = 0
for feature in record.features: #loop all features and counting ones that are CDS features
        if feature.type  == "CDS":
                num_CDS += 1
logging.info('The HCMV genome (EF999921) has ' + str(num_CDS) + ' CDS.') #write to log file
