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

#run thru script that parses through NCBI HCMV genome Genbank file
os.system("python3 genbank_parse.py")
