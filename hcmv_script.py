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
os.chdir(".."
os.system("python3 genbank_parse.py") #logs the # of CDS in genome and saves them into FASTA format txt file (hcmvCompleteGenome.txt)

#build the index of transcriptome with kallisto
os.system("mkdir index")
os.system("time kallisto index -i index/index.idx hcmvCompleteTranscriptome.txt")

#create appropriate output folders for each sample transcriptome
os.system("mkdir outputs")
os.chdir("outputs")
os.system("mkdir SRR5660030.1")
os.system("mkdir SRR5660033.1")
os.system("mkdir SRR5660044.1")
os.system("mkdir SRR5660045.1")
os.chdir("..")

#quantification of each sample with index created
os.system("time kallisto quant -i index/index.idx -o outputs/SRR5660030.1 -b 30 -t 2 SRR5660030.1_1.fastq SRR5660030.1_2.fastq")
os.system("time kallisto quant -i index/index.idx -o outputs/SRR5660033.1 -b 30 -t 2 SRR5660033.1_1.fastq SRR5660033.1_2.fastq")
os.system("time kallisto quant -i index/index.idx -o outputs/SRR5660044.1 -b 30 -t 2 SRR5660044.1_1.fastq SRR5660044.1_2.fastq")
os.system("time kallisto quant -i index/index.idx -o outputs/SRR5660045.1 -b 30 -t 2 SRR5660045.1_1.fastq SRR5660045.1_2.fastq")

#using R package sleuth to find differentially expressed genes between timepoints
os.system("Rscript sleuth_degs.R")
logging.basicConfig(filename="miniProject.log", level=logging.INFO) #loop through and log each file of the deg results
fdr05_results = open('fdr05_results.txt', 'r')
for line in fdr05_results:
    logging.info(line)

#using bowtie to map the sample reads to the hcmv transcriptome
os.system("bowtie2-build hcmvCompleteTranscriptome.txt hcmv")
os.system("bowtie2 -x hcmv -1 SRR5660030.1_1.fastq -2 SRR5660030.1_2.fastq -S hcmvmap.sam --al-conc SRR5660030.1_mapped.fq")
os.system("bowtie2 -x hcmv -1 SRR5660033.1_1.fastq -2 SRR5660033.1_2.fastq -S hcmvmap.sam --al-conc SRR5660033.1_mapped.fq")
os.system("bowtie2 -x hcmv -1 SRR5660044.1_1.fastq -2 SRR5660044.1_2.fastq -S hcmvmap.sam --al-conc SRR5660044.1_mapped.fq")
os.system("bowtie2 -x hcmv -1 SRR5660045.1_1.fastq -2 SRR5660045.1_2.fastq -S hcmvmap.sam --al-conc SRR5660045.1_mapped.fq")


