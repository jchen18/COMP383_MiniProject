# COMP383_MiniProject
## Jessie Chen

##  Pipeline for analyzing transcriptomes of HCMV post infection
### Brief Background
Human herpesvirus (abbreviated to HCMV) is a beta-herpesvirus infection that remains in the human system after initial conception. This latency can last a lifetime, but little is known on the transcriptional mechanisms of the herpesvirus latency. The Cheng lab at the Unversity of Arizona produced the transcriptomes of HCMV post infection ([Cheng et al. 2017](https://pubmed.ncbi.nlm.nih.gov/29158406/)). Here we extracted the transcriptomes of 2 patient donors 2- and 6-days post-infection (dpi) and executed a series of analyses to quantify the expression of the transcriptomes as well as compare the the transcriptome reads to other publicly available viral strains. 

### Packages to have installed before running
Python: Biopython
<br>
R: sleuth
<br>
Unix: kallisto, blast+, SPAdes, bowtie, SRA toolkit

### Files provided in repo
- test_data folder: directory with subset of full input reads
- README.md: description of script and manual on usage
- betaherp_sequences.fasta.gz: zip file database of sequences of *Betaherpesvirinae* subfamily from NCBI 
- genbank_parse.py: python script used to create transcriptome of HCMV
- \***hcmv_script.py: main python script for pipeline***
- hcmv_table.txt: text file containing samples, timestamp, and paths to kallisto results for sleuth
- \***miniProject.log: requested output from running pipeline with all input reads***
- sleuth_degs.R: r script to perform statistical analysis of kallisto output

### How to use
Once you have cloned the repository into your directory using this command,
<br>
`git clone https://github.com/jchen18/COMP383_MiniProject.git`
<br>
Move into the COMP383_MiniProject directory created,
<br>
`cd COMP383_MiniProject`
<br>
To run the script, you will need to run the hcmv_script.py file, 
<br>
`python3 hcmv_script.py`
<br>
You will be prompted with the option of either a **full** run using all input reads from the study or a **test** run using the test data available in the test_data directory. The test data consists of 10,000 pair-end reads per sample. 

### Finding the outputs
All output files will be available in the `miniProject_Jessie_Chen` directory and the `miniProject.log` file within the directory produces requested log file information.
