#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to) 
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores (here uses 2)
#SBATCH --time=24:00:00 # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: script generate gene expression profiles for each sea snake tissue 
## 	(eye, testis, liver, heart, VMO, Atenuis tail x 2 and body, Alaevis tail x 3 and body), 
## 	to be run using phoenix HPC 


# Load our modules
module load RSEM/1.2.25-foss-2015b
module load Bowtie2/Bowtie2/2.2.6-foss-2015b


## Script for using RSEM to generate FKPM from raw read counts

# Prepare reference sequence
rsem-prepare-reference /home/a1662801/ref_seq/all_tissues_contigs.fa --bowtie2 \
/home/a1662801/ref_seq/all_tissues

# Calculate expression levels

rsem-calculate-expression --bowtie2 -p 8 --paired-end \
 /data/biohub/160801_Jenna_seaSnakeRNAseq/RNAseq/Trimmed/4_HAJ15ADXX_TAGCTT_truncated1.fastq.gz \ # Trimed reads 1
 /data/biohub/160801_Jenna_seaSnakeRNAseq/RNAseq/Trimmed/4_HAJ15ADXX_TAGCTT_truncated2.fastq.gz \ # Trimed reads 2
 /home/a1662801/all_tissues /home/a1662801/rsemResults/4_VMO







