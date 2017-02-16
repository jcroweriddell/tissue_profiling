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
module load HISAT2/2.0.5-foss-2016uofa
module load HTSeq/0.6.1p1-intel-2015c-Python-2.7.11
module load Java/1.8.0_71
module load RSEM/1.2.25-foss-2015b
module load Bowtie2/Bowtie2/2.2.6-foss-2015b


## Build reference index
# Using hisat2
hisat2-build /home/a1662801/ref_seq/tissues_contigs tissues_12_ref
# Using bowtie2
# bowtie-2 /path/to/tissues_contigs tissues_12_ref

## Script for fastQC reads, align to ref, and adapter removal
working=$(pwd)

data_dir=/Volumes/seaSnakeDrive/AGRF_RNAseq
output_dir=$working/output
genome_prefix=/home/a1662801/ref_seq/tissues_12_ref

# You need to make a GFF annotation file on your denovo tissue transcriptome
gff=

mkdir -p $output_dir

# Run fastqc on each read pairs
echo "Starting fastqc"

fastqc -o $output_dir/ $data_dir/*.fastq.gz 

echo "Finished fastqc"

for FQGZ in $data_dir/*R1*.fastq.gz
 do
    # Everything indented (and before the "done") will be run
 
    # Get our raw data and trim
    echo "Starting trimming of "$FQGZ" "
    AdapterRemoval --file1 $FQGZ --file2 ${FQGZ/R1/R2} \
	 --output1  $output_dir/"$(basename $FQGZ _R1.fastq.gz)".trimed1.fq.gz \
	 --output2  $output_dir/"$(basename $FQGZ _R1.fastq.gz)".trimed2.fq.gz \
	 --gzip  --trimns --trimqualities --minlength 20
    echo "Finished trimming of "$FQGZ" "
    
    # Use RSEM to align trimmed reads against our reference and calculate expression levels (TPM and FPKM)
    echo "Starting to calculate expression lvels of "$FQGZ" "
    rsem-calculate-expression --bowtie2 -p 8 --paired-end \
    $output_dir/"$(basename $FQGZ _R1.fastq.gz)".trimed1.fq.gz \
    $output_dir/"$(basename $FQGZ _R1.fastq.gz)".trimed2.fq.gz \
    $genome_prefix \
    $output_dir/"$(basename $FQGZ _R1.fastq.gz)".rsem_12tissues
    echo "Finished calculating expression levels of "$FQGZ" "
  
    # Stats
    samtools flagstat $output_dir/"$(basename $FQGZ _R1.fastq.gz)".rsem_12tissues.bam
    
done

