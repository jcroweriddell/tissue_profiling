#!/bin/bash
#SBATCH -p batch                                            # partition (this is the queue your job will be added to) 
#SBATCH -N 1                                               # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8                                              # number of cores (here uses 2)
#SBATCH --time=24:00:00                                    # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --mem=32GB                                         # memory pool for all cores (here set to 32 GB)

# Notification configuration 
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

# Executing script (Example here is sequential script and you have to select suitable compiler for your case.)                                         # bash script used here for demonstration purpose, you should select proper compiler for your needs 

## USAGE: script generate gene expression profiles for each sea snake tissue (eye, testis, liver, heart, VMO, Atenuis tail x 2 and body, Alaevis tail x 3 and body), to be run using phoenix HPC 
module load HISTAT2
module load HTSeq

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

mkdir -p $output_dir

# Run fastqc on each read pairs
echo "starting fastqc"
fastqc -o $output_dir/ $data_dir/*.fastq.gz 
echo "finished fastqc"

for FQGZ in $data_dir/*_R1*.fastq.gz
 do
    # Everything indented (and before the "done") will be run
 
    # Get our raw data and trim
    echo "Starting Trimming of "$FQGZ" "
    AdapterRemoval --file1 $FQGZ --file2 ${FQGZ/R1/R2} \
	 --output1  $working/output/"$(basename $FQGZ _R1.fastq.gz)".trimed1.fq.gz \
	 --output2  $working/output/"$(basename $FQGZ _R1.fastq.gz)".trimed2.fq.gz \
	 --gzip  --trimns --trimqualities --minlength 20
    echo "Finishing Trimming of "$FQGZ" "
    
    # Align our trimmed reads against our reference
    echo "Starting Alignment of "$FQGZ" "
    hisat2 -x $genome_prefix \
	-1 $working/output/"$(basename $FQGZ _R1.fastq.gz)".trimed1.fq.gz \
	-2 $working/output/"$(basename $FQGZ _R1.fastq.gz)".trimed2.fq.gz | \
	samtools view -bS - > $working/output/"$(basename $FQGZ _R1.fastq.gz)".hisat2_12tissues.bam
    echo "Finishing Alignment of "$FQGZ" "
    
done

## Align reads to ref sequence
# Input variables: forward and reverse reads
R1=$1
R2=$2
reference_index=/home/a1662801/ref_seq/tissues_12_ref

# Run hisat2 on samples
# - Change -p option for the speed of your computer
hisat2 -p 2 -x $reference_index \
  -1 $R1 -2 $R2 | \
  samtools view -bS -F4 - > "$(basename $R1 _L001_R1.fastq.gz)"_mapped_hisat2_12tissues.bam

# Give me some statistics
samtools flagstat "$(basename $R1 _L001_R1.fastq.gz)"_mapped_hisat2_12tissues.bam

## Generate read counts (FPKM) for each sea snake tissue
# For paired-end data, htseq-count needs alignment to be sorted either by read name or by alignment position, use samtools to sort by position (default)
samtools sort "$(basename $R1 _L001_R1.fastq.gz)"_mapped_hisat2_12tissues.bam

htseq-count -s yes -r pos  "$(basename $R1 _L001_F1.fastq.gz)"_mapped_hitsat2_12tissues.bam <gff_file> #need an annotation file here, not sure how to do that with a custom annotation of the de novo contigs


