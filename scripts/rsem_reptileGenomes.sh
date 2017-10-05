#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to) 
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores (here uses 2)
#SBATCH --time=24:00:00 # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: script generate gene expression profiles for each snake tissue
## 	(14 TISSUES: eye, VMO x 3, skin x 7, liver, heart, testis) by aligning to different reptile genomes (lizards, snakes, chicken, alligator)

# Load our modules
# module load HISAT2/2.0.5-foss-2016uofa
# module load HTSeq/0.6.1p1-intel-2015c-Python-2.7.11
# module load Java/1.8.0_71
# module load RSEM/1.2.25-foss-2015b
# module load Bowtie2/2.2.9-GCC-5.3.0-binutils-2.25

## Script for using RSEM to calculate read counts using different reptile genomes

## Set working and output directories
working= /home/jenna
data_dir= /bg/raw_reads/kls/Alaevis_reads
output_dir= $working/output

## Download reptile genomes from ncbi
mkdir $working/ref_genomes/

# Protobothrops viper
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/527/695/GCF_001527695.2_P.Mucros_1.0/GCF_001527695.2_P.Mucros_1.0_genomic.fna.gz

# Python
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/186/305/GCF_000186305.1_Python_molurus_bivittatus-5.0.2/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.fna.gz

# Gecko
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/447/785/GCF_001447785.1_Gekko_japonicus_V1.1/GCF_001447785.1_Gekko_japonicus_V1.1_genomic.fna.gz

#Pogona bearded dragon
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/067/755/GCF_900067755.1_pvi1.1/GCF_900067755.1_pvi1.1_genomic.fna.gz

# Anolis green lizard
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/090/745/GCF_000090745.1_AnoCar2.0/GCF_000090745.1_AnoCar2.0_genomic.fna.gz

# Alligator
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/281/125/GCF_000281125.3_ASM28112v4/GCF_000281125.3_ASM28112v4_genomic.fna.gz

#Gallus chicken
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.4_Gallus_gallus-5.0/GCF_000002315.4_Gallus_gallus-5.0_genomic.fna.gz

## Build reference index for genome

# Using bowtie2
bowtie2-build $working/ref_genomes/GCF_001527695.2_P.Mucros_1.0_genomic.fna.gz $working/ref_genomes/Protobothrops

genome_prefix= $working/ref_genomes/Protobothrops

## Script for fastQC reads, align to ref, and adapter removal

# Run fastqc on each read pairs
echo "Starting fastqc"
nohup fastqc -o $output_dir/ $data_dir/*.fastq.gz

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
    echo "Starting to calculate expression levels of "$FQGZ" "
    rsem-calculate-expression --bowtie2 -p 8 --paired-end \
			      $output_dir/"$(basename $FQGZ _R1.fastq.gz)".trimed1.fq.gz \
			      $output_dir/"$(basename $FQGZ _R1.fastq.gz)".trimed2.fq.gz \
			      $genome_prefix \
			      $output_dir/"$(basename $FQGZ _R1.fastq.gz)".rsem_Alaevis
    echo "Finished calculating expression levels of "$FQGZ" "


  done
