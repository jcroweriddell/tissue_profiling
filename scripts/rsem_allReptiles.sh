#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to) 
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 32            # number of cores (here uses 2)
#SBATCH --time=3-00:00:00 # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --mem=64GB      # memory pool for all cores (here set to 32 GB)
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jimmy.breen@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: script generate gene expression profiles for each sea snake tissue 
## 	(eye, testis, liver, heart, VMO, Atenuis tail x 2 and body, Alaevis tail x 3 and body), 
## 	to be run using phoenix HPC 

TRIMDIR=/data/biohub/160801_Jenna_seaSnakeRNAseq/RNAseq/Trimmed
TRANSDIR=/data/biohub/Refs/Reptiles/gff3
OUTDIR=/home/a1650598/fastDir/tissue_profiling/All_Reptiles
db_list=/fast/users/a1650598/tissue_profiling/Reptile_rsemDB.list

# Load our modules
module load RSEM/1.2.25-foss-2015b
module load Bowtie2/Bowtie2/2.2.6-foss-2015b


# Prepare reference sequence
for gff in ${TRANSDIR}/*_genomic.gff 
 do
	rsem-prepare-reference --gff3 ${gff} ${TRANSDIR}/$(basename $gff _genomic.gff)_rna.fna \
		 --bowtie2 ${TRANSDIR}/$(basename $gff _genomic.gff)_bt2
done

# Calculate expression levels
while read line; do
	for FQGZ in ${TRIMDIR}/*_truncated1*.fastq.gz
 	 do
		rsem-calculate-expression --bowtie2 -p 8 --paired-end  \
			${FQGZ} ${FQGZ/truncated1/truncated2} \
			${TRANSDIR}/"${line}"_bt2 \
			${OUTDIR}/$(basename $FQGZ _truncated1.fastq.gz)_rsem
	done
done < ${db_list}







