#!/bin/bash
#SBATCH -p batch   # Queue
#SBATCH -N 1    # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8    # number of cores (here uses 2)
#SBATCH --time=00:02:00    # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=25GB     # memory pool for all cores 
#SBATCH -o /home/a1662801/slurm/visGenes_%j.out
#SBATCH -e /home/a1662801/slurm/visGenes_%j.err
#SBATCH --mail-type=END     
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

# Script to extract visual gene names 'XM_0000000' from ref genome fasta file
# To be run on Phoenix

# Set variables
REFGENOME=GCF_001527695.2_P.Mucros_1.0_rna.fna.gz
DATA=/data/biohub/Refs/Reptiles/gff3
GENE_LIST=visGenelist.tsv

# Extract vis genes from genome
echo "Searching for visual gene names in $ref_genome"

for FQGZ in ${DATA}/*_rna.fna.gz
	do
	while read line; do zgrep ">" GCF_001527695.2_P.Mucros_1.0_rna.fna.gz | grep "$line"; done < $GENE_LIST
done


echo "Finished searching in $ref_genome"
