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
module load R

## VISUALISATION: use R and RSEM to visualise results

# Use R to visualise each sample, order by FPKM and display top 10 hits
R
data = read.table("sample_name.genes.results", header = T, stringsAsFactors = F)
idx = order(data[,"FPKM"], decreasing = T)
data[idx[1:10], c("gene_id", "expected_count", "FPKM")]
quit(save = "sample_nameFPKM", status = 0, runLast = T)


## Use RSEM to visualise read mapping depth to reference genome
rsem-plot-model sample_name sample_name_diagnostic.pdf

## Generate wiggle plot for a particular geneID, stacked wiggle plots RSEM generate stack the expected multi-mapping read depth (shown in red) over the uniquely aligned read depth (shown in black)
rsem-plot-transcript-qiggles --gene-list --show-unique\
sample_name gene_ids.txt geneID_transcript_wiggle.pdf

## Geneate a genome-wide wiggle plot, can load this into IGV with reference sequence
rsem-bam2wig sample_name.genome.sorted.bam sample_name.wig sample_name



