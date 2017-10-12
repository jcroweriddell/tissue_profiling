#!/bin/bash -l

# Why dont we use salmon for the quantification instead of RSEM. 
# (Less wasted alignments and extra files)

TRIMDIR=/data/biohub/160801_Jenna_seaSnakeRNAseq/RNAseq/Trimmed
TRANSDIR=/data/biohub/Refs/Reptiles/gff3
OUTDIR=${HOME}/jbreen/BioinfoHub/170824_Jenna_seaSnake/tissue_profiling/salmon_test
THREADS=16

# Prepare reference sequence
for fna in ${TRANSDIR}/*_rna.fna.gz
 do
        salmon index -t $fna \
         -i ${TRANSDIR}/$(basename $fna _rna.fna.gz)_salmon_idx \
         -p ${THREADS}
done

# Calculate expression levels
for FQGZ in ${TRIMDIR}/*_truncated1*.fastq.gz
 do
    for idx in ${TRANSDIR}/*_salmon_idx
     do
        salmon quant -p ${THREADS} \
	    -1 ${FQGZ} -2 ${FQGZ/truncated1/truncated2} \
            -i ${idx} -o ${OUTDIR}
    done
done
