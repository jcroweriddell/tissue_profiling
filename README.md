
# Tissue profiling project

2017-08-24

Jenna Crowe-Riddell
Jimmy Breen (jimmymbreen@gmail.com)
Alastair Luddington

## To Do List

### Jimmy
- [x] Run BLAST search of 20 genes on King Cobra genome and identify IDs
- [x] Run same BLAST search of 20 genes on tissue-specific transcriptomes and identify IDs or contigs
- [x] Re-run RNAseq analysis on King Cobra (might have done this already) and create count table - I got TPM and RPKM as well
- [ ] Re-run RNAseq analysis on SeaSnake assembky (if possible?) and create count table (needs Alastair's gene annotation)
- [x] Run Salmon on all libraries (see below)

### Jenna
- [x] Provide Jimmy with list of 20 visual genes to BLAST
Note: Jenna has uploaded excel sheet with list of important genes in visual pathway, and uploaded a fasta file with sequences of top genes for Jimmy to use, the list includes visual opsins, phototransduction, retinoid metabolism. Missing: non-visual genes
- [X] Explore FPKM count table in R across sample (MDS plot, PCA plot, heatmap) and within samples (frequency of counts per gene). This will be useful for after Jimmy's BLASTx search
- [ ] Convert TPM to FPKM for new Salmon quant analysis, then execute R scripts for joining GeneNames, clustering analysis and heatmap for visual genes

## Tissues to Analyse

Individual 1: _A. laevis_ juvenile:
- Juvenile body skin (2-J_A_LAEVIS__BODY_DORSAL_RIGHT_CGATGT)
- Ventral Tail (2-JUV_A_LAEVIS_VENTRAL_TAIL_ACTTGA)
- Anterior Tail (2-JUB_A_LAEVIS_ANT]TAIL_RIGHT__ACAGTG)

Individual 2: _A. laevis_ KLS0468:
- VMO (4_HAJ15ADXX_TAGCTT)
- Tail A2 (5_HAJ15ADXX_GGCTAC)

Extras (potentially for Genome paper analysis):
1. _H. major_ KLS0460, Liver (2_HAJ15ADXX_ACTTGA)
2. _H. major_ KLS0460, Heart (3_HAJ15ADXX_GATCAG)

## Salmon Analysis

RSEM seemed to be fairly inefficient when it came to actually identifying the FPKM or TPM of the target genes for each genome. Not only do you need to prepare the genome using the gff3 and transcript fna files, you also need to map the fastq files to the genome to identify the level of quantification. So instead of doing that process, here I have used the "pseudo-alignment" program called [`Salmon`](https://combine-lab.github.io/salmon), which quantifies the level of expression to a reference genome _without_ running the actual alignment process. `Salmon` is one of a number of programs (`Sailfish` and `Kallisto` being two others), that offer this approach where we match k-mers (segments of your FASTQ reads of length `k`), to an index made up of target transcripts.

For this process I ran the following workflow:

INPUT: PAIRED FASTQ Files

1. Download Reptile transcriptomes (all annotated genes/transcripts)
2. Create a `salmon` index
3. Quantify FASTQ libraries against `salmon` index
4. Rename the output files (quant.sf)
5. Convert counts to FPKM (**TO DO**)

OUTPUT: Each "quant.sf" file contains information on the TPM and gene-length for each gene

For TPM to FPKM conversion in R:

```
tpmToFpkm <- function(tpm, geneLength){
	count <- apply(tpm, 2, function(e){ e/1000*geneLength })
	fpkm <- apply(count, 2, function(e){ e*1000*1000000/sum(e)/geneLength })
	return(fpkm)
}
```

