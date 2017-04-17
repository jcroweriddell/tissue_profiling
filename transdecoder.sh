#!/bin/bash
## Use transdecoder to find (predict) coding regions in transcripts assembled by Trinity

# load TransDecoder
module load TransDecoder

# Run TransDecoder
TransDecoder.Predict -t Trinity.fasta [homology options]


# View ORFs in context of the transcript structures on the de novo assembly
java -jar $GENOMEVIEW/genomeview.jar test.genome.fasta transcripts.bed transcripts.fasta.transdecoder.genome.bed