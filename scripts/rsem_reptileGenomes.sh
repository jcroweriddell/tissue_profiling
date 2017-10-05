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
