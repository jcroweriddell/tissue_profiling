#!/bin/bash
#SBATCH -p batch                                            # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes (due to the nature of sequential processing\
, here uses single node)
#SBATCH -n 32                                              # number of cores (here uses 2)
#SBATCH --time=24:00:00                                    # time allocation, which has the format (D-HH:MM), here set t\
o 1 hour
#SBATCH --mem=80GB                                         # memory pool for all cores (here set to 32 GB)
# Notification configuration
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sen\
t when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be se\
nt when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

# Executing script (Example here is sequential script and you have to select suitable compiler for your case.)          \
                               # bash script used here for demonstration purpose, you should select proper compiler for your needs\
module load Trinity

Trinity --monitoring --seqType fq --SS_lib_type FR \
--left /home/a1662801/seaSnake_reads/ALA_reads/4_HAJ15ADXX_TAGCTT_L001_R1.fastq.gz,/home/a1662801/seaSnake_reads/ALA_reads/5_HAJ15ADXX_GGCTAC_L001_R1.fastq.gz,home/a1662801/seaSnake_reads/ALA_reads/ALA_L1_1.fq.gz,/home/a1662801/seaSnake_reads/ALA_reads/C9FMBANXX-2-JUB_A_LAEVIS_ANT]TAIL_RIGHT__ACAGTG_L002_R1_001.fastq.gz,/home/a1662801/seaSnake_reads/ALA_reads/C9FMBANXX-2-JUV_A_LAEVIS_VENTRAL_TAIL_ACTTGA_L002_R1_001.fastq.gz,/home/a1662801/seaSnake_reads/ALA_reads/C9FMBANXX-2-J_A_LAEVIS__BODY_DORSAL_RIGHT_CGATGT_L002_R1_001.fastq.gz \
--right /home/a1662801/seaSnake_reads/ALA_reads/4_HAJ15ADXX_TAGCTT_L001_R2.fastq.gz,/home/a1662801/seaSnake_reads/ALA_reads/5_HAJ15ADXX_GGCTAC_L001_R2.fastq.gz,/home/a1662801/seaSnake_reads/ALA_reads/ALA_L1_2.fq.gz,/home/a1662801/seaSnake_reads/ALA_reads/C9FMBANXX-2-JUB_A_LAEVIS_ANT]TAIL_RIGHT__ACAGTG_L002_R2_001.fastq.gz,/home/a1662801/seaSnake_reads/ALA_reads/C9FMBANXX-2-JUV_A_LAEVIS_VENTRAL_TAIL_ACTTGA_L002_R2_001.fastq.gz,/home/a1662801/seaSnake_reads/ALA_reads/C9FMBANXX-2-J_A_LAEVIS__BODY_DORSAL_RIGHT_CGATGT_L002_R2_001.fastq.gz \
--CPU 32 --trimmomatic --max_memory 80G --normalize_reads
