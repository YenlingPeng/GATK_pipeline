#!/bin/bash
#PBS -q ntu192G
#PBS -l select=1:ncpus=1
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N liftover_bed
#PBS -o /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/liftover.out
#PBS -e /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/liftover.err
#PBS -M s0890003@gmail.com
#PBS -m e

# UCSC provides tools to convert BED file from one genome assembly to another

# Basic information and path
# NOTE: Use the 'chr' before each chromosome name
liftOver="/work2/lynn88065/Software/Download/LiftOver/liftOver"
INPUT_BED="/project/GP1/u3710062/AI_SHARE/Temp/20191018_DataToYL/EGFR_v1.bed"
CHAIN_FILE="/project/GP1/u3710062/AI_SHARE/reference/Liftover/hg19ToHg38.over.chain.gz"
OUTPUT_PATH="/work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline"

# Run liftOver:
# unlifted.bed file will contain all genome positions that cannot be lifted
$liftOver $INPUT_BED $CHAIN_FILE $OUTPUT_PATH/EGFR_v1_hg38.bed unlifted.bed
