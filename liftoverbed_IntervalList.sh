#!/bin/bash
#PBS -q ntu192G
#PBS -l select=1:ncpus=1
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N Liftover&Intervalbed
#PBS -o /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/liftoverInterval.out
#PBS -e /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/liftoverInterval.err
#PBS -M s0890003@gmail.com
#PBS -m e

# Basic information and path
Picard="/pkg/biology/Picard/Picard_v2.18.11/picard.jar"
LIFTOVER="/work2/lynn88065/Software/Download/LiftOver/liftOver"
REFERENCE_DICT="/project/GP1/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.dict"
CHAIN_FILE="/project/GP1/u3710062/AI_SHARE/reference/Liftover/hg19ToHg38.over.chain.gz"
OUTPUT_PATH="/work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline"
INPUT_BED="/project/GP1/u3710062/AI_SHARE/Temp/20191018_DataToYL/EGFR_v1.bed"
Panel_capture="EGFR_v1"
hg="hg38"

# UCSC provides tools to convert BED file from one genome assembly to another
# NOTE: Use the 'chr' before each chromosome name

# ==========================
# 1. Lift genome positions 
# ==========================

# Usage
$LIFTOVER $INPUT_BED $CHAIN_FILE $OUTPUT_PATH/${Panel_capture}_${hg}.bed unlifted.bed


# ===========================
# 2. BedToIntervalList (Picard)
# ===========================

# Usage
java -Xmx48g -jar $Picard BedToIntervalList \
               I=$OUTPUT_PATH/${Panel_capture}_${hg}.bed \
               O=$OUTPUT_PATH/${Panel_capture}_${hg}.interval_list \
               SD=$REFERENCE_DICT
