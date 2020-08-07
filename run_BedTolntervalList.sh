#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -q ntu192G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N interval_list
#PBS -M s0890003@gmail.com
#PBS -m e

# This tool requires a sequence dictionary
# provided with the SEQUENCE_DICTIONARY or SD argument
Picard="/pkg/biology/Picard/Picard_v2.18.11/picard.jar"
INPUT_BED="/work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/EGFR_v1_hg38.bed"
OUTPUT_PATH="/work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline"
REFERENCE_DICT="/project/GP1/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.dict"

# Usage
java -Xmx48g -jar $Picard BedToIntervalList \
               I=$INPUT_BED \
               O=$OUTPUT/Panel_targets_EGFR_v1.interval_list \
               SD=$REFERENCE_DICT
