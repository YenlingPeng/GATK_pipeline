#!/bin/bash
#PBS -q ntu192G
#PBS -l select=1:ncpus=5
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N gatk4.1_ANNOVAR
#PBS -o /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/annotation.out
#PBS -e /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/annotation.err
#PBS -M s0890003@gmail.com
#PBS -m e

SAMPLE_ID=NA12878_S38_panel
SCRIPT_PATH=/work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline
INPUT_PATH=/work2/lynn88065/GATK/Outputs/20200721_GATK_pipeline/NA12878_S38_panel
OUTPUT_PATH=/work2/lynn88065/GATK/Outputs/20200721_GATK_pipeline/NA12878_S38_panel_Annotation
CHAIN_FILE=/project/GP1/u3710062/AI_SHARE/reference/Liftover/hg38ToHg19.over.chain.gz
REF_GENOME_PATH=/project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/hg19/ucsc.hg19.fasta

mkdir -p $OUTPUT_PATH/gatk
mkdir -p $OUTPUT_PATH/mutect2

$SCRIPT_PATH/run_annotation.sh $SAMPLE_ID \
        $INPUT_PATH/${SAMPLE_ID}.bwamem.marked.recal.pass1.indexed.filtered.haplotype.SnpIndel.vcf.gz \
        $OUTPUT_PATH/gatk \
        $CHAIN_FILE \
        $REF_GENOME_PATH 

$SCRIPT_PATH/run_annotation.sh $SAMPLE_ID \
        $INPUT_PATH/${SAMPLE_ID}.bwamem.mutect2.vcf.gz \
        $OUTPUT_PATH/mutect2 \
        $CHAIN_FILE \
        $REF_GENOME_PATH 
