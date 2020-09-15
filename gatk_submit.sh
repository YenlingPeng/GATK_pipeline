#!/bin/bash
#PBS -q ntu192G
#PBS -l select=1:ncpus=10
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N gatk4.1_mutect2
#PBS -o /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/gatk.out
#PBS -e /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/gatk.err
#PBS -M s0890003@gmail.com
#PBS -m e

date="20200721"
clinical="GATK"
disease="pipeline"
sample_name="NA12878_S38_panel"
OUTPUT_PATH="/work2/lynn88065/GATK/Outputs/${date}_${clinical}_${disease}/${sample_name}"
SCRIPT_PATH=/work2/lynn88065/GATK/Inputs/${date}_${clinical}_${disease}/
INPUT_PATH=/project/GP1/u3710062/AI_SHARE/Temp/20191018_DataToYL/

$SCRIPT_PATH/run_gatk.sh $sample_name \
        $INPUT_PAHT/NA12878-D702-D508_S38_L001_R1_001.fastq.gz \
        $INPUT_PAHT/NA12878-D702-D508_S38_L001_R2_001.fastq.gz \
        "@RG\tID:${sample_name}\tLB:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA" \
        $OUTPUT_PATH \
        /project/GP1/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
        /project/GP1/reference/Homo_sapiens/GATK/hg38/dbsnp_146.hg38.vcf.gz \
        /project/GP1/reference/Homo_sapiens/GATK/Mutect2/af-only-gnomad.hg38.vcf.gz \
        /project/GP1/reference/Homo_sapiens/GATK/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
        /project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/2.8/hg38/PoN/somatic-hg38-1000g_pon.hg38.vcf.gz \
        /work2/lynn88065/GATK/Inputs/${date}_${clinical}_${disease}/Panel_targets_EGFR_v1.interval_list \
        /work2/lynn88065/GATK/Inputs/${date}_${clinical}_${disease}/Panel_targets_EGFR_v1.interval_list
