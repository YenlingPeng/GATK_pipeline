#!/bin/bash

Picard="/pkg/biology/Picard/Picard_v2.18.11/picard.jar"
VT_PATH="/project/GP1/u3710062/AI_SHARE/software/vt-0.57721/vt"
ANNOVAR="/pkg/biology/ANNOVAR/ANNOVAR_20191024/table_annovar.pl"
humandb="/project/GP1/u3710062/AI_SHARE/shared_scripts/ANNOVAR/humandb/"
TODAY=`date +%Y%m%d%H%M`
logfile=./${TODAY}_gatk_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
set -x

# Basic info
SAMPLE_ID=$1
VCF=$2
OUTPUT_PATH=$3
CHAIN_FILE=$4
REF_GENOME_PATH=$5

cd $OUTPUT_PATH

# Liftover VCF (from hg38 to hg19)
java -jar ${Picard} LiftoverVcf \
      I=${VCF} \
      O=${SAMPLE_ID}_hg19.hard.filtered.liftover.vcf \
      CHAIN=$CHAIN_FILE \
      REJECT=${SAMPLE_ID}_hg19.rejected.variants.vcf \
      R=$REF_GENOME_PATH

# VT-normalization                                   
$VT_PATH decompose -s -o ${SAMPLE_ID}_hg19.decomposed.vcf ${SAMPLE_ID}_hg19.hard.filtered.liftover.vcf
$VT_PATH normalize -o ${SAMPLE_ID}_hg19.norm.vcf -r $REF_GENOME_PATH ${SAMPLE_ID}_hg19.decomposed.vcf

# run Annotation
$ANNOVAR ${SAMPLE_ID}_hg19.norm.vcf $humandb -buildver hg19 -out ${SAMPLE_ID} \
        -remove -protocol refGene,TaiwanBiobank-official,TaiwanBiobank993WGS,avsnp150,exac03,exac03nonpsych,clinvar_20190305,gnomad211_exome,gnomad211_genome,cosmic_coding_GRCh37_v91,cosmic_noncoding_GRCh37_v91 -operation gx,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish

# Filtration records (exonic/splicing)
head -n 1 ${SAMPLE_ID}.hg19_multianno.txt > ${SAMPLE_ID}.filtered.annotation.txt
grep -e "exonic" -e "splicing" ${SAMPLE_ID}.hg19_multianno.txt >> ${SAMPLE_ID}.filtered.annotation.txt

