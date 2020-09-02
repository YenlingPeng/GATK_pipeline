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

date="20200721"
clinical="GATK"
disease="pipeline"
sample_name="NA12878_S38_panel"
OUTPUT_PATH="/work2/lynn88065/GATK/Outputs/${date}_${clinical}_${disease}/${sample_name}"

mkdir -p ${OUTPUT_PATH}/${sample_name}_hg19
cd ${OUTPUT_PATH}/${sample_name}_hg19
VCF=${OUTPUT_PATH}/${sample_name}.bwamem.marked.recal.pass1.indexed.filtered.haplotype.SnpIndel.vcf.gz

#liftover from hg38 to hg19 and normalized
Picard="/pkg/biology/Picard/Picard_v2.18.11/picard.jar"
VT_PATH="/project/GP1/u3710062/AI_SHARE/software/vt-0.57721/vt"
TODAY=`date +%Y%m%d%H%M`
logfile=./${TODAY}_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
set -x

java -jar ${Picard} LiftoverVcf \
      I=${VCF} \
      O=${sample_name}_hg19.hard.filtered.liftover.vcf \
      CHAIN=/project/GP1/u3710062/AI_SHARE/reference/Liftover/hg38ToHg19.over.chain.gz \
      REJECT=${sample_name}_hg19.rejected.variants.vcf \
      R=/project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/hg19/ucsc.hg19.fasta
                                   
$VT_PATH decompose -s -o ${sample_name}_hg19.decomposed.vcf ${sample_name}_hg19.hard.filtered.liftover.vcf
$VT_PATH normalize -o ${sample_name}_hg19.norm.vcf -r /project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/hg19/ucsc.hg19.fasta ${sample_name}_hg19.decomposed.vcf

# run Annotation
para="ANNOVAR_${sample_name}"
ANNOVAR="/pkg/biology/ANNOVAR/ANNOVAR_20191024/table_annovar.pl"
normalized_VCF=${sample_name}_hg19.norm.vcf
humandb="/project/GP1/u3710062/AI_SHARE/shared_scripts/ANNOVAR/humandb/"

$ANNOVAR $normalized_VCF $humandb -buildver hg19 -out ${para} -remove -protocol refGene,TaiwanBiobank-official,TaiwanBiobank993WGS,avsnp150,exac03,exac03nonpsych,clinvar_20190305,gnomad211_exome,gnomad211_genome,cosmic_coding_GRCh37_v91,cosmic_noncoding_GRCh37_v91 -operation gx,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish

# filter---exonic/splicing and nonsynonymous
head -n 1 ANNOVAR_${sample_name}.hg19_multianno.txt > ${sample_name}.filtered.annotation.txt
grep -e "exonic" -e "splicing" ANNOVAR_${sample_name}.hg19_multianno.txt >> ${sample_name}.filtered.annotation.txt

