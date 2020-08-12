#!/bin/bash
#PBS -q ntu192G
#PBS -l select=1:ncpus=5
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N gatk4.1_mutect2
#PBS -o /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/gatk.out
#PBS -e /work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/gatk.err
#PBS -M s0890003@gmail.com
#PBS -m e

BAIT="/work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/Panel_targets_EGFR_v1.interval_list"
COVERED="/work2/lynn88065/GATK/Inputs/20200721_GATK_pipeline/Panel_targets_EGFR_v1.interval_list"

# Mutect2
GATK_PATH=/pkg/biology/GATK/GATK_v4.1.8.0
PICARD_PATH=/pkg/biology/Picard/Picard_v2.18.11/picard.jar
BWA_PATH=/pkg/biology/BWA/BWA_v0.7.17
SAMTOOLS_PATH=/pkg/biology/SAMtools/SAMtools_v1.10/bin

TUMOR_FASTQ_1_PATH="/project/GP1/u3710062/AI_SHARE/Temp/20191018_DataToYL/NA12878-D702-D508_S38_L001_R1_001.fastq.gz"
TUMOR_FASTQ_2_PATH="/project/GP1/u3710062/AI_SHARE/Temp/20191018_DataToYL/NA12878-D702-D508_S38_L001_R2_001.fastq.gz"
TUMOR_READGROUP='@RG\tID:NA12878-D702-D508_S38\tLB:NA12878-D702-D508_S38\tSM:NA12878-D702-D508_S38\tPL:ILLUMINA\'
OUTPUT_PATH="/work2/lynn88065/GATK/Outputs/20200721_GATK_pipeline"
REF_GENOME_PATH="/project/GP1/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta"
HUMAN_DBSNP_PATH="/project/GP1/reference/Homo_sapiens/GATK/hg38/dbsnp_146.hg38.vcf.gz"
GERMLINE_RESOURCE_PATH="/project/GP1/reference/Homo_sapiens/GATK/Mutect2/af-only-gnomad.hg38.vcf.gz"
GERMLINE_RESOURCE_FOR_PILEUP_PATH="/project/GP1/reference/Homo_sapiens/GATK/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz"
NUM_THREAD=40                      
Panel_Of_Normals="/project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/2.8/hg38/somatic-hg38-1000g_pon.hg38.vcf.gz"

bwa mem -t 16 -R $TUMOR_READGROUP \
        $REF_GENOME_PATH $TUMOR_FASTQ_1_PATH $TUMOR_FASTQ_2_PATH > $OUTPUT_PATH/tumor.sam &&

java -Xmx48g -jar $PICARD_PATH SortSam \
         I=$OUTPUT_PATH/tumor.sam \
         O=$OUTPUT_PATH/tumor.bam \
         SORT_ORDER=coordinate \
         VALIDATION_STRINGENCY=LENIENT \
         CREATE_INDEX=true &&

# ValidateSamFile
java -Xmx48g -jar $PICARD_PATH ValidateSamFile \
         I=$OUTPUT_PATH/tumor.bam \
         O=$OUTPUT_PATH/tumor_validate_sam.log \
         REFERENCE_SEQUENCE=$REF_GENOME_PATH \
         MODE=SUMMARY &&

# MarkDuplicates (Picard)
java -Xmx48g -jar $PICARD_PATH MarkDuplicates \
         I=$OUTPUT_PATH/tumor.bam \
         O=$OUTPUT_PATH/tumor_marked.bam \
         M=$OUTPUT_PATH/tumor_metrics.txt > $OUTPUT_PATH/picard_tumor.log 2>&1 \
         VALIDATION_STRINGENCY=LENIENT \
         CREATE_INDEX=true &&

# Base Quality Score Recalibration (BQSR) score
### BQSR first pass
$GATK_PATH/gatk BaseRecalibrator \
         -I $OUTPUT_PATH/tumor_marked.bam \
         -R $REF_GENOME_PATH \
         --known-sites $HUMAN_DBSNP_PATH \
         -O $OUTPUT_PATH/tumor_recal_pass1.table > $OUTPUT_PATH/tumor_BQSR_first_pass.log 2>&1

$GATK_PATH/gatk ApplyBQSR \
         -I $OUTPUT_PATH/tumor_marked.bam \
         -R $REF_GENOME_PATH \
         --bqsr-recal-file $OUTPUT_PATH/tumor_recal_pass1.table \
         -O $OUTPUT_PATH/tumor_marked.recal.pass1.bam > $OUTPUT_PATH/tumor_apply_BQSR.log 2>&1

### BQSR second pass
$GATK_PATH/gatk BaseRecalibrator \
        -I $OUTPUT_PATH/tumor_marked.recal.pass1.bam \
        -R $REF_GENOME_PATH \
        --known-sites $HUMAN_DBSNP_PATH \
        -O $OUTPUT_PATH/tumor_recal_pass2.table > $OUTPUT_PATH/tumor_BQSR_second_pass.log 2>&1

### Analyze covariates
$GATK_PATH/gatk AnalyzeCovariates \
         --before-report-file $OUTPUT_PATH/tumor_recal_pass1.table \
         --after-report-file $OUTPUT_PATH/tumor_recal_pass2.table \
         --plots-report-file $OUTPUT_PATH/tumor_covariates.pdf > $OUTPUT_PATH/tumor_analyze_covariates.log 2>&1

java -Xmx48g -jar $PICARD_PATH SortSam \
         I=$OUTPUT_PATH/tumor_marked.recal.pass1.bam \
         O=$OUTPUT_PATH/tumor_marked.recal.pass1.indexed.bam \
         SORT_ORDER=coordinate \
         VALIDATION_STRINGENCY=LENIENT \
         CREATE_INDEX=true &&

# Picard coverage per exon
java -Xmx48g -jar $PICARD_PATH CollectHsMetrics \
         I=$OUTPUT_PATH/tumor_marked.recal.pass1.indexed.bam \
         O=$OUTPUT_PATH/tumor_marked.recal.pass1.indexed.metrics.txt \
         R=$REF_GENOME_PATH \
         BAIT_INTERVALS=${BAIT} \
         TARGET_INTERVALS=${COVERED} > $OUTPUT_PATH/tumor_collecthsmetrics.log 2>&1 \
         PER_TARGET_COVERAGE=$OUTPUT_PATH/tumor_marked.recal.pass1.indexed.per.target.cov \
         PER_BASE_COVERAGE=$OUTPUT_PATH/tumor_marked.recal.pass1.indexed.per.base.cov \
         THEORETICAL_SENSITIVITY_OUTPUT=$OUTPUT_PATH/tumor_marked.recal.pass1.indexed.sens.metrics.txt &&
         
# HaplotypeCaller & Filteration
$GATK_PATH/gatk HaplotypeCaller \
         -R $REF_GENOME_PATH \
         -I $OUTPUT_PATH/tumor_marked.recal.pass1.indexed.bam \
         -ERC GVCF \
         --dbsnp $HUMAN_DBSNP_PATH \
         -O $OUTPUT_PATH/tumor_marked.recal.pass1.indexed.haplotype.SnpIndel.g.vcf.gz &&

$GATK_PATH/gatk GenotypeGVCFs \
         -R $REF_GENOME_PATH \
         -V $OUTPUT_PATH/tumor_marked.recal.pass1.indexed.haplotype.SnpIndel.g.vcf.gz \
         -O $OUTPUT_PATH/tumor_marked.recal.pass1.indexed.haplotype.SnpIndel.vcf.gz &&

$GATK_PATH/gatk VariantFiltration \
         -R $REF_GENOME_PATH \
         --variant $OUTPUT_PATH/tumor_marked.recal.pass1.indexed.haplotype.SnpIndel.vcf.gz \
         -O $OUTPUT_PATH/tumor_marked.recal.pass1.indexed.filtered.haplotype.SnpIndel.vcf.gz \
         --cluster-window-size 10 \
         --filter-expression "DP < 5" \
         --filter-name "LowCoverage" \
         --filter-expression "QUAL < 30.0" \
         --filter-name "VeryLowQual" \
         --filter-expression "QUAL > 30.0 && QUAL < 50.0" \
         --filter-name "LowQual" \
         --filter-expression "QD < 1.5" \
         --filter-name "LowQD"

# Mutect2
## get sample names that will be used for later in several command line calls 

$GATK_PATH/gatk GetSampleName -I $OUTPUT_PATH/tumor_marked.recal.pass1.bam -O $OUTPUT_PATH/tumor_sample_name.txt
#$GATK_PATH/gatk GetSampleName -I $OUTPUT_PATH/normal_marked.recal.pass1.bam -O $OUTPUT_PATH/normal_sample_name.txt

TUMOR_SAMPLE_NAME=$(cat $OUTPUT_PATH/tumor_sample_name.txt)
#NORMAL_SAMPLE_NAME=$(cat $OUTPUT_PATH/normal_sample_name.txt)

## Call somatic short variants and generate a BAM with Mutect2
$GATK_PATH/gatk Mutect2 \
        -R $REF_GENOME_PATH \
        -I $OUTPUT_PATH/tumor_marked.recal.pass1.bam \
        -tumor $TUMOR_SAMPLE_NAME \
        --germline-resource $GERMLINE_RESOURCE_PATH \
        --f1r2-tar-gz $OUTPUT_PATH/f1r2.tar.gz \
        -pon $Panel_Of_Normals \
        -O $OUTPUT_PATH/Mutect2.vcf.gz \
        -L $BAIT 
        -bamout $OUTPUT_PATH/Mutect2.bam \
        --native-pair-hmm-threads $NUM_THREAD > $OUTPUT_PATH/mutect2.log 2>&1 
        

$GATK_PATH/gatk LearnReadOrientationModel \
        -I $OUTPUT_PATH/f1r2.tar.gz \
        -O $OUTPUT_PATH/read-orientation-model.tar.gz\

# Calculating Contamination        
$GATK_PATH/gatk GetPileupSummaries \
        -I $OUTPUT_PATH/tumor_marked.recal.pass1.bam \
        -V $GERMLINE_RESOURCE_FOR_PILEUP_PATH \
        -L $GERMLINE_RESOURCE_FOR_PILEUP_PATH \
        -O $OUTPUT_PATH/getpileupsummaries.table > $OUTPUT_PATH/get_pileup_summaries.log 2>&1

$GATK_PATH/gatk CalculateContamination \
        -I $OUTPUT_PATH/getpileupsummaries.table \
        -tumor-segmentation $OUTPUT_PATH/tumor_segments.table \
        -O $OUTPUT_PATH/contamination.table > $OUTPUT_PATH/calculate_contamination.log 2>&1

# FilterMutectCalls
$GATK_PATH/gatk FilterMutectCalls \
        -R $REF_GENOME_PATH \
        -V $OUTPUT_PATH/Mutect2.vcf.gz \
        --tumor-segmentation $OUTPUT_PATH/tumor_segments.table \
        --contamination-table $OUTPUT_PATH/contamination.table \
        --ob-priors $OUTPUT_PATH/read-orientation-model.tar.gz \
        -O $OUTPUT_PATH/filtered.vcf > $OUTPUT_PATH/filter_mutect_calls.log 2>&1

