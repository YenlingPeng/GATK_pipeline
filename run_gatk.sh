#!/bin/bash

GATK_PATH=/pkg/biology/GATK/GATK_v4.1.8.0
PICARD_PATH=/pkg/biology/Picard/Picard_v2.18.11/picard.jar
BWA_PATH=/pkg/biology/BWA/BWA_v0.7.17
SAMTOOLS_PATH=/pkg/biology/SAMtools/SAMtools_v1.10/bin

sample_name=$1
FASTQ_1_PATH=$2
FASTQ_2_PATH=$3
FASTQ_READGROUP=$4
OUTPUT_PATH=$5
REF_GENOME_PATH=$6
HUMAN_DBSNP_PATH=$7
GERMLINE_RESOURCE_PATH=$8
GERMLINE_RESOURCE_FOR_PILEUP_PATH=$9
Panel_Of_Normals=$10
BAIT=$11
COVERED=$12

# Ideally, you don't need to modify following lines

NUM_THREAD=40                      

mkdir -p $OUTPUT_PATH

# Mapping
$BWA_PATH/bwa mem -t 16 -R $FASTQ_READGROUP \
        $REF_GENOME_PATH $FASTQ_1_PATH $FASTQ_2_PATH > $OUTPUT_PATH/${sample_name}.bwamem.sam 

java -Xmx48g -jar $PICARD_PATH SortSam \
         I=$OUTPUT_PATH/${sample_name}.bwamem.sam \
         O=$OUTPUT_PATH/${sample_name}.bwamem.bam \
         SORT_ORDER=coordinate \
         VALIDATION_STRINGENCY=LENIENT \
         CREATE_INDEX=true 

# ValidateSamFile
java -Xmx48g -jar $PICARD_PATH ValidateSamFile \
         I=$OUTPUT_PATH/${sample_name}.bwamem.bam \
         O=$OUTPUT_PATH/${sample_name}.validate.sam.log \
         REFERENCE_SEQUENCE=$REF_GENOME_PATH \
         MODE=SUMMARY 

# MarkDuplicates (Picard)
java -Xmx48g -jar $PICARD_PATH MarkDuplicates \
         I=$OUTPUT_PATH/${sample_name}.bwamem.bam \
         O=$OUTPUT_PATH/${sample_name}.bwamem.marked.bam \
         M=$OUTPUT_PATH/${sample_name}.metrics.txt > $OUTPUT_PATH/${sample_name}.picard.log 2>&1 \
         VALIDATION_STRINGENCY=LENIENT \
         CREATE_INDEX=true 

# Base Quality Score Recalibration (BQSR) score
### BQSR first pass
$GATK_PATH/gatk BaseRecalibrator \
         -I $OUTPUT_PATH/${sample_name}.bwamem.marked.bam \
         -R $REF_GENOME_PATH \
         --known-sites $HUMAN_DBSNP_PATH \
         -O $OUTPUT_PATH/${sample_name}.recal.pass1.table > $OUTPUT_PATH/${sample_name}.BQSR.first.pass.log 2>&1

$GATK_PATH/gatk ApplyBQSR \
         -I $OUTPUT_PATH/${sample_name}.bwamem.marked.bam \
         -R $REF_GENOME_PATH \
         --bqsr-recal-file $OUTPUT_PATH/${sample_name}.recal.pass1.table \
         -O $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.bam > $OUTPUT_PATH/${sample_name}.apply.BQSR.log 2>&1

### BQSR second pass
$GATK_PATH/gatk BaseRecalibrator \
        -I $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.bam \
        -R $REF_GENOME_PATH \
        --known-sites $HUMAN_DBSNP_PATH \
        -O $OUTPUT_PATH/${sample_name}.recal.pass2.table > $OUTPUT_PATH/${sample_name}.BQSR.second.pass.log 2>&1

### Analyze covariates
$GATK_PATH/gatk AnalyzeCovariates \
         --before-report-file $OUTPUT_PATH/${sample_name}.recal.pass1.table \
         --after-report-file $OUTPUT_PATH/${sample_name}.recal.pass2.table \
         --plots-report-file $OUTPUT_PATH/${sample_name}.covariates.pdf > $OUTPUT_PATH/${sample_name}.analyze.covariates.log 2>&1

java -Xmx48g -jar $PICARD_PATH SortSam \
         I=$OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.bam \
         O=$OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.bam \
         SORT_ORDER=coordinate \
         VALIDATION_STRINGENCY=LENIENT \
         CREATE_INDEX=true 

# Picard coverage per exon
java -Xmx48g -jar $PICARD_PATH CollectHsMetrics \
         I=$OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.bam \
         O=$OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.metrics.txt \
         R=$REF_GENOME_PATH \
         BAIT_INTERVALS=${BAIT} \
         TARGET_INTERVALS=${COVERED} > $OUTPUT_PATH/${sample_name}.collecthsmetrics.log 2>&1 \
         PER_TARGET_COVERAGE=$OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.per.target.cov \
         PER_BASE_COVERAGE=$OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.per.base.cov \
         THEORETICAL_SENSITIVITY_OUTPUT=$OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.sens.metrics.txt 
         
# HaplotypeCaller & Filteration
$GATK_PATH/gatk HaplotypeCaller \
         -R $REF_GENOME_PATH \
         -I $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.bam \
         -ERC GVCF \
         --dbsnp $HUMAN_DBSNP_PATH \
         -O $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.haplotype.SnpIndel.g.vcf.gz 

$GATK_PATH/gatk GenotypeGVCFs \
         -R $REF_GENOME_PATH \
         -V $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.haplotype.SnpIndel.g.vcf.gz \
         -O $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.haplotype.SnpIndel.vcf.gz

$GATK_PATH/gatk VariantFiltration \
         -R $REF_GENOME_PATH \
         --variant $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.haplotype.SnpIndel.vcf.gz \
         -O $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.indexed.filtered.haplotype.SnpIndel.vcf.gz \
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

$GATK_PATH/gatk GetSampleName -I $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.bam -O $OUTPUT_PATH/${sample_name}.sample_name.txt

TUMOR_SAMPLE_NAME=$(cat $OUTPUT_PATH/${sample_name}.sample_name.txt)

## Call somatic short variants and generate a BAM with Mutect2
$GATK_PATH/gatk Mutect2 \
        -R $REF_GENOME_PATH \
        -I $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.bam \
        -tumor $TUMOR_SAMPLE_NAME \
        --germline-resource $GERMLINE_RESOURCE_PATH \
        --f1r2-tar-gz $OUTPUT_PATH/f1r2.tar.gz \
        -pon $Panel_Of_Normals \
        -O $OUTPUT_PATH/${sample_name}.bwamem.mutect2.vcf.gz \
        -bamout $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.mutect2.bam \
        --native-pair-hmm-threads $NUM_THREAD > $OUTPUT_PATH/mutect2.log 2>&1 
        
$GATK_PATH/gatk LearnReadOrientationModel \
        -I $OUTPUT_PATH/f1r2.tar.gz \
        -O $OUTPUT_PATH/read-orientation-model.tar.gz\

# Calculating Contamination        
$GATK_PATH/gatk GetPileupSummaries \
        -I $OUTPUT_PATH/${sample_name}.bwamem.marked.recal.pass1.bam \
        -V $GERMLINE_RESOURCE_FOR_PILEUP_PATH \
        -L $GERMLINE_RESOURCE_FOR_PILEUP_PATH \
        -O $OUTPUT_PATH/getpileupsummaries.table > $OUTPUT_PATH/get_pileup_summaries.log 2>&1

$GATK_PATH/gatk CalculateContamination \
        -I $OUTPUT_PATH/getpileupsummaries.table \
        -tumor-segmentation $OUTPUT_PATH/${sample_name}.segments.table \
        -O $OUTPUT_PATH/contamination.table > $OUTPUT_PATH/calculate.contamination.log 2>&1

# FilterMutectCalls
$GATK_PATH/gatk FilterMutectCalls \
        -R $REF_GENOME_PATH \
        -V $OUTPUT_PATH/${sample_name}.bwamem.mutect2.vcf.gz \
        --tumor-segmentation $OUTPUT_PATH/${sample_name}.segments.table \
        --contamination-table $OUTPUT_PATH/contamination.table \
        --ob-priors $OUTPUT_PATH/read-orientation-model.tar.gz \
        -O $OUTPUT_PATH/${sample_name}.bwamem.mutect2.filtered.vcf > $OUTPUT_PATH/filter.mutect.calls.log 2>&1

