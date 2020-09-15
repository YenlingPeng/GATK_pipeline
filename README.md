# GATK_pipeline

- GATK 4.1.8 version
- Call somatic short variants and generate a bamout with Mutect2
### Processing Steps

#### 1. Lift genome positions
UCSC provides binary liftOver tools to convert BED file from one genome assembly to another.
> NOTE: Use the 'chr' before each chromosome name

#### 2. BedToIntervalList (Picard)
Converts from BED to the Picard interval_list format.
**This tool requires a sequence dictionary, provided with the SEQUENCE_DICTIONARY or SD argument.**

```liftoverbed_IntervalList.sh```

### Running pipeline

The script for pipeline ```run_gatk.sh```
Submit a Job: An example is available [here](https://github.com/YenlingPeng/GATK_pipeline/blob/master/gatk_submit.sh)

### Annotation 
liftover VCF from hg38 to hg19 & ANNOVAR

The script for ANNOVAR ```run_annotation.sh```
Submit a Job: An example is available [here](https://github.com/YenlingPeng/GATK_pipeline/blob/master/annotation_submit.sh)

***

### Reference Links
- [LiftOver](https://genome.sph.umich.edu/wiki/LiftOver#Binary_liftOver_tool0)

- [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)

- [BedToIntervalList (Picard)](https://gatk.broadinstitute.org/hc/en-us/articles/360036883931-BedToIntervalList-Picard-)
