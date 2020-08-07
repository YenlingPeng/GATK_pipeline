# GATK_pipeline

- GATK 4.1.8 version
- Call somatic short variants and generate a bamout with Mutect2
### Processing Steps

#### 1. Lift genome positions
UCSC provides binary liftOver tools to convert BED file from one genome assembly to another.
> NOTE: Use the 'chr' before each chromosome name

```
liftoverbed.sh
```

#### 2. BedToIntervalList (Picard)
Converts from BED to the Picard interval_list format.
**This tool requires a sequence dictionary, provided with the SEQUENCE_DICTIONARY or SD argument.**

```
run_BedTolntervalList.sh
```

### Running pipeline

```
run_gatk.sh
```

***

### Reference Links
- [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)

- [LiftOver](https://genome.sph.umich.edu/wiki/LiftOver#Binary_liftOver_tool0)

- [BedToIntervalList (Picard)](https://gatk.broadinstitute.org/hc/en-us/articles/360036883931-BedToIntervalList-Picard-)
