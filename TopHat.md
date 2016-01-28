## Introduction ##

TopHat is a component of the Tuxedo suite of tools developed for the analysis of short sequence reads from various high-throughput sequencing platforms. TopHat works in conjunction with the short read mapper Bowtie to align sequence reads to the genome, identify and quantify junction expression (known and novel).

## Details ##

Current version details on TopHat can be found [here](http://tophat.cbcb.umd.edu/).

### Usage with AltAnalyze ###

AltAnalyze can import junction alignment data from TopHat junction BED files and exon alignment data from the associated accepted\_hit.bam file. See TopHat documentation for the latest details on usage with the sequencing technology you using. Below is an example for Illumina single read data.

**Analyze a single non-paired end sample FASTQ file**

_Illumina single read_
```
tophat --GTF /work/Common/Data/Reference_RNASeq_Assembly_GTF/Homo_sapiens.GRCh37.56.chr.gtf --keep-tmp --output-dir /home/user/tophat_data/NP902 hg19 /home/user/hESC-NP/SRS011902-NP.fq
```

_Additional options: Paired-end data_

`--mate-inner-dist 104 --mate-std-dev 17 --solexa1.3-quals --segment-length 18`

_Additional options: Exclude novel junctions_

`--no-novel-juncs`

**Obtaining exon expression estimates from TopHat BAM files**

Install [BEDTools](BEDTools.md) and call the utility bamToBed (recognized on Unix systems once BEDTools has been added to the local or global .bashrc file). The file accepted\_hits.bam is produced with each TopHat run in the same output directory as the junction BED file.

In the below example, "hESC\_differentiation\_exons.bed" is produced by AltAnalyze prior to running BEDTools ([see instructions here](BAMtoBED.md)), containing all known mRNA exon region coordinates from Ensembl/UCSC and all novel exon coordinates indicated from the TopHat junction BED results. These methods should work equivalently for non-TopHat produced BAM files, however, additional sorting of the BAM file may be required (e.g., SAMTools).

_Build Exon BED file from BAM_
```
bamToBed -i accepted_hits.bam -split| coverageBed -a stdin -b /home/user/BAMtoBED/hESC_differentiation_exons.bed > /home/user/RNASeqStudy/Sample1/day0_s1__exons.bed
```

**File naming**

If one or multiple junction and exon alignment files are produced, specific naming conventions are needed to ensure that different files for the the same biological samples are properly matched up. Below is an example of multiple junction and exon files in a single directory for import into AltAnalyze.

|`Cancer_s1__canonical-junction.bed`|
|:----------------------------------|
|`Cancer_s1__noncanonical-junction.bed`|
|`Cancer_s1__exon.bed`              |
|`Wt_s1__canonical-junction.bed`    |
|`Wt_s1__noncanonical-junction.bed` |
|`Wt_s1__exon.bed`                  |

When imported into AltAnalyze, any file name proceeding a double underscore will be recognized as the same sample and merged with the other files. Exon and junction BED files are differentiated based on the number of columns present in the two files (12 for junction, 10 for exon).