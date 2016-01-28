# Instructions for Obtaining RNASeq Input Files #

For RNA-Seq analyses, AltAnalyze accepts either unaligned RNA-sequence read files (fastq), genome/transcriptome alignments (BAM), extracted reads per exon and junction (BED, BioScope, TCGA) and already normalized gene expression values (RPKM, FPKM, CPM, TPM). For deep splicing analyses, junction-level results are required which are only provided from external alignment programs (e.g., TopHat2, STAR).

### Downloading Sample FASTQ Data ###

You are welcome to try with your own FASTQ data, as this is typically easiest for most users rather than downloading large datasets. For optimal junction analysis results, we recommend >50nt reads, however, we have achieved acceptable results from 35nt reads. If you currently do not have your own data, example files and download instructions can be found [here](DownloadFASTQ.md).

### Direct Sequence Alignment in AltAnalyze ###

Included in AltAnalyze version 2.0.9.1, is the software [kallisto](http://pachterlab.github.io/kallisto/), that is capable of aligning an RNA-Seq FASTQ or paired-end FASTQ file in under 3 minutes, on a standard laptop or desktop. Please read the Kallisto license before using, as some restrictions apply. This analysis will produce Ensembl isoform TPM (transcript per million) estimates from which AltAnalyze will calculate gene-level TPM values for all Ensembl genes associated with the Ensembl database version used by AltAnalyze. As such, analysis is compatible with any AltAnalyze supported species (species support varies with different database versions). Please check which database versions from AltAnalyze are most appropriate or contact us for assistance.

### Using TopHat and Bowtie ###

Currently, we recommend two pipelines for producing input data for AltAnalyze: 1) TopHat and 2) BioScope. TopHat should work with most all RNA-seq datasets, while BioScope is compatible with ABI SOLiD data. As TopHat is open-source and runs on typical hardware configurations, this solution will be more appropriate for most users.

For TopHat, junction BED files can be easily obtained using simple command-line arguments to TopHat, once Bowtie and indexed genome files have been downloaded. Junction.bed files produced by TopHat can be directly used in AltAnalyze or the aligned BAM format files can now be directly analyzed. If BAM processing compatibility issues exist, other options for obtaining exon and junction expression files (BED format) can be found [here](BAMtoBED.md).

### Using BioScope ###

Details on BioScope usage can be found [here](http://www3.appliedbiosystems.com/cms/groups/global_marketing_group/documents/generaldocuments/cms_074971.pdf). Once result files are produced, AltAnalyze will accept the junction and exon reads files produced from this pipeline as input. Prior to import, the extension of these two files must be changed to .tab in order for AltAnalyze to recognize these from other text files. To ensure the exon and junction files for each sample are properly matched, make sure that the exon and junction file have the same name, followed by a double underscore, followed by a distinguishing annotation. For example:

|`Cancer_s1__canonical-junction.tab`|
|:----------------------------------|
|`Cancer_s1__noncanonical-junction.tab`|
|`Cancer_s1__exon.tab`              |
|`Wt_s1__canonical-junction.tab`    |
|`Wt_s1__noncanonical-junction.tab` |
|`Wt_s1__exon.tab`                  |