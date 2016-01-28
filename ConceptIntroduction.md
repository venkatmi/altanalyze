The sections below are meant to be an introduction to concepts surrounding the analysis of modern microarray and RNA-sequencing data from biological experiments. This is recommended for students and science professionals with minimal background in the field of genomics.

# What is AltAnalyze Useful For? #

The software described on this website is used to understand how changes in the concentration of different RNAs lead to meaningful differences in cells that effect their function and identity. Since the amount of RNA is an indicator of both the activity of genes and the amount of protein produced, technologies that measure RNA levels can be useful in understanding how cells adapt to disease, genetic modifications or to external stimuli. These technologies include microarrays ([GeneChips](CompatibleArrays.md)), RNA-sequencing (RNASeq) and protein mass-spectrometry (proteomics). While the measurement methods differ and on their own can be complex, the analysis of this data is also challenging and often requires different perspectives and analytical methods (e.g., bioinformatics).

Below we provide additional details on these methods:

### [Gene Activity Measurement Technologies](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#Gene_Activity_Measurement_Technologies) ###
  * [Proteomics](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#Proteomics)
  * [Microarrays](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#Microarrays)
  * [RNA-Sequencing](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#RNA-Sequencing)

### [Distinct Forms of Gene Regulation](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#Distinct_Forms_of_Gene_Regulation) ###
  * [Transcriptional Regulation](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#Transcriptional_Regulation)
  * [Alternative Splicing](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#Alternative_Splicing)

### [Bioinformatics](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#Bioinformatics) ###
  * [Analyzing Transcriptome Profiling Data](http://code.google.com/p/altanalyze/wiki/ConceptIntroduction#Analyzing_Transcriptome_Profiling_Data)

# Gene Activity Measurement Technologies #

### Proteomics ###

Proteomics allows researches to directly measure protein peptides with different mass-to-charge ratios. While proteomics can identify the specific protein isoforms produced, RNA-Seq and more often microarrays are used due to cost considerations and sensitivity.

### Microarrays ###

Microarrays, in their most simple form, are glass slides that have small pieces of DNA attached to them to bind to cDNA produced from your RNA sample. Since DNA readily will bind to its antisense strand with high fidelity, this technique allows researchers to quickly and cheaply assay for the amounts of specific RNA in the cell. Each cDNA made from RNA is labelled with a florescent marker that allows these slides to scanned to determine the amount of cDNA that binds a probe to estimate the amount of RNA present in the cell. Based on the design of the array, researchers can detect every gene and potentially every mRNA produced from those genes (see Alternative Splicing below). However, most arrays currently focus on one isoform for each gene and hence may miss important gene regulation data. This is why the approach taken by AltAnalyze, which aims to separate out the effects of gene transcription and gene splicing as separate pieces of information.

### RNA-Sequencing ###

Microarrays are biased in their design, since they only contain probes to known RNA regions. An alternative and unbiased method is RNA sequencing, in which the cDNA produced from the RNA in a cell (using common laboratory reagents) can be directly sequenced and quantified. We still make cDNA from the RNA (using a technique called reverse transcription) since DNA is more stable that RNA and can be more easily amplified if needed. Using RNASeq we can determine the amounts of each RNA in the cell relatively well, assuming the sequencing reads are long enough (number of base-pairs), we sequence most of the transcripts in the cell (depth of sequencing) and have good software to predict all of the exons and junctions that are present (see Alternative Splicing below).

# Distinct Forms of Gene Regulation #


### Transcriptional Regulation ###

The major way in which the cell controls control the amount of protein made from a gene is through a process called transcription. Transcription occurs when the DNA on a chromosome is accessible to binding by transcription factors (proteins) to produce RNA. Other factors that influence the amount of protein produced include regulation by non-coding RNAs, such as microRNAs and alternative splicing. Bioinformatically, we can differentiate gene transcription from other modes of gene regulation by measuring the amount of RNA present for specific gene features, notably the exons of genes which are most common to all mRNA isoforms produced for that gene (see Alternative Splicing below).

### Alternative Splicing ###

Genes are typically composed of very long stretches of DNA sequence, only a small fraction of which actually encodes for protein sequence. The regions of the gene that encode for proteins are known as exons and the interspersing non-coding sequence is referred to as introns. Introns are removed through the process of splicing in the cell nucleus, whereby exon sequences are joined together thereby removing the intron sequences. The resulting small sequence is exported from the nucleus and translated within ribosomes to proteins. Alternative splicing is the process in which not only introns but exons are retained or excluded from the processed mRNAs in varying amounts. Which exons are retained is determined by factors (proteins and RNAs) that bind to the pre-processed RNA and influence exon inclusion. These factors include kinases which are activated or inhibited by cell signaling changes that occur due to the activity of receptors on the surface of the cell. Often, different cell types produce different combinations of exons in mRNAs that result in different protein sequences and often with different protein function.

# Bioinformatics #

Bioinformatics is a broad field that applies to various biological disciplines and computational approaches. Bioinformaticians can have broad expertise in multiple disciplines or can be highly specialized experts in one or more fields. Specific examples can be found [here](http://theconversation.edu.au/explainer-what-is-bioinformatics-9911#).

### Analyzing Transcriptome Profiling Data ###

When it comes to RNA expression data, also known as transcriptome profiling, there are a large variety of different analytical approaches. Most of these center on the statistical analysis of gene data with the aim of determining significant differences in two biological sample groups (e.g., treatment and control). Some statistical analysis environments, such as [R](http://cran.org), allow users to perform nearly all transcriptome related bioinformatics analyses, from processing raw microarray image files to visualizing RNA differences as colors on biological pathways. Most researchers using R has a background in statistics or computer programming, hence, the associated tools can be restrictive for non-computationally oriented biologists.

AltAnalyze was designed to be a powerful analytical resource for computational biologists and non-computational biologists alike. The software is available as both a simple to use graphical interface as well as an optional command-line version for remote use. After clicking through a few simple screens, the user will receive a highly structured set of results which often meet the bioinformatics needs of the researcher. In addition to implementing many of the most common transcriptome analyses, such as microarray data [normalization](Tutorials.md), [quality control](QualityControl.md), [RNASeq](RNASeqExpressionNorm.md) gene expression analyses, [clustering/heatmap visualization](Heatmaps.md), moderated [statistical tests](AltAnalyzeStatistics.md), Gene Ontology [over-representation](PathwayAnalysis.md) and pathway [visualization](http://code.google.com/p/go-elite/wiki/Tutorial_GUI_version#Data_Visualization), AltAnalyze performs a number of incredibly sophisticated analyses not available in any other tool. These include the: (1) determination of alternative splicing from underlying [exon and junction](AltAnalyzeStatistics.md) expression changes, (2) a [functional assessment](DomainAnalysis.md) of these differences in the context of known protein isoforms and functional protein/RNA sequences, (3) analysis of likely cell types present in a given sample (LineageProfiler) and (4) analysis of possible [transcription factor](http://code.google.com/p/go-elite/wiki/GeneSetRepository#Transcription_Factor_to_Target_Genes) regulation. As a result, AltAnalyze is the most streamlined and analytically comprehensive single analysis workflow that we are aware of.