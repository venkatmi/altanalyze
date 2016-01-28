# What To Do After Getting Alternative Exon Results? #
Below are some suggestions.

## Interpreting the Results ##
When AltAnalyze was running it produced a number of output files, most to the folder AltResults/AlternativeOutput in the user output directory. These include:

  1. ExpressionOutput/DATASET-YourExperiment.txt
  1. AltResults/AlternativeOutput/YourComparison-splicing-index-exon-inclusion-results.txt
  1. AltResults/AlternativeOutput/YourComparison-splicing-index-exon-inclusion-GENE-results.txt
  1. AltResults/AlternativeOutput/YourComparison-splicing-index-ft-domain-zscores.txt
  1. AltResults/AlternativeOutput/YourComparison-splicing-index-miRNA-zscores

These files are tab-delimited text files that can be opened in a spreadsheet program like Microsoft Excel, OpenOffice or Google Documents.

**File #1** reports gene expression values for each sample and group in your probeset input expression file. The values are derived from probe sets that align to regions of a gene that are common to all transcripts and thus are informative for transcription (unless all probe sets are selected – see “Select expression analysis parameters”, above) and expressed above specified background levels.  Along with the raw gene expression values, statistics for each indicated comparison (mean expression, folds, t-test p-values) will be included along with gene annotations for that array, including putative microRNA binding sites.  This file is analogous to the results file you would have with a typical, non-exon microarray experiment and is saved to the folder “ExpressionOutput”.

Results from **files #2-5** are produced from all probe sets that may suggest alternative splicing, alterative promoter regulation, or any other variation relative to the constitutive gene expression for that gene (derived from comparisons file). Each set of results correspond to a single pair-wise comparison (e.g., cancer vs. normal) and will be named with the group names you assigned (groups file).

**File #2** reports probe sets that are alternatively regulated, based on the user defined splicing-index score and p-value. For each probe set several statistics, gene annotations and functional predictions are provided. A detailed description of all of the columns in this file is provided [here](http://altanalyze.org/image/AltAnalyze_column_annotations.xls).

**File #3** is a summarization of probeset data at the gene level from file #2. In addition to this summary, Gene Ontology terms for that gene are reported.

**Files #4 and #5** report over-representation results for protein domains (or other protein features) and microRNA-binding sites, predicted to be regulated by AltAnalyze.  These files include over-representation statistics and genes associated with the different domains or features¸ predicted to be regulated.

More information about these files can be found in the [AltAnalyze ReadMe](http://www.altanalyze.org/help_main.htm) (section 2.3).

## Down-Stream Analyses ##

At this point you have many options but some of the most common are:

  * Visualize AltAnalyze results in DomainGraph
  * Filter and sort the data in MS-Excel to find interesting genes.
  * Look for Gene Ontology terms and over-represented pathways.
  * Load the data in a pathway analysis program to see your data on pathways.
  * Find novel gene interaction networks using Cytoscape.
  * Cluster alternative exon data in Cytoscape.

### Visualizing AltAnalyze Results in DomainGraph ###
The text file results produced by AltAnalyze can be directly used as input in the protein domain and microRNA binding site visualization program, DomainGraph. DomainGraph is a plugin for the Java program Cytoscape which can be immediately opened from AltAnalyze.

### Filter and sort the data in MS-Excel to find interesting genes ###
When looking at differentially expressed genes, AltAnalyze will have exported specific comparisons, such as cancer versus normal. If you have a dataset with at least 2 replicates in each group, a ttest p-value will also be calculated for each probeset. Selecting the menu option Data>Filter in Excel will let you search for specific criterion. Looking for a fold change >2 and p<0.05 will give you an initial list to examine in more detail. You will also want to sort by p-value or fold change by going to Data>Sort.

When looking at alternative exons, you will likely wish to identify alternative exons for validation. You may wish to prioritize based on those probesets that are predicted to alter the inclusion of protein domains/motifs or microRNA binding sites or select those that correspond to specific splicing annotations (e.g., alternative cassette-exons, intron-retention, alternative promoters). A useful way to prioritize genes is to also look at the gene level alternative exon file (file #3). Since this file has all of the exon-level data agglomerated to the gene level, you can find genes with relatively few or many alternative exons. You may also wish to sort or filter by the alternative exon score (dI) or splicing p-values (e.g., MiDAS or SI).

### Visualizing Over-Represented GO terms and Pathways ###
Once over-represented pathways have been found or before doing this analysis, you can see which genes on which pathways are alternatively regulated in the program PathVisio or GenMAPP 2.1. PathVisio is a cross-platform analysis program, while GenMAPP is restricted to Windows. Both tools are easy use and have access to a large archive of curated pathways. An input file for either PathVisio or GenMAPP is found in the directory "ExpressionOutput" with the prefix "GenMAPP-". For making pathways, PathVisio or WikiPathways is recommended, since these resources produce superior pathway content (valid interactions between genes and metabolite IDs) in the same format (gpml). PathVisio can also export pathways to the GenMAPP format. A PathVisio tutorial can be [found here](http://www.pathvisio.org/wiki/PathVisioTutorials), while a GenMAPP tutorial can be [found here](http://www.genmapp.org/tutorial_v2.html).

### Find novel gene interaction networks using Cytoscape ###
You can also use the program Cytoscape to create literature based networks and view your data on networks. Tutorials are present at Cytoscape.org.

### Cluster the data ###
If you have multiple comparisons in your dataset, you may be interested in looking for global similarities in alternative expression. One way to do this is to select the option to export normalized intensities for an alternative exon analysis, select the significant alternative exon results from one or more comparison analyses, use the normalized intensities to calculate alternative exon scores for all samples (relative to an appropriate control group) and then analyze these data in a clustering program, such as [clusterMaker](http://www.cgl.ucsf.edu/cytoscape/cluster/clusterMaker.html) in Cytoscape.