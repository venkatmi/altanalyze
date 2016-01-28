![http://altanalyze.org/graphics/pcbc.png](http://altanalyze.org/graphics/pcbc.png)

## Overview of PCBC C4 Compendium Genomic Profiles in Synapse ##

Information on the NHLBI Progenitor Cell Biology Consortium (PCBC) can be found at:
http://www.progenitorcells.org

AltAnalyze was used to generate signature sets (differentially expressed genes) for a large number of pluripotency and iPS associated comparisons. These comparisons include human embryonic stem cells, derived differentiation products and control profiled cell types. The same data exists for a large number of induced pluripotent stem (iPS) cell lines, profiled at the level of mRNA-seq, miRNA-seq, copy-number variation and methylation based assays. The associated raw data will be posted for public access in the [Sage Bionetworks](http://sagebase.org/synapse-overview/) open data repository Synapse, following internal review.

This page describes the criterion used in AltAnalyze to generate automated comparisons lists, processed data, alternative exon results, lineage predictions, enrichment analysis results and associated graphical outputs.

### Version 1.0 Compendium Results ###

The first version of comprehensive comparison analyses posted in Synapse were generated on 9/22/13 using AltAnalyze version 2.0.8 using the human EnsMart72 database.

**Parameters**: All parameters were set to their defaults except for the following:
  1. --geneRPKM 3 (gene-level RPKM > 3)
  1. --GEcutoff 3 (maximum gene expression fold change of 3 for alternative exon analysis)
  1. --GEeliteptype adjp (BH adjusted, moderated t-test p-value for all comparisons)

These defaults include a 2-fold gene expression cutoff.