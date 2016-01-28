# Introduction #

This is the folder which contains all AltAnalyze relationship files used to:
  * perform microarray normalization (library files, [APT](APT.md))
  * Annotate gene expression results
  * Identify which probesets align to which exons and annotations
  * Associate proteins, domains and motifs with probesets
  * Determine reciprocal junctions to compare for junction arrays

Nearly all AltAnalyze annotations are pre-built as opposed to extracted directly from array manufacturer files (Affymetrix CSV files). Files are downloaded from the main menu by selecting a species to analyze.

# Details #

Most users will not need to ever open this directory unless they need to do the following:
  * Look up array annotations in batch (e.g., AltDatabase/EnsMart55/Hs/exon/Hs\_Ensembl\_probesets.text).
  * Customize database files (developers and advanced users only)
  * Copy Affymetrix library and annotation files for other purposes

Multiple database versions will be stored in this directory for different versions of Ensembl (e.g., Ensembl version 55 - EnsMart55).