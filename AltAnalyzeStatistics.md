## Statistical Methods Implemented in AltAnalyze ##

AltAnalyze implements various previously described algorithms for gene expression quantification, differential gene expression, alternative exon expression and feature/pathway over-representation. These algorithms are implemented in the context of a rigorous exon/transcript annotation framework, based on data from Ensembl and UCSC mRNAs. Details on all of these methods is **[provided here](http://www.altanalyze.org/help.htm)**.

In general, implemented methods include:

### Gene Expression Analyses ###

  * Robust Multichip Analysis (**Affymetrix**)
  * Detection Above Background p-values (**Affymetrix**)
  * Quantile normalization (**RNASeq**)
  * RPKM (gene, exon, junction) (**RNASeq**)
  * Geometric fold calculation
  * [Moderated and conventional two-sample tests and f-test](ModeratedTest.md)
  * Benjamini-Hochberg adjusted FDR test

### [Alternative Exon Analyses](ScoringMethod.md) ###

  * Reciprocal-junction analyses (ASPIRE, Linear regression)
  * Alternative exon (splicing-index, FIRMA, MiDAS, adjusted f-tests)

### [Over-representation](http://genmapp.org/go_elite/help.htm#s4) ###

  * Fisher Exact Test
  * Z score calculation
  * Permutation-based p-values
  * Benjamini-Hochberg adjusted FDR test

### LineageProfiler analysis ###

  * Pearson/Spearman's correlation
  * Z score calculation