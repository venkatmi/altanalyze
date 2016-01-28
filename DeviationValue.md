# Deviation Value Splicing-Index Statistic #

A splicing-index statistic is a common algorithm used for detecting alternative exon expression. While this statistic provides one useful method for detecting alternative exons, it can be used in conjunction with other statistics, filters and other parameters to increase specificity of predictions prior to validation.

One of these additional statistical measures is called a [deviation value](http://www.nature.com/jhg/journal/v57/n6/full/jhg201237a.html#bib23). In AltAnalyze, a deviation value is also included from the splicing-index statistic (Yamashita et. al, 2012 Journal of Human Genetics) across the entire gene to evaluate overall gene-level variation in exon-expression relative to the exon-region of interest. The authors of the study that introduce the deviation value or DV, recommend a DV>3.0 be considered as one variable contributing to a greater likelihood of alternative exon differential expression.

### Additional Analysis Considerations ###

While the deviation value or DV is a useful statistic to consider, we recommend also running AltAnalyze with the [FIRMA](FIRMA.md) algorithm where available (Affymetrix data only), as this statistic is likely a more rigorous measure of overall significance at the gene-level for a single regulated exon.

In addition to DV, additional parameters recommended by Yamashita et. al include:

  1. Remove genes with less than four exons.
  1. DABG p<0.01 in one or both groups
  1. At least 15% of all analyzed probesets for a gene should meet #2 (minimum of 3 per gene).
  1. Average expression should be greater than 150
  1. Utilize conservative gene annotation predictions

In AltAnalyze these are addressed in the following way:
  1. After getting alternative exon results, in the alternative-exon-results.txt file, find the column "distal exon-region-ID", find and replace "E" in this field with nothing, filter out any line where the exon region number is less than 4.
  1. Expression Analysis Parameters: Remove probesets with a DABG p-value above (default is 0.05 but you can set to 0.01). Command-line option: --dabgp 0.01
  1. This is not exactly replicated in AltAnalyze, but for any constitutive probesets to make it into the gene expression analysis, they must be match the defaults specified in #2 and #4. However, if none exist, then any predicted constitutive probesets are used. This is actually not too useful for exon array analysis, since gene expression is recalculated in the splicing set of functions, but is important for RNA-Seq where this exact value is reused in the splicing functions. I will consider this for splicing array analyses in the future as a new option.
  1. Expression Analysis Parameters: Remove probesets expressed below (default is 1 and use to be 70). Command-line option: --rawexp 100
  1. AltAnalyze does not use Affymetrix transcript clusters as a gene model but rather derives its own for each build of Ensembl using the Ensembl exon gene structure as a scaffold and augmenting with UCSC mRNAs that share at least one exon.