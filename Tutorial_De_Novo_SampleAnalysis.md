# Tutorial 4 - De Novo Identification of Gene-Sample Clusters #

A challenging bioinformatics challenge is the unbiased identification of both major and sub-populations of bulk or single-cell transcriptomes where no obvious reference gene or splice signatures exist. In addition, while transcriptome samples are often annotated as belonging to one or more groups, miss-classified samples or other unexpected covariants may also result in atypical variation that is difficult to characterize. A new set of algorithms have been implemented in AltAnalyze to identify the most coherent gene and sample expression signatures, prior to group comparison analyses. These new methods become particularly important when considering unknown single-cell RNA-Seq profiles and tumor subtype classification. In addition to gene expression, the software allows for the highly accurate and unbiased identification of highly similar alternative splicing signatures.

## Algorithms ##

There are currently two major methods implemented in the de novo sample discovery workflow:
  * [Unsupervised ICGS](ICGS.md)
  * Supervised ICGS

Following both ICGS analyses, the user will be presented with the option to select one of the clustered outputs for the standard pairwise comparison workflow (all possible sample cluster in comparisons). Both workflows can optionally exclude genes with low overall expression, low differential expression and ncRNAs. The expression thresholds are used defined. Differential expression is carried out between a set of the lowest and highest expressed measurements for a gene across samples (user defined cut point). Both allow for the additional removal of cell-cycle associated gene-signatures using either a stringent or method or conservative method.

### Unsupervised ICGS ###
Details on the unsupervised ICGS workflow can be found [here](ICGS.md). This analysis performs iterative rounds of filtering and clustering to obtain novel supervising genes to pull in the highly correlated genes sets for clustering discovery. Different thresholds will likely be required (see below).

### Supervised ICGS ###
The supervised workflow accepts a single gene, multiple genes, or gene sets available from GO-Elite (e.g., disease ontologies, gene ontology, pathways, cell-type specific genes) or prior signatures stored by the user (Additional Analyses > heatmap or PCA analyses). The workflow is identical to the unsupervised workflow, except following the Pioneering Round ([Step 1](ICGS.md)), the initially filtered dataset (RPKM, gene counts and fold change) is further filtered for any genes in the supplied genes or gene gets, where at least 4 genes are intra-correlated those genes. When a custom gene list is supplied, the software finds any genes correlated to those in the expression filtered dataset. The resulting genes are used as input for all remaining ICGS steps.

## Input Data Formats ##

These analyses are compatible with all regular supported file formats in AltAnalyze (e.g., RNA-Seq bed, Affymetrix CEL, Agilent feature extraction, FPKM/TPM/RPKM txt files). For RNA-Seq we recommend using either junction.bed files directly (TopHat, SpliceMap), TCGA junction\_expression.txt files and/or BedTools produced exon.bed files. Already normalized FPKM/TPM/RPKM files, can be loaded as the "vendor" "Other ID" and "platform" corresponding to the associated gene identifier (e.g., Ensembl, Symbol). For more information on these inputs, go http://code.google.com/p/altanalyze/wiki/Tutorial_AltExpression_RNASeq here].

## Running ICGS in the Graphical Interface ##

The instructions for running the de novo analysis are identical to all other workflow analyses in AltAnalyze (see [tutorials](https://code.google.com/p/altanalyze/wiki/Tutorials)), however, an additional step allows the user to discover groups by clustering.

  * Open the AltAnalyze program folder and open the file “AltAnalyze". In Windows, this file has the extension “.exe”. If you are working on a Linux machine or are having problems starting AltAnalyze, you can also start the program directly from the [code](source.md).
  * Follow the steps in Tutorials 1, 2 or 3, until you arrive at the [AltAnalyze: Assign CEL files to a Group](http://altanalyze.org/image/AssignCELGroups.jpg) interface. Select the option at the top of the menu named: _Run de novo cluster prediction (ICGS) to discover groups_.
  * Accept the default parameters listed or modify to customize. For single-cell analyses, we recommend an RPKM >= 1, fold > 4, Pearson correlation cutoff > 0.4 and at least 3 samples differing. For column clustering, the HOPACH algorithm (requires that R is installed, see below) is recommended but not required. For supervised analyses, input space delimited (or Excel pasted lists) of standard accepted gene symbols for the species analyzed in the field "(optional) Enter genes to build clusters from". To display these genes in the larger cluster heatmap, enter these or other genes into the field "(optional) Display selected gene IDs in results". Alternatively, you can select a GeneSet to perform your supervised analysis upon.
  * If you wish to remove genes clusters that are driven by cell-cycle dependent changes, select the option: "(optional) Eliminate cell cycle effects". Please note, if the changes in your samples in your set are primarily distinct due to a major cell cycle effect (e.g., cancer), you may not want to add this option. However, you can try this and then simply re-run when prompted to see both possible outcomes.
  * Once completed, you will be presented with a series of small heatmap thumbnails. Cell and tissue-specific BioMarkers enriched in each gene cluster will be shown to the left of each heatmap when selected. Red genes to the right of the heatmap indicate the software selected guide genes.
  * Selecting a cluster number at the top of the interface will direct AltAnalyze to use the associated sample cluster groups for performing all group pairwise comparison analyses, using the options the user established in earlier menus.
  * Re-run the analysis when prompted using different gene sets,  different filtering parameters or different clustering options (e.g., HOPACH) to select the most informative sample groupings.

_**Note**_: For HOPACH clustering analyses, you must have the R programming environment installed. On Windows, R must also be registered in your system path. Once R is registered by your system, AltAnalyze can install any necessary packages locally in the AltAnalyze Config directory automatically.

## Running ICGS by Command-Line ##

**Unsupervised ICGS from FPKM/TPM/RPKM tab-delimited text file**
```
python AltAnalyze.py --runICGS yes --expdir "/Users/SingleCell/Myoblast/myoblast.txt" --platform RNASeq --species Hs --column_method hopach --rho 0.4 --ExpressionCutoff 4 --FoldDiff 4 --SamplesDiffering 3 --excludeCellCycle conservative
```


**Supervised ICGS from FPKM/TPM/RPKM tab-delimited text file**
```
python AltAnalyze.py --runICGS yes --expdir "/Users/SingleCell/Myoblast/myoblast.txt" --platform RNASeq --species Hs --GeneSetSelection BioMarkers --PathwaySelection Myotube --PathwaySelection "Mesenchymal Stem Cells" --column_method hopach --rho 0.4 --ExpressionCutoff 4 --justShowTheseIDs "NKX2-5 T TBX5 NANOG SOX9" --FoldDiff 4 --SamplesDiffering 3 --excludeCellCycle conservative
```