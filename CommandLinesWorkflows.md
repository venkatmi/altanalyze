# Workflow Command Line Options in AltAnalyze #

Full analysis pipelines can be run in AltAnalyze from the command line (source code and binary AltAnalyze versions). Below are example command line examples and detailed descriptions of available flags. For more information and links see: http://code.google.com/p/altanalyze/wiki/CommandLineMode

### Workflow Examples ###

_**Analyzing CEL files – Affymetrix 3’ array using default options and GO-Elite**_
```
python AltAnalyze.py --species Mm --arraytype "3'array" --celdir "C:/CELFiles" --groupdir "C:/CELFiles/groups.CancerCompendium.txt" --compdir "C:/CELFiles/comps.CancerCompendium.txt" --output "C:/CELFiles" --expname "CancerCompendium" --runGOElite yes --returnPathways all
```

_**Analyzing RNA-Seq (RNASeq) data – BED files using default options**_
```
python AltAnalyze.py --species Mm --platform RNASeq --bedDir "C:/BEDFiles" --groupdir "C:/BEDFiles/groups.CancerCompendium.txt" --compdir "C:/BEDFiles/comps.CancerCompendium.txt" --output "C:/BEDFiles" --expname "CancerCompendium" 
```

_**Analyzing CEL files – Exon 1.0 array using default options**_
```
python AltAnalyze.py --species Mm --arraytype exon --celdir "C:/CELFiles" --groupdir "C:/CELFiles/groups.CancerCompendium.txt" --compdir "C:/CELFiles/comps.CancerCompendium.txt" --output "C:/CELFiles" --expname "CancerCompendium"
```

_**Analyzing CEL files – Affymetrix 3’ array using default options and GO-Elite**_
```
python AltAnalyze.py --species Mm --arraytype "3'array" --celdir "C:/CELFiles" --groupdir "C:/CELFiles/groups.CancerCompendium.txt" --compdir "C:/CELFiles/comps.CancerCompendium.txt" --output "C:/CELFiles" --expname "CancerCompendium" --runGOElite yes --returnPathways all
```

_**Analyzing Filtered Expression file – RNA-Seq using custom options**_
```
python AltAnalyze.py --species Mm --platform RNASeq –-filterdir "C:/BEDFiles" --altpermutep 1 --altp 1 --altpermute yes --additionalAlgorithm none --altmethod linearregres --altscore 2 --removeIntronOnlyJunctions yes
```

_**Analyzing CEL files – Exon 1.0 array using custom options**_
```
python AltAnalyze.py --species Hs --arraytype exon --celdir "C:/CELFiles" --output "C:/CELFiles" --expname "CancerCompendium" --runGOElite no --dabgp 0.01 --rawexp 100 --avgallss yes --noxhyb yes --analyzeAllGroups "all groups" --GEcutoff 4 --probetype core --altp 0.001 --altmethod FIRMA --altscore 8 --exportnormexp yes --runMiDAS no --ASfilter yes --mirmethod "two or more" --calcNIp yes
```

_**Analyzing CEL files – HJAY array using custom options**_
```
python AltAnalyze.py --species Hs --arraytype junction --celdir "C:/CELFiles" --output "C:/CELFiles" --expname "CancerCompendium" --runGOElite no --dabgp 0.01 --rawexp 100 --avgallss yes --noxhyb yes --analyzeAllGroups "all groups" --GEcutoff 4 --probetype core --altp 0.001 –altmethod "linearregres" --altscore 8 --exportnormexp yes --runMiDAS no --ASfilter yes --mirmethod "two or more" --calcNIp yes --additionalAlgorithm FIRMA --additionalScore 8
```

_**Analyzing Expression file – Gene 1.0 array using default options, without GO-Elite**_
```
python AltAnalyze.py --species Mm --arraytype gene --expdir "C:/CELFiles/ExpressionInput/exp.CancerCompendium.txt" --groupdir "C:/CELFiles/groups.CancerCompendium.txt" --compdir "C:/CELFiles/comps.CancerCompendium.txt" --statdir "C:/CELFiles/ExpressionInput/stats.CancerCompendium.txt" --output "C:/CELFiles"
```

_**Analyzing Filtered Expression file – Exon 1.0 array using default options**_
```
python AltAnalyze.py --species Hs --arraytype exon --filterdir "C:/CELFiles/Filtered/Hs_Exon_prostate_vs_lung.p5_average.txt" --output "C:/CELFiles"
```

_**Annotate External Probe set results – Exon 1.0 array using default options**_
```
python AltAnalyze.py --species Rn --arraytype exon --annotatedir "C:/JETTA_Results/Hs_tumor_progression.txt" --output "C:/JETTA_Results" --runGOElite yes
```

_**Filter AltAnalyze results with predefined IDs using default options**_
```
python AltAnalyze.py --species Mm --arraytype gene --celdir "C:/CELFiles" --output "C:/CELFiles" --expname "CancerCompendium" --returnAll yes
```


_**Operating System Example Folder Locations**_

  * **PC**: "C:/CELFiles"

  * **Mac OSX**: "/root/user/admin/CELFiles

  * **Linux**: "/hd3/home/admin/CELFiles

When the analysis is finished, the results will be stored just as with the GUI (e.g., alternative exon results in AltResults and gene expression results in ExpressionOutput).

### Primary Analysis Variables ###

No default value for these variables is given and must be supplied by the user if running an analysis. For example, if analyzing CEL files directly in AltAnalyze, you must include the flags ` --species --arraytype --celdir --expname ` and `--output `, with corresponding values. Likewise, when analyzing an existing expression file you must include the flags `--species --arraytype --expdir` and `---output`. Most of the variable values are file or folder locations.

### Universally Required Variables ###

`--platform or --arraytype`
> Platform being analyzed (enter "3'array" if none of the below apply). Values include RNASeq, exon, gene, junction, AltMouse and “3'array”. This variable indicates the general array type correspond to the input CEL files or expression file. An example exon array is the Mouse Affymetrix Exon ST 1.0 array, an example gene array is the Mouse Affymetrix Gene ST 1.0 array and example 3’array is the Affymetrix Mouse 430 version 2.0 array. See the files ArrayFileInfo.txt and array.txt in AltAnalyze/Config/ for details.
`--species`
> Species codes are provided for this variable (e.g., Hs, Mm, Rn). Additional species can be added through the graphic user interface.
`--output`
> Required for all analyses. This designates the directory which results will be saved to.

**Analysis Specific Required Variables**

`--expname`
> Required when analyzing CEL files. This provides a name for your dataset. This name must match any existing groups and comps files that already exist. The groups and comps file indicate which arrays correspond to which biological groups and which to compare. These files must exist in the designated output directory in the folder “ExpressionInput” with the names “groups.**expname**.txt” and “comps.**expname**.txt” where expname is the variable defined in this flag. Alternatively, the user can name their CEL files such that AltAnalyze can directly determine which group they are (e.g., wildtype-1.CEL, cancer-1.CEL, cancer-2.CEL). Go [here](ManualGroupsCompsCreation.md) for more info.
`--celdir --bedDir `
> Required when analyzing CEL files. This provides the path of the CEL files to analyze. These must all be in a single folder.
`--expdir`
> Required when analyzing a processed expression file. This provides the path of the expression file to analyze.
`--statdir`
> Optional when analyzing a processed expression file. This provides the path of the DABG p-value file for the designated expression file to analyze (see –expdir).
`--filterdir`
> Required when aalyzing an AltAnalyze filtered expression file. This provides the path of the AltAnalyze filtered expression file to analyze.
`--cdfdir`
> Required when directly processing some CEL file types. This variable corresponds to the location of the Affymetrix CDF or PGF annotation file for the analyzed array. If you are analyzing an exon, gene, junction, AltMouse or 3’arrays, AltAnalyze has default internet locations for which to download these files automatically, otherwise, you must download the compressed CDF file from the Affymetrix website (support), decompress it (e.g., WinZip) and reference it’s location on your hard-drive using this flag. If you are unsure whether AltAnalyze can automatically download this file, you can try to exclude this variable and see if annotations are included in your gene expression results file.
`--csvdir`
> Required when analyzing some expression files or CEL file types. This variable corresponds to the location of the Affymetrix CSV annotation file for the analyzed array. If you are analyzing an exon, gene, junction, AltMouse or 3’arrays, AltAnalyze has default internet locations for which to download these files automatically, otherwise, you must download the compressed CSV file from the Affymetrix website (support), decompress it (e.g., WinZip) and reference it’s location on your hard-drive using this flag. If you are unsure whether AltAnalyze can automatically download this file, you can try to exclude this variable and see if annotations are included in your gene expression results file.
`--annotatedir`
> Required when annotating a list regulated probe sets produced outside of AltAnalyze. This variable corresponds to the location of the directory containing one or more probe set files. These files can be in the standard JETTA export format, or otherwise need to have probe set IDs in the first column. Optionally, these files can have an associated fold change and p-value (2nd and 3rd columns), which will be reported in the results file.
`--groupdir`
> Location of an existing group file to be copied to the directory in which the expression file is located or will be saved to.
`--compdir`
> Location of an existing comps file to be copied to the directory in which the expression file is located or will be saved to.

**Optional Analysis Variables**

These variables are set as to default values when not selected. The default values are provided in the configurations text file in the Config directory of AltAnalyze (default-.txt) and can be changed by editing in a spreadsheet program.

**GO-Elite Analysis Variables**

AltAnalyze can optionally subject differentially or alternatively expressed genes (AltAnalyze and user determined) to an over-representation analysis (ORA) along Gene Ontology (GO) and pathways (WikiPathways) using the program GO-Elite. GO-Elite is seamlessly integrated with AltAnalyze and thus can be run using default parameters either the graphic user interface or command line. To run GO-Elite using default parameters in command line mode, include the first flag below with the option yes.
`--runGOElite`
> Default is no. Used to indicate whether to run GO-Elite analysis following AltAnalyze. Indicating yes would prompt GO-Elite to run.
`--mod`
> Default is Ensembl. Primary gene system for Gene Ontology (GO) and Pathway analysis to link Affymetrix probe sets and other output IDs to. Alternative values are EntrezGene.
`--elitepermut`
> > Default is 2000. Number of permutation used by GO-Elite to calculate an over-representation p-value. When set to FisherExactTest, Fisher Exact test p-values will be calculated instead.
`--method`

> Default is z-score. Sorting method used by GO-Elite to compare and select the top score of related GO terms. Alternative values are “gene number” and combined
`--zscore`
> Default is 1.96. Z-score threshold used following over-representation analysis (ORA) for reported top scoring GO terms and pathways.
`--elitepval`
> Default is 0.05. Permutation p-value threshold used ORA analysis for reported top scoring GO terms and pathways.
`--dataToAnalyze`
> Default is both. Indicates whether to perform ORA analysis on pathways, Gene Ontology terms or both. Alternative values are Pathways or Gene Ontology
`--num`
> Default is 3. The minimum number of genes regulated in the input gene list for a GO term or pathway after ORA, required for GO-Elite reporting.
`--GEelitepval`
> Default is 0.05. The minimum t-test p-value threshold for differentially expressed genes required for analysis by GO-Elite.
`--GEeliteptype`
> Default is rawp. Indicates whether to run rawp or adjp (Benjamini-Hochberg) p-value for filtering.
`--GEelitefold`
> Default is 2. The minimum fold change threshold for differentially expressed genes required for analysis by GO-Elite. Applied to any group comparisons designated by the user.
`--ORAstat`
> Default is Fisher Exact Test. When the value is set to “Permute p-value”, a permutation p-value will be calculated using the default of provided number of permutations (`--elitepermut`). The adjusted p-value will be calculated from the selected type of ORAstat.
`--additional`
> Default is None. When the value is set to one of a valid resource or “all”, GO-Elite will download and incorporate that resource along with the default downloaded (WikiPathways and Gene Ontology). Additional resources include: “miRNA Targets”, ”GOSlim”, ”Disease Ontology”, ”Phenotype Ontology”, ”KEGG”, “Latest WikiPathways”, ”PathwayCommons”, ”Transcription Factor Targets”, ”Domains” and ”BioMarkers” (include quotes).
`--denom`
> Default is None. This is the folder location containing denominator IDs for corresponding input ID list(s). This variable is only supplied to AltAnalyze when independently using the GO-Elite function to analyze a directory of input IDs (`--input`) and a corresponding denominator ID list.
`--returnPathways`
> Default is None. When set equal to “yes” or “all”, will return all WikiPathways as colored PNG or PDF files (by default both) based on the input ID file data and over-representation results. Default value is “None”. When equal to “top5”, GO-Elite will only produce the top 5 (or other user entered number – e.g., top10) ranking WikiPathways.

**AltAnalyze Expression Filtering and Summarization**

These variables are used to determine the format of the expression data being read into AltAnalyze, the output formats for the resulting gene expression data and filtering thresholds for expression values prior to alternative exon analysis. Since AltAnalyze can process both convention (3’array) as well as RNA-Seq data and splicing arrays (exon, gene, junction or AltMouse), different options are available based on the specific platform.

**Universal Array Analysis Variables**

`--logexp`
> Default is log for arrays and non-log for RNA-Seq. This is the format of the input expression data. If analyzing CEL files in AltAnalyze or in running RMA or GCRMA from another application, the output format of the expression data is log 2 intensity values. If analyzing MAS5 expression data, this is non-log.
`--inclraw`
> Default is yes. When the value of this variable is no, all columns that contain the expression intensities for individual arrays are excluded from the results file. The remaining columns are calculated statistics (groups and comparison) and annotations.
`--vendor`
> Default is Affymetrix. This variable can be set to **"Other ID"** when analyzing data from proteomics, metabolomics or other data not explicitely listed in the vendor/data-type menu. When entering "Other ID", also include the `--platform` as the specific ID system (e.g., PubChem).

**RNASeq, Exon, Gene, Junction or AltMouse Platform Specific Variables**

`--dabgp`
> Default is 0.05. This p-value corresponds to the detection above background (DABG) value reported in the “stats.” file from AltAnalyze, generated along with RMA expression values. A mean p-value for each probe set for each of the compared biological groups with a value less than this threshold will be excluded, both biological groups don’t meet this threshold for a non-constitutive probe set or if one biological group does not meet this threshold for constitutive probe sets.
`--rawexp`
> Default is 0 for microarrays and 1 for RNASeq reads. More stringent examples for Affymetrix arrays are 70 and 3 for RNASeq (RPKM) . For Affymetrix arrays, this value is the non-log RMA average intensity threshold for a biological group required for inclusion of a probe set. The same rules as the `--dabgp` apply to this threshold accept that values below this threshold are excluded when the above rules are not met.
`--avgallss`
> Default is no. For RNA-Seq analyses, default is yes. Indicating yes, will force AltAnalyze to use all exon aligning probe sets or RNA-Seq reads rather than only features that align to predicted constitutive exons for gene expression determination. This option applies to both the gene expression export file and to the alternative exon analyses.
`--runalt`
> Default is yes. Designating no for this variable will instruct AltAnalyze to only run the gene expression analysis portion of the program, but not the alternative exon analysis portion.
`--groupStat`
> Default is moderated t-test. To designate an alternate statistic, this variable can be set to one of the following options are paired t-test, Kolmogorov Smirnov, Mann Whitney U, Rank Sums (these variable should be placed in quotes).

**AltAnalyze Alternative Exon Statistics, Filtering and Summarization**

**Universal Array Analysis Variables**

`--altmethod`
> Default is splicing-index (exon and gene) and ASPIRE (RNASeq, junction and AltMouse). For exon, gene and junction arrays, the option FIRMA is also available and for RNASeq, junction and AltMouse platforms the option linearregres is available.
`--altp`
> Default is 0.05. This variable is the p-value threshold for reporting alternative exons. This variable applies to both the MiDAS and splicing-index p-values.
`--probetype`
> Default is core (exon and gene) and all (junction and AltMouse). This is the class of probe sets to be examined by the alternative exon analysis. Other options include: extended and full (exon and gene) and “exons-only”, “junctions-only”, “combined-junctions” (RNASeq, junction and AltMouse).
`--altscore`
> Default is 2 (splicing-index) and 0.2 (ASPIRE).This is the corresponding threshold for the default algorithms listed under {{{--altmethod.
`--GEcutoff`
> Default is 3. This value is the non-log gene expression threshold applied to the change in gene expression (fold) between the two compared biological groups. If a fold change for a gene is greater than this threshold it is not reported among the results, since gene expression regulation may interfere with detection of alternative splicing.
`--analyzeAllGroups`
> Default is pairwise. This variable indicates whether to only perform pairwise alternative exon analyses (between two groups) or to analyze all groups, without specifying specific comparisons. Other options are “all groups” and both.
`--altpermutep`
> Default is 0.05. This is the permutation p-value threshold applied to AltMouse array analyses when generating permutation based alternative exon p-values. Alternative exon p-values can be applied to either ASPIRE or linregress analyses.
`--altpermute`
> Default is yes. This option directs AltAnalyze to perform the alternative exon p-value analysis for the AltMouse array (see –altpermutep).
`--exportnormexp`
> Default is no. This option directs AltAnalyze to export the normalized intensity expression values (feature expression/constitutive expression) for all analyzed features (probe sets or RNA-Seq reads) rather than perform the typical AltAnalyze analysis when its value is yes. For junction-sensitive platforms, rather than exporting the normalized intensities, the ratio of normalized intensities for the two reciprocal-junctions are exported (pNI1/pNI2). This step can be useful for analysis of exon array data outside of AltAnalyze and comparison of alternative exon profiles for many biological groups (e.g., expression clustering).
`--runMiDAS`
> Default is yes. This option directs AltAnalyze to calculate and filter alternative exon results based on the MiDAS p-value calculated using the program Affymetrix Power Tools.
`--calcNIp`
> Default is yes. This option directs AltAnalyze to filter alternative exon results based on the t-test p-value obtained by comparing either the normalized intensities for the array groups examined (e.g., control and experimental) (splicing-index) or a t-test p-value obtained by comparing the FIRMA scores for the arrays in the two compared groups.
`--mirmethod`
> Default is one. This option directs AltAnalyze to return any microRNA binding site predictions (default) or those that are substantiated by multiple databases (two or more).
`--ASfilter`
> Default is no. This option directs AltAnalyze to only analyze probe sets or RNA-Seq reads for alternative expression that have an alternative-splicing annotation (e.g., mutually-exclusive, trans-splicing, cassette-exon, alt-5’, alt-3’, intron-retention), when set equal to yes.
`--returnAll`
> Default is no. When set to yes, returns all un-filtered alternative exon results by setting all associated filtering parameters to the lowest stringency values. This is equivalent to providing the following flags: `--dabgp 1 --rawexp 1 --altp 1 --probetype full --altscore 1 --GEcutoff 10000 `. Since this option will output all alternative exon scores for all Ensembl annotated RNA-Seq reads or probe sets, the results file will be exceptionally large (>500,000 lines), unless the user has saved previously run alternative exon results (e.g., MADS) to the directory “AltDatabase/filtering” in the AltAnalyze program directory, with a name that matches the analyzed comparison. For example, if the user has a list of 2,000 MADS regulated probe sets for cortex versus cerebellum, then the MADS results should be saved to “AltDatabase/filtering” with the name “Cortex\_vs\_Cerebellum.txt” and in AltAnalyze the CEL file groups should be named Cortex and Cerebellum and the comparison should be Cortex versus Cerebellum. When the filename for a file in the “filtering” directory is contained within the comparison filename (ignoring “.txt”), only these AltAnalyze IDs or probe sets will be selected when exporting the results. This analysis will produce a results file with all AltAnalyze statistics (default or custom) for just the selected features, independent of the value of each statistic.
`--additionalAlgorithm`
> Default is splicing-index. For Affymetrix arrays, setting this flag equal to FIRMA changes the individual probe set analysis algorithm from splicing-index to FIRMA for junction arrays. This method is applied to RNA-Seq data and junction arrays following reciprocal-junction analysis (e.g., ASPIRE) in a second run. To exclude this feature, set variable equal to none.
`--additionalScore`
> Default is 2. Setting this flag equal to another numeric value (range 1 to infinity) changes the non-log fold change for the additional\_algorithms.
`--removeIntronOnlyJunctions`
> Default is no. Indicates whether to remove junctions where both splice sites align to outside of annotated exons. Setting this value to yes will remove these putative junctions prior to analysis.
`--normCounts`
> Default is RPKM. Indicates whether to normalize exon and/or junctions counts using the methods RPKM. Setting this value to none will use the original counts for gene expression and alternative exon analyses. Setting this variable to none, will directly analyze counts with no normalization.
`--buildExonExportFile`
> Default is no. Indicates whether to halt an RNA-Seq analysis after reading in the junction.bed files and export a new file with the suffix exon.bed in the folder BAMtoBED. This file is used as described [here](BAMtoBED.md) to obtain exon.bed files for all experimental BAM files. Setting this value to yes will use the original counts for gene expression and alternative exon analyses.
`--groupStat`
> Default is “moderated t-test”. Indicates the algorithm to employ for all pairwise group comparisons in AltAnalyze. Other options include: “paired t-test”, “Kolmogorov Smirnov”, “Mann Whitney U” and “Rank Sums”.
`--exonRPKM`
> Default is 0.5. Numerical filter used to threshold exons as expressed when using RPKM normalization for mean RPKM values for each biological group. At least one biological group for a pairwise comparison must meet this threshold to be included in further analyses. Increasing values of this filter will increase the stringency of exon/gene expression.
`--geneRPKM`
> Default is 1. Numerical filter used to threshold exons as expressed when using RPKM normalization for mean RPKM values for each biological group. At least one biological group for a pairwise comparison must meet this threshold to be included in further analyses. Increasing values of this filter will increase the stringency of exon/gene expression.
`--exonExp`
> Default is 10. Numerical filter used to threshold genes or exons as expressed based on the mean read counts values for each biological group. At least one biological group for a pairwise comparison must meet this threshold to be included in further analyses. Increasing values of this filter will increase the stringency of exon/gene expression.

**AltAnalyze Database Updates**

**Universal Array Analysis Variables**

`--update`
> Default is empty. Setting this flag equal to Official, without specifying a version, will download the most up-to-date database for that species. Other options here are used internally by AltAnalyze.org for building each new database. See the method “commandLineRun” in AltAnalyze.py for more details.
`--version`
> Default is current. Setting this flag equal to a specific Ensembl version name (e.g. EnsMart49) supported by AltAnalyze will download that specific version for the selected species, while setting this to current will download the current version.
`--specificArray`
> Default is none. Indicates a sub-type of a particular array platform (e.g., junction – HJAY or hGlue) when building the database. This variable only needs to be set when currently building the hGlue junction database (see [BuildingDatabases](BuildingDatabases.md)).
`--ignoreBuiltSpecies`
> Default is none. Indicating yes will only build species databases for species without already built species directories. This is used during internal database release building to simultaneously build multiple species databases.

**Additional Analysis, Quality Control and Visualization Options**

In addition to the core AltAnalyze workflows (e.g., normalization, gene expression summarization, evaluation of alternative splicing), additional options are available to evaluate the quality of the input data (quality control or QC), evaluate associated cell types and tissues present in each biological sample (Lineage Profiler), cluster samples or genes based on overall similarity (expression clustering) and view regulation data on biological pathways (WikiPathways). These options can be run as apart of the above workflows or often independently using existing AltAnalyze results or input from other programs.

`--outputQCPlots`
> Default is yes. Instructs AltAnalyze to calculate various QC measures specific to the data type analyzed (e.g., exon array, RNASeq) when a core workflow is run. This will include hierarchical clustering and PCA plots for genes considered to be differentially expressed (see `--GEelitefold` and `--GEelitepval`). Outputs various QC output plots (PNG and PDF) to the folder “DataPlots” in the user defined output directory. If run from python source code, requires installation of Scipy, Numpy and Matplotlib.
`--runLineageProfiler`
> Default is yes. Instructs AltAnalyze to calculate Pearson correlation coefficients for each analyzed user sample relative to all cell types and tissues in the BioMarker database. Resulting z-scores for calculated from the coefficients are automatically visualized on a WikiPathways Lineage network and are hierarchically clustered. Outputs various tables to the folder ExpressionOutput and plots (PNG and PDF) to the folder “DataPlots” in the user defined output directory. Can be run as apart of an existing workflow or independently with the option `--input`. If run from python source code, requires installation of suds, Scipy, Numpy and Matplotlib.
`--input`
> Default is None. Including this option indicates that the user is referencing an expression file location that is supplied outside of a normal AltAnalyze workflow. Example analyses include only performing hierarchical clustering, PCA or GO-Elite.
`--image`
> Default is None. Including this option with any of the variables: WikiPathways, hierarchical, PCA, along with an input expression file location (`--input`), will prompt creation and export of the associated visualization files. These analyses are outside of the typical AltAnalyze workflows, only requiring the designated input file. If run from python source code, requires installation of suds, Scipy, Numpy and Matplotlib.

**Hierarchical Clustering Variables**

`--row_method`
> Default is average. Indicates the cluster metric to be applied to rows. Other options include: average, single, complete, weighted and None. None will result in no row clustering.
`--column_method`
> Default is single. Indicates the cluster metric to be applied to columns. These options are the same as row\_method.
`--row_metric`
> Default is “row\_metric”, cosine. Indicates the cluster distance metric to be applied to rows. Other options include: braycurtis, canberra, chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard, kulsinski, mahalanobis, matching, minkowski, rogerstanimoto, russellrao, seuclidean, sokalmichener, sokalsneath, sqeuclidean and yule (not all may work). If the input metric fails during the analysis (unknown issue with Numpy), euclidean will be used instead.
`--column_metric`
> Default is euclidean. Indicates the cluster distance metric to be applied to rows. These options are the same as row\_metric.
`--color_gradient`
> Default is red\_white\_blue. Indicates the color gradient to be used for visualization as up-null-down. Other options include: red\_black\_sky, red\_black\_blue, red\_black\_green, yellow\_black\_blue, green\_white\_purple, coolwarm and seismic.
`--transpose`
> Default is False. Will transpose the matrix of columns and rows prior to clustering, when set to True.

`--compendiumPlatform`
`--compendiumType`
`--exonMapFile`
`--geneExp`
`--labels`
`--contrast`
`--plotType`
`--runMarkerFinder`
`--update_interactions`
`--includeExpIDs`
`--degrees`
`--genes`
`--inputType`
`--interactionDirs`
`--GeneSetSelection`
`--PathwaySelection`
`--OntologyID`
`--dataType`
`--combat`
`--channelToExtract`
`--showIntrons`
`--display`
`--join`
`--uniqueOnly`
`--accessoryAnalysis`
`--inputIDType`
`--outputIDType`
`--FEdir`
`--channelToExtract`
`--AltResultsDir`
`--geneFileDir`
`--modelSize`
`--geneModel`
`--reference`

`python AltAnalyze.py --expdir /MyExperiment/ExpressionInput/exp.tumors.txt --update markers --species Hs --platform RNASeq --genesToReport 200`