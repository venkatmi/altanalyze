# News and Updates #

Proceed [here](UpdateHistory.md) for detailed software version update information.

## 8-23-15 ##

We have migrated this Wiki and code repository to [Github](https://github.com/nsalomonis/altanalyze). However, documentation on this site remains accessible (specifically for version less than AltAnalyze 2.0.9.1).

## 8-21-15 ##

A radically improved version of AltAnalyze, version 2.0.9.1, is now available. This software includes improved algorithms for single-cell expression analysis [ICGS](ICGS.md), a new standalone [Interactive Results Viewer](InteractiveResultsViewer.md), direct pseudo-alignment and quantification from FASTQ files (single or paired-end) using [kallisto](https://github.com/pachterlab/kallisto), a direct BAM file interface, support for Sashimi-Plots built in, new statistical methods (stringent PSI algorithm, q-value), pseudo-temporally ordering of samples with [Monocle](http://bioconductor.org/packages/release/bioc/html/monocle.html) via a connection to R and more. For details please see the UpdateHistory.

## 3-22-15 ##
A new version of AltAnalyze, version 2.0.9, is now available with new algorithms, improved functionality and significantly enhanced features. This version represents a significant advance upon previous version. Major highlights of this version are 1) multiple independent tools and workflows for single-cell gene expression and splicing analyses (Iterative Clustering and Driver Selection or ICDS), 2) major improvements to the heatmap viewer with integrated pathway enrichment, visualization options, advanced correlation methods and interactive gene and TreeView linkouts, 3) a new splicing analysis algorithm called Reciprocal-isoform Percent Spliced In (Ri-PSI) provides improved accuracy in alternative splicing and promoter identification., 4) support for Affymetrix HTA2.0 arrays, 4) integrated support for TCGA junction expression files, 5) interactive Venn Diagrams and the ability to store gene signatures, 6) improved results menu, 7) improved PCA, MarkeFinder, network analysis, enrichment analysis, alternative exon visualization and clustering features. A number of important bug fixes have also been made.

Additional details, [tutorials](Tutorials.md) and [YouTube videos](YouTubeVidoes.md) are available from http://www.altanalyze.org and in the software documentation directory.

## 9-18-14 ##
We have released a prototype version of AltAnalyze that is still in development for the de novo identification of sample groups (aka single-cell analysis). A challenging bioinformatics challenge is the unbiased identification of both major and sub-populations of bulk or single-cell transcriptomes where no obvious reference gene or splice signatures exist (see [tutorial 4](https://code.google.com/p/altanalyze/wiki/Tutorial_De_Novo_SampleAnalysis)). In addition, heatmap visualization of splicing gene data is now available along with support for TCGA input junction files and multiterm selection options from the heatmap viewer and NetPerspective. Please let us know if you have challenges with this version as we are still testing it. To download go to:
http://sourceforge.net/projects/altanalyze/

## 8-13-13 ##
We are happy to announce several new features to AltAnalyze for the advanced analysis of gene expression and splicing data. AltAnalyze version 2.0.8 provides a several new analyses for functional interpretation of your results.  These include a new tool analysis called NetPerspective, batch effects normalization with combat, splicing graph visualization, supervised clustering from genes or GO-Elite gene sets, new sample classification analyses, support for the import and normalization of Agilent Feature Extraction format files, Venn diagram visualization and parallel processing options.  NetPerspective allows for creation of networks from differentially expressed set of genes, pathways, ontologies, existing networks and input gene lists. These results can be visualized with existing data (GO-Elite input format files) using multiple interaction algorithms, via directed interactions from multiple useful resources (WikiPathways, KEGG, PAZAR, Amadeus, DrugBank, BioGRID). PyCombat provides support for batch effect normalization where significant confounding variables exist. For more information on these resources see our [tutorials page](https://code.google.com/p/altanalyze/wiki/Tutorials) and [Wiki](https://code.google.com/p/altanalyze/w/list).

## 7-7-12 ##
We are happy to announce the release of a substantially improved version of AltAnalyze. AltAnalyze version 2.0.7 provides a wide-array of new supporting analysis methods to help improve biological interpretation of gene-level results. New methods include hierarchical clustering, principal component analysis, various QC metrics, expanded gene-set enrichment features and cell lineage analysis. Each of these analyses includes associated data visualization (plots and pathways). A new moderated pairwise t-test method is introduced, replicating the popular method found in the R package limma. Lineage analysis is performed using a newly integrated method called LineageProfiler, capable of assessing relative populations of cell types. For more information see the program [documentation](http://www.altanalyze.org/help_main.htm), [update history](UpdateHistory.md) and our [new blog](http://altanalyze.blogspot.com/).

## 2-19-12 ##
Multiple Ensembl release 65 databases (human, mouse, rat and chicken) were modified and re-posted. Each was missing alternative polyadenylation annotations that were not properly included in the 1-18-12. To update we recommend deleting the folder EnsMart65 in the folder AltDatabase. We apologize for any inconvenience.

## 1-19-12 ##
Posted AltAnalyze version 2.0.6. This version includes additional improvements to RNA-Seq gene-level RPKM calculations to provide more accurate gene expression estimates, new filtering options to differentiate between expressed and non-expressed genes, exons or junctions, support for the Affymetrix hGlue splicing platform, inclusion of all analyzed splicing scores, new fully automated methods for AltAnalyze database creation and fixes to several bugs introduced in version 2.0.5. For more information see the program [documentation](http://www.altanalyze.org/help_main.htm) and [update history](UpdateHistory.md).

## 1-18-12 ##
New species databases posted for Ensembl release 65. Introduced in this version are improved methods for constitutive exon identification, inclusion of additional Affymetrix probesets and addition of splicing factor annotations.

## 12-19-11 ##
A [new study](CitingArticles.md) uses AltAnalyze to identify alternative mRNA isoform biomarkers unique to Parkinson's disease blood leukocytes.

## 11-27-11 ##
Posted AltAnalyze version 2.0.5. This version includes improvements to RNA-Seq gene-level RPKM calculations to provide more accurate gene expression estimates, additional annotations for conventional array analyses and miRNA array, selection of multiple statistical p-value methods added and critical bug fixes. For more information see the program [documentation](http://www.altanalyze.org/help_main.htm) and [update history](UpdateHistory.md).

## 7-9-11 ##
Posted AltAnalyze version 2.0.4 with critical bug fixes to version 2.0.3. For more information see the program [documentation](http://www.altanalyze.org/help_main.htm) and [update history](UpdateHistory.md).

## 6-17-11 ##
Posted AltAnalyze version 2.0.3. This version is compatible with high-throughput RNA sequencing (RNASeq) exon alignment results in addition to junction results. Several new database and annotation features have been added (e.g., alternative polyadenylation, nonsense mediated decay) along with new analysis methods, program defaults, bug fixes and memory optimization. For more information see the program [documentation](http://www.altanalyze.org/help_main.htm) and [update history](UpdateHistory.md).

## 6-7-11 ##
A [new study](CitingArticles.md) uses AltAnalyze to assess the effect of alternative splicing on protein domains and miRNA binding sites in Kaposi's sarcoma-associated herpesvirus-infected lymphatic endothelial cells.

## 4-18-11 ##
A [new study](CitingArticles.md) uses AltAnalyze to understand the role of alternative splicing in distinct neuroblastoma stages.

## 3-17-11 ##
New Windows 32 bit and 64 bit AltAnalyze version 2.0 beta were posted. These replace archives which had an outdated Affymetrix Power Tools directory structure (version 1.4) that will result in a program error when analyzing Affymetrix CEL files on Windows. See the following thread for [details](http://groups.google.com/group/alt_predictions/browse_thread/thread/31ab010a489ac4d9).

## 3-5-11 ##
Posted AltAnalyze version 2.0 beta, with improved data import, increased efficiency for datasets with many novel junction predictions (100k plus) and options for removing junctions with both ends mapping to introns. To improve performance and prediction success, reciprocal junctions where both junctions align only to introns or UTRs and no exons are no longer considered.

## 2-10-11 ##
Introducing AltAnalyze version 2.0 alpha! AltAnalyze 2.0 is compatible with the analysis of high-throughput RNA sequencing (RNASeq) experimental data. AltAnalyze will import and analyze BED format junction files containing genomic alignment coordinates and read counts from two or more biological groups. These files can be produced from multiple software tools ([TopHat](http://tophat.cbcb.umd.edu/), [HMMSplicer](http://derisilab.ucsf.edu/index.php?software=105), [SpliceMap](http://www.stanford.edu/group/wonglab/SpliceMap)) to identify reciprocal exon-junctions that display evidence of alternative exon inclusion in two or more experimental conditions. In addition to known splice junctions (Ensembl and UCSC), novel reciprocal junctions and trans-splicing events are also analyzed and annotated (e.g., alternative cassette-exon). For more information see the associated [ReadMe](http://altanalyze.org/AltAnalyze_Manual.pdf) or the online analysis [tutorial](Tutorial_AltExpression_RNASeq.md).

## 1-09-11 ##
Posted an update to AltAnalyze version 1.15 (version 1.155). Provides solutions to reported bugs that effected program usability and performance. Issues addressed include: 1) Improved performance of FIRMA algorithm, 2) prevents crash due to conflicting parameters, 3) fixes missing junction array protein/domain annotations, 4) incorrect constitutive fold reported for multigroup FIRMA and SI analysis, 5) better handling of downloading errors, 6) corrected t-test method for splicing-index. None of these issues adversely effect previous AltAnalyze results but do enhance performance and annotation. Also posted a beta version (available under downloads) of version 1.16 which analyzes mouse and human junction (JAY) arrays.

## 12-28-10 ##
An AltAnalyze database release for Ensembl version 60 is now available. In this version, comparison between thousands of new mRNA and protein isoforms have been added relative to the previous release, providing new functional predictions. Such annotations result in reassignment of previous extended or full probeset annotations to core. Probeset-to-miRNA binding site ssociations from the microRNA databases miRanda and TargetScan have also been updated in this release. Corresponding GO-Elite database updates for Ensembl version 60 are included. To update your existing database, simply select "Get new species/vendor/array databases online" in the initial AltAnalyze Main Dataset Parameter Window.

## 10-6-10 ##
The AltAnalyze website has been updated to include links to a new support site hosted at http://code.google.com/p/altanalyze. This site includes wiki documentation, an SVN repository, downloads and a bug tracking system.

## 10-1-10 ##
Posted an updated version of 1.15 that had a bug resulting in problems with the ["Annotate External Results"](http://code.google.com/p/altanalyze/wiki/AnnotateExternal) option on 9-27-10, however, a subsequent update was performed on 10-1-10, since an important line of code was missing from the 9-27-10 version. Thus, AltAnalyze 1.15 PC and Linux installations downloaded during that period are faulty.

## 9-10-10 ##
AltAnalyze was the featured application for the September issue of the [Canadian Bioinformatics Help Desk Newsletter](http://gchelpdesk.ualberta.ca/news/15sep10/cbhd_news_15sep10.php#spotlight). This monthly newsletter is intended to keep Genome Canada researchers and other Help Desk users informed about new software, events, job postings, conferences, training opportunities, interviews, publications, awards, and other newsworthy items concerning bioinformatics, computational biology, genomics, proteomics, metabolomics, transcriptomics, systems biology, and synthetic biology.

## 6-05-10 ##
AltAnalyze is featured in the list of Science's [ST NetWatch: Bioinformatics Resources](http://stke.sciencemag.org/cgi/ul/sigtransUl%3BCAT_2).

## 4-12-10 ##
Introducing AltAnalyze version 1.15 with expanded features for alternative exon array analysis. New to AltAnalyze; 1) Alternative exon analysis of Affymetrix Gene 1.0 ST arrays, 2) FIRMA analysis, 3) Bundled version of Cytoscape and DomainGraph included, 4) AltAnalyze DomainGraph start option, 5) Import and filtering based on 3rd party alternative exon analysis results and 6) Linux exectuable version. Please contact us with any questions.

## 2-01-10 ##
Posted a significant update to AltAnalyze (beta 1.14). Although this version also uses database version EnsMart54 (optional), the referenced database has also been updated specifically for this version (version 1.13 downloads a distinct EnsMart54 database which has not changed). This new database includes new constitutive probeset predictions that will cause the alternative exon scores to differ than version 1.13, as a result of constitutive fold change estimates based on a slightly different set of constitutive probesets. These constitutive probesets are now derived in AltAnalyze, by determining which exon-regions are most common to all Ensembl and UCSC database mRNAs rather than relying on Affymetrix probeset mRNA counts (see Section 3.4 and Section 6.2 of the online documentation). Additional changes include; (1) Improved filtering methods for gene expression data, (2) Support for core probe set only gene expression calculation, 3) Addition of advanced gene expression analysis statistics, 4) Support for multi-group (>2) alternative exon analyses, 5) GO-Elite support for adjusted p-value filtering optional settings, 6) Improved automatic array type identification and 7) More streamlined graphical user interface options. For more details, click here.

## 12-23-09 ##
A single Mac OS X installer is now posted. Due to python default installation differences on Mac OS X, separate OS X version installers were previously required, however, a compiling configureation on Leopard (10.5) with Python 2.4.3 was found sufficient to build a cross Mac OS X build. A previously un-identified OS X 10.5 bug crashing or preventing folder browsing was also repaired.

## 12-10-09 ##
Posted an update to AltAnalyze (release 1.13) that includes; (1) versioned database releases based on Ensembl-Affymetrix release cycle and archival updates, (2) improved annotation support for non-Affymetrix arrays for gene expression analysis, (3) addressed unix file related errors, (4) improved user interface (remembering prior directory selections, error handling), (5) introudced a command-line version of AltAnalyze for python source code remote server use, (6) google groups issues resolved.s been causing problems.

## 6-23-09 ##
Posted an update to AltAnalyze (beta 1.12) that includes; (1) built-in GO-Elite pathway over-representation analysis options, (2) fixed linux directory read error and UI window size issues, (3) fixed 'back' button error that caused AltAnalyze to exit after establishing options, (4) posted a more intuitive linux version, (5) exception for result folder opening in Linux is properly handled.

## 5-29-09 ##
Posted an update to AltAnalyze (beta 1.11) that includes additional improvements to the new UI, bug fixes, integrated help, a new database build (synched with version 52 of Ensembl) and automated download of Affymetrix library and annotation files.

## 4-29-09 ##
Posted a substantial update to AltAnalyze (beta 1.1) that adds new functionality and improved methods to AltAnalyze. These include running on RMA on CEL files (using bundled Affymetrix Power Tools modules), an improved interface, ability to define groups and comparisons in the UI, specialized export files for GenMAPP and DomainGraph, improved error handling and log files and improved domain-level predictions and more (see Feature Request page).

## 10-29-08 ##
Posted an update to AltAnalyze (beta 1.01) that addresses specific bugs. These include processing of conventional gene expression microarrays (e.g., Affymetrix 3' and gene arrays). Also, better handling of errors has been implemented for comparison ("comps.") files where the headers do not match those in the expression ("exp.") file. Posted a stand-alone Mac version that no longer requires a Terminal interface to run, but is rather a standard Mac OS X application. NOTE: The update function is still not implemented within the GUI, however, is available through the command line via the "update.py" module. Update should be implemented as apart of the GUI in the next version.

## 7-20-08 ##
Posted the first public version of AltAnalyze (beta 1.0). This software currently supports Affymetrix Exon 1.0 ST array analyses (exon-level) and the custom Affymetrix splicing array AltMouse (exon and junction-level). See Overview and Instruction information or the Google Groups page for more information. Database is based on build 49 of Ensembl.