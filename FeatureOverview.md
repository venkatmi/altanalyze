# Major Features in AltAnalyze #

The below items are a current list of available features in AltAnalyze as of version 2.0.9.

## Interfaces ##

  * Graphical User-Interface - AltAnalyze
  * Graphical User-Interface - [AltAnalyze Results Viewer](InteractiveResultsViewer.md)
  * [Command-line AltAnalyze](CommandLineMode.md)

## Major Features ##

  * Easy-to-use interfaces with support for many species and platforms
  * RNA-Seq direct sequence and aligned results analysis
  * Single-cell transcriptome discovery analysis
  * Conventional and splicing microarray (Affymetrix) raw file and processed analysis
  * Comprehensive annotation databases (Ensembl, UCSC, miRBase, WikiPathways, GeneOntology, UniProt, others)
  * Basic and advanced statistics for gene expression analysis (empirical Bayes limma t-test, Benjamini-Hochberg, q-value, others)
  * Basic and advanced normalization methods (quantile, RMA, Combat, RPKM, TPM)
  * Advanced methods for splicing determination
  * Advanced methods for splicing functional inference prediction
  * Cell and tissue-type sample prediction analysis (LineageProfiler)
  * Basic and advanced expression clustering and visualization
  * Basic and advanced gene-set enrichment analysis (GO-Elite)
  * Protein, metabolite and gene-regulatory network analysis (NetPerspective)
  * Group specific marker identification analysis and visualization (MarkerFinder)
  * Sample classification and feature selection (LineageProfilerIterate)
  * Quality control analysis and visualization
  * Principal component (2D and 3D) and t-SNE analysis and visualization
  * Splicing event visualization (AltExonViewer, SashimiPlot)
  * Interfaces to R (HOPACH clustering, Monocle, others)
  * Identify translation, between any two supported systems
  * Venn Diagram visualization and interactive gene selection
  * WikiPathway visualization (via REST interface)
  * Interactive merging of text files
  * Interactive results visualization of AltAnalyze produced tables and image files (InteractiveResultsViewer)

### Databases ###

  * Over 60 species databases (varying per database version), centralized around Ensembl versions (refereed to as EnsMart).
  * Multiple versions of Ensembl supported (updated 1-2 times per year)
  * Ensembl database augmented with annotations from UCSC (mRNA isoform and splicing), UniProt (isoform, domain and functional motif), microRNA (various databases, sequence-level targets), WikiPathways and others.
  * LineageProfiler gene marker and expression database.
  * GO-Elite gene marker database.
  * Downloadable versions from the GUI or command-line.

### RNA-Seq Initial Processing ###

  * Read pseudo-alignment and gene/isoform quantification (TPM) via direct integration Kallisto (local version included with the software). FASTA files automatically downloaded for the associated version of Ensembl.
  * BAM file import support for TopHat, STAR and others. If the XS tag is missing from the BAM, AltAnalyze will use the locally installed database to filter and identify junction reads in a strand-aware manner.
  * BED file import support. ExonBed and JunctionBed formats (TopHat, SpliceMap).
  * BioScope exon and junction text file support.
  * TCGA junction count files (see documentation).
  * Exon and junction-level quantification and gene-level RPKM calculation for output from TopHat or STAR (BED or BAM), BioScope or TCGA.
  * Gene and isoform (Ensembl transcript, RefSeq ID) support for pre-processed expression files produced from outside tools (e.g., RSEM, Express, Cufflinks).

### Microarray Initial Processing ###

  * Affymetrix splicing sensitive microarray support (CEL files) with pre-built gene and isoform annotations included. RMA normalization is applied for exon, junction and gene-level quantification along with DABG p-value calculation.
  * Gene-level quantification for Affymetrix and Agilent mRNA arrays.
  * Pre-normalized text file analysis.
  * Expression normalization for any tabular dataset using quantile normalization or Combat.
  * Log2 transformation.

### Gene Expression Analysis ###

  * Calculation of fold and comparison statistics. Non-adjusted pair-wise group p-values calculated using multiple available options and adjusted p-values using Benjamini-Hochberg or q-values. ANOVA group statistics calculated comparing all annotated groups.
  * Expression filtering options at the gene or exon-level.
  * Gene annotations using multiple databases (ExpressionOutput/DATASET file).
  * Single-cell transcriptome filtering, clustering, guide-gene selection and visualization with the novel algorithm [ICGS](ICGS.md). Integrated automated cell-type prediction analyses. Sample-group and pattern discovery from non-single cell data also supported.
  * Single-cell transcriptome prediction analysis from the AltAnalyze [Heatmaps](Heatmaps.md) interface.
  * Single-cell transcriptome and loading gene prediction analysis from the AltAnalyze [PCA](PCA.md) interface.
  * Gene-set, pathway and Ontology enrichment analyses using GO-Elite (automated and manual run options).
  * Marker prediction for individual all genes and annotated groups and associated heatmap visualization using the novel algorithm MarkerFinder.
  * Gene-level expression graphs using the InteractiveResultsViewer.
  * Network visualization analysis from a priori identified gene lists, regulated genes identified from AltAnalyze, interaction files, pathways, Ontologies and more using the novel NetPerspective algorithm. Multiple algorithms and interaction databases provided.
  * Venn Diagram analysis of produced result files (e.g., GO-Elite folder differential expressed gene inputs or custom files.
  * Pathway visualization upon any existing WikiPathways, including those you just made (automatically from [GO-Elite](GOElite.md) or manually).


### Alternative Splicing Analysis ###

  * Exon-level splicing-index/FIRMA analysis for any exon-sensitive profiling platform (Affymerix exon array, junction array, RNA-Seq)
  * Reciprocal junction splicing analyses for RNA-Seq and junction arrays (ASPIRE, LinearRegression, junction-centric Percent Spliced In (PSI)).
  * Splicing-event type annotation (e.g., cassette-exon, intron retention).
  * Domain, protein motif, protein length and microRNA binding site disruption prediction analyses.
  * Integrative splicing-index and reciprocal junction results to identify extremely high-confidence events.
  * MIDAS p-value calculation.
  * Domain and microRNA-binding site disruption enrichment analyses
  * GO-Elite analysis.
  * Automated heatmap visualization.
  * DomainGraph visualization.
  * Splicing graph (AltExonViewer) and Sashimi-Plot visualization.
  * NetPerspective analysis.
  * Isoform differential expression analysis (Kallisto output).
  * ICGS analysis of PSI results.

### Gene-Set Enrichment Analysis ###
  * Pruned ontology results sets (GeneOntology, DiseaseOntology, PhenotypeOntology).
  * Multiple pathway sources (WikiPathways, Reactome, PathwayCommons, KEGG).
  * Multiple regulatory prediction analyses (PAZAR/Amadeus experimental transcription factor, WholeGenomeRVista transcription factor, microRNA binding sites detected by two or more algorithms).
  * LineageProfiler BioMarker enrichment analysis.
  * Automated visualization of enriched WikiPathways.
  * Automated gene-to-term network visualization.
  * Automated comparison of enrichment results via hierarchical clustering and heatmap visualization.

### Advanced Clustering Analyses ###
  * Hierarchical clustering using all default SciPy associated methods and metrics.
  * Automated connection to R and use of HOPACH clustering for improved cluster delineation.
  * Multiple color gradients, scaling and normalization options.
  * Drill-down analysis within java TreeView (must have JAVA installed).
  * Automated coloring of annotated biological groups (groups. file in ExpressionInput folder).
  * Gene-set filtering based on all GO-Elite available sets, in addition to different protein-types (transcription factor, kinase, GPCR, cell surface, etc.).
  * Gene-set driven pattern discovery using these gene-sets or custom supplied gene(s) and correlation cutoffs. Can identify correlated and anti-correlated, or just correlated.
  * Filtering of clustered inputs for all highly intra-correlated genes only.
  * Interactive web linkouts to gene annotation databases from the heatmap.
  * Within cluster automated GO-Elite enrichment analysis reporting with interactive linkouts to the gene-to-term interaction network and associated gene lists.
  * Automated [Monocle](Monocle.md) pseudotemporal ordering analysis via a connection to R.
  * Store results into your local gene-set library for future clustering or GO-Elite analyses.
  * Perform [ICGS](ICGS.md) on an input expression dataset (unsupervised), with or without additionally supplied genes or gene-sets (supervised).


### Other Analyses ###
  * See [Tutorials](Tutorials.md).