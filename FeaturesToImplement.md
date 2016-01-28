## Features to Implement in AltAnalyze version 2+ ##

### Core Algorithms ###

  * **Expression Import** - Add probe-level support from Ensembl, for non-splicing Affy arrays to analyze similar to RNA-Seq data. Allow users to supply their own probe-level annotation file for Affy and non-Affy arrays when probe-level expression values are provided (genomic coordinates required). This will effectively build the current exon.bed files. Include the original IDs as annotations in the final tables (or create simple translation tables).
  * **Lineage Profiler** - Add additional methods and UI control for LineageProfiler. Methods include: presence calls for all splicing-sensitive platforms (gene-level), optional exon-level analysis, possibly average tissue signature and more advanced statistics. Add documentation for building the databases (for users to add species support). Add translation support across species.
  * **Built-in BAMtoBED via C** - Integrate the BAMtoBED conversion step for creating exon.bed files into AltAnalyze. This would work by selecting a directory of directories containing TopHat results for different samples.
  * **Paired-end support** - Pre-compute and add junction counts based on paired ends that flank junctions (assuming they flank a junction within a known range given the mate-pair distance and no known conflicts). This would be done at the BAMtoBED-level.
  * **Transcript-Level analysis** - Integrate predictions from other applications (Cufflinks, NSMAP, etc.) and develop a unique approach for AltInteract that uses unique exons to determine likely isoforms that are expressed or not in a given sample or differentially expressed.
  * **Motif over-representation** - Include ESE/ESS/ISE/ISS prediction enrichment methods.
  * **Code re-design** - more efficient methods throughout (using external libraries optionally to improve performance), module reorganization, more efficient database access (MySQL, Derby, other), python 3+ compatibility (http://docs.python.org/2/glossary.html#term-to3).
  * **Robust cross-platform network visualization.**
  * **miRNA and splicing** - Extension of miRNA exon predictions to alternative isoforms (not just the regulated exon).
  * **ncRNA database integration** - Include a comprehensive ncRNA and lincRNA database for gene expression studies.
  * **Automatic genome version detection** - Detect which genome version is being analyzed or look for potential mismatches and flag the user immediately if the % of mapped reads are not isoform aligned.

### Database Builds ###

  * Address [issue 16](https://code.google.com/p/altanalyze/issues/detail?id=16)
  * Add reference isoform annotations (RefSeq? UCSC? Other?)

### Visualization ###

  * **Clustering** - Add additional clustering metrics found in Scipy (agglomerative) are less likely to crash for gene-clustering as well as hopach using PIPE or Java equivalent. Add GUI options and the option for uses to visualize group associations with the cluster. Output results for TreeView.
  * **Comparison List Tool** - Allows for simple comparison of lists, recognizing known file types (e.g., alternative exon) and generates Venn Diagrams.
  * **Gene Correlation Tool Finder** - Input a gene and get the top correlated genes for GE or AS.
  * **View Gene, Pathway or Alternative Exon Expression Patterns** - Hierarchical clustering of alt.exons for all groups and samples or pathways. Also allow for gene-specific HC view (all exons or all junctions) or pathway view (all alt exons or diff. exp. genes).
  * **Lineage Profiler** - For significant cell types, display regulated genes associated with lineages as a network.
  * **SubgeneViewer** - Support for all splicing array platforms via SQLite relationship access. Assess feasibility of viewing as networks as opposed to Matplotlib.
  * **Robust RNA-Seq QC** - Identify and integrate more robust methods (FASTQC, EverSeq).

### Distribution ###
  * Evaluate virtualenv and buildout egg based distributions
  * For MatPlotLib, Numpy and Scipy, only make import calls to specific functions to decrease dependencies and distribution size.
  * Java wrapper of pure-python code in Cytoscape to allow for module-by-module migration from python to Java. Eliminates OS specific build issues. Requires using Cytoscape based visualization as an alternative (long term goal).

### Bug Fixes ###
  * Ontology annotation and nesting error within for non-GeneOntology categories (fixed in 2.0.8).
  * Excessive wait times for hierarchical clustering (reduced default allowed cluster size from 10k to 7k).
  * PAZAR GO-Elite annotation parser updated to the new PAZAR file format.
  * RNA-Seq counts filtering error for lines starting with 0 counts (resulted in inclusion of some gene fold changes that otherwise would have been excluded due to low expression).
  * For Affymetrix splicing-sensitive arrays, mitochondrial genome encoded probesets were missing from databases EnsMart65 and before due to Affymetrix and UCSC annotating these chromosomes as chrM rather than chrMT (Ensembl).
  * RPKM values were actually reported as a factor of 10 low. Hence, existing analyses will have had 10 times lower reported RPKM values than approaches (e.g., Cufflinks). This will not effect analyses except for RPKM level filtering which by default has been set very low (0.1 equivalent to a real RPKM equal to 1). The new RPKM filter is set to 0.5 for genes and exons (equivalent to 0.05 from AltAnalyze version 2.0.7).