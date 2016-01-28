# Feature Requests #
### Version 1.12 (release date - June 2009) ###

  * Directly process CEL files and run MiDAS via APT remotely (no user command line interface required - cross-platform).
  * Add an a new interface for creating groups and comparison files using just an expression dataset.
  * User-defined input and results directories.
  * Export a GenMAPP/PathVisio import text file.
  * Include direct domain aligningment for probesets annotations in the primary alternative exon results file (similar to DomainGraph).
  * Export a new DomainGraph input file containing analyzed and not-regulated probesets with all primary splicing statistics.
  * Select files and folders for analysis using the browse option.
  * Better error handling and information on errors (or re-routing to previous menus).
  * Warn user if no related Affymetrix CSV annotations are present and instruct how to download.
  * Create a log file for each run for users to view possible issues (saved to user-defined result folder).
  * More sophisticated protein and domain level prediction with the following options:
    * Multiple transcript isoform comparison - compare all possible gene transcripts at the level of exon composition to find minimal differences, for probe aligning and not aligning transcripts in public databases (using genomic coordinates). Store protein and domain associations as a table for input rather than generate on the fly.
    * Multiple protein isoform comparison - similar as above, but instead of finding transcripts with minimal differences in exon structure, identify minimal protein composition differences between all possible isoforms (N-terminal/C-terminal differences, domain changes, protein length) to find a pair with the smallest impact on protein sequence. Store protein and domain associations as a table for input (note: please contact us if you want to use this function - files not inlcuded with standard download)
    * Direct domain alignment - Analyze all gene linked proteins and associated domains by matching domain and probeset  genomic coordinates, similiar to DomainGraph. This method is less biased and is better at identifying predicted domain changes for uncharacterized isoforms. However, it does not include UniProt sequence regions with functional annotations, just Ensembl InterPro (e-score<1). Store protein and domain associations as a table for input.
  * Automatically download and organize all Ensembl source files from SQL files on the Ensembl website (currently requires BioMart and custom Perl script accessing Ensembl API).
  * Add custom annotation layers to the gene-expression and alternative gene files - protein class (TF, GPCR) and protein compartment (membrane, nucleus).
  * Introduction of a command line interface allowing AltAnalyze to be run remotely or locally through a terminal or from other local applications.
  * GO-Elite functionality option for gene expression and alternative exon analysis. Requires a user input option gene-expression analyses to define cut-offs to establish criterion.

### Version 1.13 ###

  * Improved GUI - Handles version specific database updates and selection of different array types and array manufacturers.
  * Complete user controlled update - Specific database versions, based on Ensembl builds will be introduced, including previously provided species databases. This includes download of GO-Elite annotatation and relationship species databases.
  * Command-line run options - Rather than run the program through the UI, all primary functionality is available through command line flags passed to the python source code, allowing remote access to AltAnalyze by command-line. This includes database all major update features.
  * Improved support for non-Affymetrix arrays - Comprehensive annotation of Agilent, Illumina, Codelink and other arrays supported by the Ensembl database are now supported.
  * Fixed Linux related bugs - AltAnalyze should be able to read and write files propperly to different Linux operating systems.
  * Affymetrix CSV annotation files no longer needed - Although these annotations can still be used in the same way as before, these will also automatically be extracted from the local GO-Elite databases downloaded with each species installation.
  * Version specific exports to DomainGraph - Output files contain flags to which version of Ensembl was used, ensuring the same version is used in DomainGraph.

### Version 1.2 ###

  * Introduce generic methods for linking genomic coordinates of sequences and probes to known exons, introns, UTRs and junctions.
  * Optional alternative splicing gene expression calculation using either custom or all probe sets.
  * Support for multiple group analyses rather than pairwise comparison of two biological groups.
  * Support for RNAseq data as well as other exon, junction or combination arrays.
  * Introduction of novel methods for alternative splicing domain and motif analysis.

### Version 2.0 ###

  * Option for webstart and automatic loading of data in DomainGraph.
  * Table browser for results with linkout to relevant resources (probeset, gene, domain, microRNA binding site).
  * Progress bar for update process.
  * Graph display of expression changes for all exon or exon/intron probe sets for a gene.
  * Low-level probe set analyses and splicing-index probe level calculation (similiar to or incorporating JETTA integration via R).
  * Additional sophisticated protein and domain level prediction as options with the existing:
    * De novo isoform prediction - Using the longest characterized isoform aligning to the probeset, identify the most likely aligning exon sequence (containing the probeset sequence, but accounting for alternative 3' and 5' AS) and predict the corresponding mRNA sequence of the isoform missing this exon. Allows you to identify truly novel protein predictions that should have minimal impact on the reference mRNA sequence. Allows for analysis of novel splice events and the ability to assess downstream domain-level change that would not directly overlap with the probeset. Store associations as a table for input.
    * Combination - Use all of the above methods (version 1.1 and 2.0) and choose the ones that are the most reliable prediction methods (most conservative).  Otherwise, let the user choose which method to use.