# Building a Dataset Exon Database for BEDTools #

While junction-level input files can be directly produced by several junction finding algorithms (e.g., TopHat, HMMSplicer, SpliceMap) for AltAnalyze, exon-level input files must currently be produced through an intermediate set of steps. These steps are:
  1. Obtain a BAM alignment file (e.g., TopHat)
  1. Output known and novel exons from junction.bed files with AltAnalyze
  1. Produce exon-level bed read counts with BEDTools
  1. Name the exon and junction BED files so that AltAnalyze matches them up

## Instructions ##

We will illustrate this pipeline for use with TopHat.
  * Install and run TopHat using the appropriate parameters as described [here](TopHat.md).
  * Once finished, rename the TopHat junction BED files (e.g., sample1) and move all sample files to a single directory.
  * Download and extract the latest version of AltAnalyze (http://www.altanalyze.org). If using by [command-line](RNASeqCommandLine.md), download the source-code version and ensure Python is installed (from command-line, type "python").
  * Download and install BEDTools as described [here](BEDTools.md).
  * At this point you can choose to use AltAnalyze via the graphical user interface or command-line.

### AltAnalyze Graphical User Interface Instructions ###

  * When prompted, choose to connect to the internet and choose your species database of interest.
  * When finished installing, ensure the correct species and RNA-seq dataset option are chosen and select "Continue".
  * If mouse, human or rat, you will be prompted to download additional files. Select the download option and continue.
  * Select the option "Process RNA-seq reads" and "Continue".
  * Enter a descriptive dataset name (e.g., hESC\_differentiation), select the folder containing your junction.bed files and select the same folder as an output directory.
  * At the bottom of this window, select the option "Build exon coordinate bed file to obtain BAM file exon counts". This will instruct AltAnalyze to select novel junction coordinates that correspond to novel exons to add to known Ensembl/UCSC exon regions.
  * A new window will appear, confirming that you will find the AltAnalyze generated exon database in the BAMtoBED folder of your output directory. Select "Continue".
  * A status window will appear indicating the progress of the junction.bed file import and annotation. Once finished, a final indicator window will appear with the path of the exon database (BED file) indicated.
  * A new file will be found in your output directory under the folder BAMtoBED (e.g., Hs\_hESC\_differentiation\_EnsMart62\_exons.bed).

### AltAnalyze Command-line Instructions ###

  * Download the latest species database as demonstrated [here](RNASeqCommandLine.md) under "Downloading and installing a species specific database".
  * Once installed, construct a new set of command-line arguments to select the folder containing your junction.bed files and output the BAMtoBED output as described [here](RNASeqCommandLine.md) under "Build a BAMtoBED exon database".
  * A new file will be found in your output directory under the folder BAMtoBED (e.g., Hs\_hESC\_differentiation\_EnsMart62\_exons.bed).

### Run BEDTools ###

Once BEDTools is installed, the "bamToBed" utility should be callable from the command-line. If not, you may need to edit your .bashrc file. With this utility, you will instruct BEDTools to extract read counts for exon coordinates in your BAMtoBED exon database from AltAnalyze. Example instructions are shown [here](BEDTools.md).

### Re-run AltAnalyze ###

  * Once the exon-level BED files have been produced, save these to the same folder as your junction.bed files. The junction and exon files will need to have the same beginning name, separated from a unique name by a double underscore. For example:

|`Cancer_s1__canonical-junction.bed`|
|:----------------------------------|
|`Cancer_s1__noncanonical-junction.bed`|
|`Cancer_s1__exon.bed`              |
|`Wt_s1__canonical-junction.bed`    |
|`Wt_s1__noncanonical-junction.bed` |
|`Wt_s1__exon.bed`                  |

  * Multiple junction or exon files for a sample can be included as AltAnalyze will combine these. AltAnalye will differentiate between exon and junction BED files based on the number of columns automatically present (see ExonBED and JunctionBED).
  * Run AltAnalyze again, as shown above, without the BAMtoBED option to obtain gene expression, alternative exon and over-representation statistics.
  * Additional instructions and tips can be found [here](Tutorial_AltExpression_RNASeq.md).