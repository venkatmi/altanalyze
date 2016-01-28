## _Can AltAnalyze work with exon/junction data for an unsupported platform (e.g., custom array)?_ ##

_**Answer:**_ It is actually possible to analyze alternative exon expression in AltAnalyze 2.0 for any splicing sensitive expression data. To do this, you must treat your input expression data the same as RNA-Seq data, submitted in the UCSC BED format. The BED format files for exons and junctions differ slightly, with the exon BED format described [here](ExonBED.md) and the junction BED format [here](JunctionBED.md).

### Formatting your Expression Data ###

You will first need the genomic coordinates for all sequences probed on your platform (e.g., microarray). You will also want to match the expression values from your experiment to these coordinates. You will ultimately need to create a BED file for each sample (two if exon and junction data are present).

It is required that you know the genome assembly for that species to proceed, e.g., Human genome build hg19/GRCh37, and determine which genome builds correspond to which Ensembl builds. This can be determined by going to Ensembl genome information page (http://uswest.ensembl.org/Homo_sapiens/Info/Index) for your species and determining the genome build number (you may need to look at archived versions of this website to find the matching version). Since we always recommend matching your genome coordinates to the latest AltAnalyze supported Ensembl build, we recommend updating your BED files to the most recent genome build using UCSC's genome Lift Over website (http://genome.ucsc.edu/cgi-bin/hgLiftOver). Once your exon and/or junction BED files are properly formatted, you are ready to go (NOTE: errors during processing may be due to manual errors introduced during formatting).

### Running the analysis ###

Once you have created your BED files, you are ready to load this data into AltAnalyze and analyze this data as if it were RNA-seq data. See the [RNA-seq tutorial](Tutorial_AltExpression_RNASeq.md) for more information or [contact us](ContactUs.md).