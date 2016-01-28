# Download Instructions for FASTQ Data #

As input to AltAnalyze, you can currently upload exon and/or junction coordinate and read count data obtained from various RNA-seq alignment tools (e.g., TopHat, BioScope, HMMSplicer, SpliceMap). If you have your own sequence data, we recommend starting with this rather than downloading demo data, however, several sample datasets are listed below. You can download one of these studies from [SRA](http://www.ncbi.nlm.nih.gov/sra):

  * ERP000570 GSE21860 GSE27218 (Mus musculus)
  * GSE24447 (Homo sapiens)

_**Note: SRA is in the process of being discontinued, hence, FASTQ files may be unavailable.**_

These studies contain multiple replicates for multiple in vivo/in vitro conditions. To obtain FASTQ sequence:
  1. Search for the above accession numbers at SRA.
  1. For each sample, select the experiment accession number on the right side of the page (under Total:).
  1. Select the button "Download" half way down the page.
  1. On the next page select "FASTQ", de-select "filtered" and select "save".
  1. Repeat for each sample (often 1-7GB files) and rename after download.

These files will need to be moved to a suitable machine for (ideally multiple processors with 8GB plus RAM) for sequence alignment. For paired-end data, both paired ends are stored sequentially in a file, for example:

@SRR037165.2.1 FC30B11:1:1:256:1595.1 length=35

@SRR037165.2.2 FC30B11:1:1:256:1595.2 length=35

You can either directly load this data in a sequence alignment program or divide the file into two files with each paired-end stored in a different file using other software. If multiple read files for a single sample are downloaded, you can combine these using the Linux command **cat**.