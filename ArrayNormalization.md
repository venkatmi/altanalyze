## Introduction ##
Before using AltAnalyze, the user must have CEL files from their experiment or have processed their Affymetrix microarray CEL files to produce a normalized and background corrected expression value file. In addition to this file, it is recommended that the user produce a detection p-value file for alternative splicing analyses. It is recommended that you use AltAnalyze to perform RMA and obtain detection p-values on your CEL files, but other external options the has are also provided in this tutorial.

## Using Your Own Data ##
Save all of the Affymetrix .CEL files from your study to a new folder on your computer. Proceed to one of the below options.

## Using Sample Data ##
  * (**Exon array analysis**) If you would like to use sample exon array data, you can download CEL files posted at NCBI's Gene Expression Omnibus from [here](ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/GSE13297/GSE13297_RAW.tar). This study is a comparison of human embryonic stem cells and differentiated cardiomyocytes.This will save a compressed file (GSE13297\_RAW.tar - 147MB) to your hard-drive.
  * (**Non-exon array analysis**) If will only be doing gene expression analyses with conventional arrays you can download a mouse sample dataset [here](http://www.genmapp.org/AltAnalyze/mESC_timecourse_A.tar).
  * Move this file to a new folder and extract the compressed files within this .tar file to this folder. For PCs, we recommend the free program [WinRAR](http://winrar.d0wnloadz.net/).
  * The extacted files themselves will be compressed using a gzip algorithm. Select all of these files, right-click and choose to extract to the current location.

## Option 1) AltAnalyze (Cross-Platform) ##
Using AltAnalyze to directly process your CEL files is the simplest option for most users. Instructions for using AltAnalyze to perform RMA analysis to obtain expression and detection p-value files, can be found here for gene-expression only and here for exon arrays. Before processing your CEL files, you will need to have the relevant library files from the specific microarray as well as the appropriate NetAffx CSV annotaiton file. Most mouse, human and rat arrays can be automatically downloade and installed by AltAnalyze but otherwise will need to be downloaded from the Affymetrix website. CEL files are processed using modules from Affymetrix Power Tools (APT), distributed with AltAnalyze. Further instructions are provided within AltAnalyze.

Option 2) Affymetrix Expression Console (Windows Only)
An alternative to processing CEL files in AltAnalyze is the program Affymetrix Expression Console. Users using this program have additional options for array annotation and summarization (e.g. PLIER summarization). This program supports multiple microarray platforms, normalization algorithms and is able to download necessary library files. A pictorial view of these directions, including most of these steps is available [here](http://domaingraph.bioinf.mpi-inf.mpg.de/docu/affy_pre.php). To begin:

  * Download Expression Console from [here](http://www.affymetrix.com/products_services/software/specific/expression_console_software.affx). Prior to downloading this software, you must sign up for a free Affymetrix account.
  * Extract the compressed download to a temporary new folder. Open the file with the extension ".exe" to install the software.
  * Go to the Windows Start Menu>All Programs>Affymetrix>Expression Console
  * Expression Console will prompt to you enter a new user profile name. Do so.
  * Select the Menu option File>New Study.
  * Select the "Download library files..." option from the File menu and enter your Affymetrix website username and password. Select the array type you will be analyzing data from (e.g., HuEx-1\_0-st-v2 or MOE430A). Expression Console will download the necessary library files.
  * After the download is complete, in the new window, select Add Intensity Files. Select the folder that contains your .CEL files (note: if you are looking for the Desktop, you'll need to go to C:\Documents and Settings\<your user name>\Desktop) and shift-select all .CEL files you wish to analyze.
  * Select the button "Run Analysis" and under "Available Analyses" select the RMA analysis for your array (recommended). For an Exon 1.0 ST array, select "Exon Level - All: RMA-Sketch (including speculative content)" (recommended - full is also alright, but excludes some content).  For all other arrays, choose RMA analysis option.
  * When prompted to add a suffix name for the study type any name or no name.
  * Expression Console will now run. This may take 30 minutes or longer, depending on the number and type of arrays. If analyzing exon array data, when finished, two new files will be saved to the folder containing your CEL files named "rma-exon-all.summary.txt' (expression file) and "dabg.summary.txt" (detection p-value file), otherwise, only an expression file will be saved. If you like you can also perform QC on your data by checking out the QC options on the lower right hand part of the Expression Console screen.
  * (**Exon array analyses**) Re-name these two files from "rma-exon-all.summary.txt' to "exp.hESC\_differentiation.txt" and "dabg.summary.txt" to "stats.hESC\_differentiation.txt" or whatever you would like to call your experiment.
  * (**Non-exon array analyses**) Re-name the file "rma.summary.txt" to "exp.mESC\_differeniation.txt".
  * Save these files into the appropriate AltAnalyze\_v1beta folder directories ("AltAnalyze\_v1beta/ExpressionInput/exon" or AltAnalyze\_v1beta/ExpressionInput/3'array").
  * Run AltAnalyze.

## Option 3) Affymetrix Power Tools (Cross-Platform) ##
Affymetrix Power Tools (APT) provides similar analysis options to those found in Expression Console, but is not restricted to Windows. Unlike Expression Console, APT has a command line interface as opposed to a graphical user interface (GUI). Two modules for APT are bundled with AltAnalyze that allow users to easily process CEL files using the APT via the AltAnalyze GUI as well as calculate MiDAS p-values for pair-wise comparisons (see AltAnalyze ReadMe or tutorial 2).

### Installing APT ###
Download APT for your operating system at: http://www.affymetrix.com/partners_programs/programs/developer/tools/powertools.affx#1_2. Unzip the program archive and install APT from the extracted installer.

**Additional Files**

  * Create a new text file with a list of Cel files, named “celfiles.txt” and save to the directory containing your CEL files. This file should look like this:

> cel\_files
> GSM335817.CEL
> GSM335818.CEL

### Running APT ###

**Exon Arrays**

Go to the Affymetrix website and download the corresponding Library files for your array. For example, the for the Human Exon 1.0 ST Array, download the non-CDF library files zip archive (Human Exon 1.0 ST Array Analysis). Extract this archive to a new folder that also contains your CEL files to analyze.

> _Windows_

  * Open APT Command Prompt. On a PC, this can be found from your Start menu under All Programs>Affymetrix Power Tools>apt-"version">APT Command Prompt.
  * Change directories to the one containing your CEL files. For example, if the files are saved in a folder named “Affy” in your root C drive, type in: “cd C:\” and then “cd Affy”.
  * Type in the following APT command to have APT create an expression and detection p-value files, in this example, for the human exon array:

apt-probeset-summarize  -p HuEx-1\_0-st-v2.[r2](https://code.google.com/p/altanalyze/source/detail?r=2).pgf  -c HuEx-1\_0-st-v2.[r2](https://code.google.com/p/altanalyze/source/detail?r=2).clf  -b HuEx-1\_0-st-v2.[r2](https://code.google.com/p/altanalyze/source/detail?r=2).antigenomic.bgp  --qc-probesets HuEx-1\_0-st-v2.[r2](https://code.google.com/p/altanalyze/source/detail?r=2).qcc  -s HuEx-1\_0-st-v2.[r2](https://code.google.com/p/altanalyze/source/detail?r=2).dt1.hg18.full.ps  -a rma-sketch -a dabg -o output-gene  --cel-files celfiles.txt

**Gene Arrays**

Go to the Affymetrix website and download the Library files for your array. For example, the for the Human Gene 1.0 ST Array, download the non-CDF library file (Human Exon 1.0 ST Array Analysis).

> _Windows_

  * Change directories to the one containing your CEL files. For example, if the files are saved in a folder named “Affy” in your root C drive, type in: “cd C:\” and then “cd Affy”.
  * Type in the following APT command to have APT create an expression file, in this example, for the human exon array:

apt-probeset-summarize -p HuGene-1\_0-st-v1.[r3](https://code.google.com/p/altanalyze/source/detail?r=3).pgf -c HuGene-1\_0-st-v1.[r3](https://code.google.com/p/altanalyze/source/detail?r=3).clf -b HuGene-1\_0-st-v1.[r3](https://code.google.com/p/altanalyze/source/detail?r=3).bgp --qc-probesets HuGene-1\_0-st-v1.[r3](https://code.google.com/p/altanalyze/source/detail?r=3).qcc -m HuGene-1\_0-st-v1.[r3](https://code.google.com/p/altanalyze/source/detail?r=3).mps -a rma-sketch -o output-gene --cel-files celfiles.txt

**3' Arrays**

Go to the Affymetrix website and download the Library files for your array. For example, the for the Mouse 430 A Array, download the library file zip archive. This contains several directories. You will need to move the file named "MOE430A.CDF" from the folder "CD\_MOE430/Full/MOE430A/LibFiles" to the directory containing your CEL files.

> _Windows_

  * Change directories to the one containing your CEL files. For example, if the files are saved in a folder named “Affy” in your root C drive, type in: “cd C:\” and then “cd Affy”.
  * Type in the following APT command to have APT create an expression file, in this example, for the human exon array:

apt-probeset-summarize -d MOE430A.CDF -a rma -o output-dir cel-files celfiles.txt

In the folder named “output-gene”, there will be an expression file (e.g., rma-sketch.summary.txt) and optionally, a detection p-value file (e.g., dabg.summary.txt). Re-name these two files from " rma-sketch.summary.txt” to "exp.hESC\_differentiation.txt" and from “dabg.summary.txt” to "stats.hESC\_differentiation.txt" or whatever you would like to call your experiment.

## Option 4) easyExon (Cross-Platform - requires Java SE 6) ##

The program easyExon is a Java based application that you can run over the web to process Affymetrix Exon 1.0 arrays (but not other platforms), perform down-stream statistics and visualize probe set data at the gene-level. This program requires pre-installation of APT along with the supporting array annotation files (e.g., PGF, BGP) to process Affymetrix CEL files. Since it requires Java SE 6, operating systems, such as Mac OS X 10.4 (Tiger) are not compatible with easyExon at this time, since this version of Java SE is not supported.