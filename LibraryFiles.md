## Introduction ##

The directory LibraryFiles is an AltAnalyze specific folder contained with the AltAnalyze program directory "AltAnalyze\_v1release/affymetrix/LibraryFiles". The folder is where AltAnalyze stores Affymetrix library files (e.g., CDF, PGF, CLF and BGP). AltAnalyze will automatically download the appropriate files to this directory when analyzing an microarray that requires these library files for expression summarization. Expression summarization is performed using the algorithm [RMA](RMA.md) by the program [APT](APT.md). The files saved here by AltAnalyze or by the user will be used for all subsequent analyses.


## Details ##

AltAnalyze automatically downloads library files from the following archive site:
http://altanalyze.org/archiveDBs/LibraryFiles/

or from http://www.affymetrix.com