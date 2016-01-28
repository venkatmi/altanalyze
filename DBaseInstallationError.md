### _Why Does AltAnalyze crash when trying to download the database when first run?_ ###

_**Answer:**_: The first time AltAnalyze is opened, the user will be prompted to allow AltAnalyze to establish an internet connection to download a species database. Since some locations have restricted proxy settings, an error may appear here, preventing an internet connection (note: this error is properly reported in version 1.16 of AltAnalyze, but still prevents proper download).

Although to properly fix this issue, http://www.altanalyze.org will need to be allowed proxy access, the species databases can also be downloaded manually using the following steps:

_Note: the following steps are for human exon array analysis. If this is not the case, see the "Alternative Species/Array Instruction" first._

#### Download and extract following zip files: ####
  1. Download the alternative exon database for your species [here](http://altanalyze.org/archiveDBs/AltDatabase/updated) (e.g., EnsMart55/Hs.zip for the latest human Exon 1.0 and Gene 1.0 database) to "AltAnalyze\_v2release/AltDatabase". If not analyzing human, mouse or rat, download the species archive with the extension "RNASeq.zip".
  1. Download the GO-Elite database for your species [here](http://altanalyze.org/archiveDBs/Databases/) (e.g., EnsMart55/Hs.zip for the latest human annotations) to the folder "AltAnalyze\_v1release/AltDatabase/EnsMart55/goelite" (EnsMart55 is used as an example).
  1. Download the GO-Elite OBO.zip file from [here](http://altanalyze.org/archiveDBs/Databases/) to the same goelite folder as in #2 and extract.
  1. Download the configuration files "species\_all.txt", "array\_versions.txt" and "source\_data.txt" from [here](http://www.altanalyze.org/archiveDBs/Config/) to the folder "Config". Although these files are already present, they are over-written as blank files when trying to connect online.

#### _OPTIONAL_ ####

#### Download and extract Affymetrix library files to "AltAnalyze\_v1release/affymetrix/LibraryFiles": ####
1) Download the .pgf, .clf and .bgp zip gzip files for your array from [here](http://altanalyze.org/archiveDBs/LibraryFiles). Human arrays have the prefix "Hu", mouse "Mo" and rat "Ra". Exon arrays .contain "Ex" in the filename, whereas gene arrays have "Gene".

Once completed, you should be ready to run AltAnalyze and will not be prompted to download any files. For assistance, ContactUs.