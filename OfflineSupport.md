## Offline Support for AltAnalyze ##

When starting AltAnalyze 2.0 and greater, it will require the download of multiple zipped database files. If you do not have an internet connection on that machine or AltAnalyze is unable to download the necessary files due to a proxy, firewall or other issue, follow the below directions.

On a machine with support for downloading files through http, download the following files. This example assumes the user is working with human or Hs data.

### Base Download Files ###

**Primary Gene Database**

http://altanalyze.org/archiveDBs/AltDatabase/updated/EnsMart65/Hs.zip (save to AltDatabase folder and extract)

**Lineage Profiler Database**

http://altanalyze.org/archiveDBs/AltDatabase/updated/EnsMart65/Hs_LineageProfiler.zip (save to AltDatabase/EnsMart65/ensembl and extract)

**GO-Elite Annotation Database**
http://altanalyze.org/archiveDBs/Databases/EnsMart65Plus/Hs.zip (save to AltDatabase/EnsMart65/goelite and extract

**GO-Elite Ontology Database**

http://altanalyze.org/archiveDBs/Databases/EnsMart65Plus/OBO.zip (save to AltDatabase/EnsMart65/goelite and extract)


### Support for Junction Arrays and RNA-Seq ###

If you want to perform a Affymetrix junction array or an RNA-Seq analysis, additionally download the following to the folder AltDatabase and extract:

**HJAY array**

http://altanalyze.org/archiveDBs/AltDatabase/updated/EnsMart65/Hs_junction.zip

**HGlue array**

http://altanalyze.org/archiveDBs/AltDatabase/updated/EnsMart65/Hs_junction_Glue.zip (Glue array)

**RNA-Seq**

http://altanalyze.org/archiveDBs/AltDatabase/updated/EnsMart65/Hs_RNASeq.zip

### Microarray Library Files ###

Affymetrix arrays require library files to be downloaded. These can be obtained from the manufacturer or obtained from the below directory:
http://www.altanalyze.org/archiveDBs/LibraryFiles/

Download ALL of the files for your array to AltAnalyze/AltDatabase/affymetrix/LibraryFiles and extract.

### Other dependencies ###

When running the python source-code directly, other Python specific dependencies will be need to be installed to support data visualization, clustering and QC analyses. See
[here](http://code.google.com/p/altanalyze/wiki/StandAloneDependencies) for more details.