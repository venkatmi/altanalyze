## Creating an AltAnalyze Database From Scratch ##

Although most users will wish to wait for the official AltAnalyze database release, users and developers can also build and customize these databases themselves.

### Adding Support for a New Ensembl Build ###

New AltAnalyze databases can be built in a semi-automated fashion using builtin [command-line](CommandLineMode.md) arguments. Below is the exact method that AltAnalyze developers use to build each new version of the AltAnalyze database (excluding code updates). For microarray specific databases (e.g., Affymetrix Exon 1.0, Gene 1.0 and Junction JAY arrays), users will need to download two files for each array type and species before running the update. Files can be found [here](http://www.affymetrix.com/support/technical/byproduct.affx?).
  1. Probeset annotation file
  1. Probeset sequence file (Affymetrix Exon arrays only)

An example probeset annotation file is the "MoEx-1\_0-st-v1.na31.mm9.probeset.csv.zip". These files will need to be saved to the appropriate species and array type destination directory for that species and extracted to that directory. For example, if downloading the human exon array annotations for build 61 of Ensembl, save these to the folder AltDatabase/EnsMart61/Hs/exon.

Once installation is complete, additional files will have to be added to the final database installation package, however, these can be found with existing database installations (probes\_to\_remove and probeset-probes association - e.g., AltDatabase/EnsMart55/Hs).

**Updating Configuration File**

Default files to download are specified for UniProt, Ensembl, miRbase, TargetScan and the UCSC Genome Browser. For the custom AltMouse array, additional files are downloaded during the build process from our server if not locally found. For UniProt, UCSC, Ensembl databases, the user will not need to update the download path specified in the configuration file. When miRanda and TargetScan predictions are updated, you can change the path to the most recent version indicated on their websites (see the listed target URL) by opening the file Config/default-files.csv and changing the paths. Otherwise, no changes need to be made by the user.

**Build Commands**

_Build all databases for all platforms (RNASeq and arrays) and all species_
```
python AltAnalyze.py --species all --platform all --update all --version 61
```

_Build specific databases for select arrays and species_

```
AltAnalyze.py --update Ensembl --update UniProt --update Probeset --update Domain --update miRBs --species Mm --arraytype exon --arraytype gene --arraytype junction --version 61
```

_Build specific databases for RNASeq and species_

```
AltAnalyze.py --update Ensembl --update UniProt --update ExonAnnotations --update Domain --update miRBs --species Mm --platform RNASeq --version 61
```

_Note_: --arraytype and --platform are interchangeable as is Probeset and ExonAnnotations.

The machine these arguments are run on must have Python installed (recommended is version 2.7.1) and BioPython. To ensure BioPython is installed correctly, [follow these steps](setBioPython.md). These build commands, when run on a server will produce a complete AltAnalyze database with additional build files. To create a downloadable distribution archive for users, you can enter the below command to write species database zip files to the folder 'ArchiveDBs' in the AltAnalyze program directory.

```
python AltAnalyze.py --update package --version 61
```

### Adding Support for a New Species/Array/Data Type ###

Explicitly adding support for new species, array or data type (e.g., RNA-sequencing exon and junction data), requires either modification of the AltAnalyze source code or direct modification of the AltAnalyze database files.

## Under Construction ##