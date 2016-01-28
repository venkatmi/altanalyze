## Introduction ##

To run AltAnalyze from command-line, you will need to have installed AltAnalyze and make sure the source code is in the AltAnalyze program main directory. If you have downloaded the python-source code or Linux version, this will be the case, otherwise, you will need to copy the contents of the folder "Source\_code" to the parent AltAnalyze directory. Before supplying the command-line argument to this program, you will need to open a command prompt and change to the directory with the AltAnalyze source code.

Examples:

_**Downloading and installing a species specific database (human)**_

```
python AltAnalyze.py --species Hs --update Official --version EnsMart72
```


_**Exporting a Exon BED reference file for BedTools**_

```
python AltAnalyze.py --species Hs --platform RNASeq --bedDir "C:/BEDFiles" --output "C:/BEDFiles" --expname "hESC_differentiation" --buildExonExportFile yes
```


_**Analyzing  BED files – RNASeq Junction Bed files using default options**_

```
python AltAnalyze.py --species Hs --platform RNASeq --bedDir "C:/BEDFiles" --groupdir "C:/BEDFiles/groups.hESC_differentiation.txt" --compdir "C:/BEDFiles/comps.hESC_differentiation.txt" --output "C:/BEDFiles" --expname "hESC_differentiation"
```

_**Analyzing  BED files – RNASeq Junction Bed files using custom options**_

_Note: before running the below command, you will need to have created and saved the groups and comps files to the ExpressionInput directory. See [here](ManualGroupsCompsCreation.md) for instructions._

```
python AltAnalyze.py --species Hs --platform RNASeq --bedDir "C:/BEDFiles" --groupdir "C:/BEDFiles/ExpressionInput/groups.hESC_differentiation.txt" --compdir "C:/BEDFiles/ExpressionInput/comps.hESC_differentiation.txt" --output "C:/BEDFiles" --expname "hESC_differentiation"
```

For details on how to format generate and format the input text files, go [here](ObtainingRNASeqInputs.md).