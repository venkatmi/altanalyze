## _Running AltAnalyze from Command Line Interface_ ##

In addition to the graphical user interface (GUI), AltAnalyze can be easily run by command-line. This includes jobs run locallly, on a remote Linux server or cluster. This works fine given that the user knows the file paths of the directories containing input files, the output directory and has already created files containing the groups and comparisons for all samples analyzed.

**Creating Groups and Comparison Files** -
Creating groups and comparison files is needed beforehand, but is fairly easy. Just follow the directions listed [here](ManualGroupsCompsCreation.md). This can be done in an automated fashion as well, if input files have a [defined naming structure](PredictGroupsComps.md).

**Running Command Line from Source or Compiled Versions** -
The command-line can be run from the source code or OS-specific binaries. The binaries are recommended since these already contain graphical, statistical and webservice dependencies that need to be separately installed for the source code (see more information [here](StandAloneDependencies.md)).

When running with OS-specific binaries of AltAnalyze directly call the binary files themselves:
  * Windows OS _**` AltAnalyze.exe `**_
  * Mac OS X _**` AltAnalyze.app/Contents/MacOS/AltAnalyze `**_
  * Ubuntu OS _**` ./AltAnalyze `**_
  * Python source code _**` python AltAnalyze.py `**_

### Example Options ###

_**Downloading and installing a species specific database (mouse)**_

```
python AltAnalyze.py --species Mm --update Official --version EnsMart65 --additional all
```

_**Analyzing RNA-Seq files – BED format exons and junction results using default options and GO-Elite**_

```
python AltAnalyze.py --species Hs --platform RNASeq --bedDir "C:/BEDFiles" --groupdir "C:/BEDFiles/groups.YourExperiment.txt" --compdir "C:/BEDFiles/comps.YourExperiment.txt" --output "C:/BEDFiles" --expname "YourExperiment --runGOElite yes" --returnPathways all
```

_**Analyzing CEL files – Affymetrix 3’ array using default options and GO-Elite**_

```
python AltAnalyze.py --species Mm --platform "3'array" --celdir "C:/CELFiles" --groupdir "C:/CELFiles/groups.YourExperiment.txt" --compdir "C:/CELFiles/comps.YourExperiment.txt" --output "C:/CELFiles" --expname "YourExperiment" --runGOElite yes --returnPathways all
```

### Details ###

Many more additional example workflow analysis options and detailed option descriptions for various AltAnalyze functions are provided in the below links.

[Full AltAnalyze Workflows](CommandLinesWorkflows.md)

[Pathway Enrichment Analysis and Visualization](CommandLinesGOElite.md)

[Clustering, QC, and Alternative Exons Visualization](CommandLinesClustering.md)

[File comparison, ID translation and visualization](CommandLinesFileComps.md)

[LineageProfiler and Sample Classification](CommandLinesClassification.md)