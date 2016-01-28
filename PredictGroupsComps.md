### Introduction ###

Typically, groups and comparisons for array samples are assigned by the user in the graphical user interface (GUI) or built outside of AltAnalyze in a spreadsheet program. These files are used by AltAnalyze to (A) group samples according to experimental type and (B) establish pairs of groups to compare to each other for various statistical analyses (e.g., fold, p-value). When using the [command-line interface](CommandLineMode.md) of AltAnalyze, as opposed to the GUI, users have to specify the location of an groups and comparison file when beginning their run or place these with the [proper name](http://code.google.com/p/altanalyze/wiki/GroupsAndComps) in the ExpressionInput directory beforehand. In addition to these mechanisms, there is a way to have AltAnalyze automatically determine the groups and comparisons for your analyses. Because of the generic nature of this method, it is only recommended when a small set of groups are present in the analysis.

### Details ###

A specialized method was written in AltAnalyze named predictGroupsAndComps (UI.py), that predicts which Affymetrix CEL files correspond to which groups based on their filename and then builds a comparison file based on all possible pairwise group comparisons. The naming convention of these files is the group name first and the sample name second, delimited by a period, dash or an underscore.

|_**Proper Example Names**_|
|:-------------------------|
|tumor.sample1.CEL         |
|tumor-sample1.CEL         |
|tumor\_sample1.CEL        |
|Day 1 tumor.sample1.CEL   |

|_**Improper Example Names**_|
|:---------------------------|
|tumor-day1.sample1.CEL      |
|tumor.day1.sample1.CEL      |
|tumor\_day1.sample1.CEL     |

Only one delimiter type should be be present in single CEL file name. Note, that ".CEL" and ".cel" are ignored by this function, thus only reading one period delimiter when present elsewhere in the name. The groups and comparison files are created based on these file names when running AltAnalyze in the [command-line interface](CommandLineMode.md). If many groups are present, your analysis however can be very long since all possible pairwise comparisons are created. To avoid this, restrict your analysis to a few groups are analyze all groups as opposed to pairwise comparisons (flag --analyzeAllGroups " all groups").