# Creating Groups and Comps Outside AltAnalyze #
## Creating a Comparison and Groups File ##
Two other files, besides expression are also required by AltAnalyze to instruct the program which samples correspond to which groups (groups file) and which groups you want to compare (comparison file). While these files can be produced within AltAnalyze through the user interface by the user, for very large datasets with hundreds of samples and dozens of groups or users running AltAnalyze by command-line, the user may need to do this outside of AltAnalyze in Excel.  When made by the AltAnalyze GUI, the program will save two files with the prefix "groups." and "comps" to the same folder with the expression file (must have the same suffix as the expression file). These files are saved to the same directory as the expression file (ExpressionInput).

The below example is with Affymetrix CEL files, but the same format applies to using .bed or .tab files (RNASeq). To make the groups file outside of AltAnalyze:

## Starting from CEL Files: ##
  * Open AltAnalyze, proceed to analyze CEL files and decide on a name for your dataset. Remember this name or copy it, as it will be needed later.
  * Proceed to the next window and QUIT AltAnalyze.
  * Open the file "cel\_files.txt" in the folder containing your CEL files in Microsoft Excel, OpenOffice or Google Docs. With these CEL files names proceed to Setting Up the "Groups File section below".

## Starting from a Processed Expression File: ##
  * Open the expression file in a spreadsheet application such as Microsoft Excel, OpenOffice or Google Docs. The program should warn you that the file is not loaded completely. This is OK, just make sure not to save this file.
  * Once opened, copy the row containing names of the CEL files. Note: this column may not be the first column in the file, depending on the program used to summarize the expression values. Scroll down until you see a row with column headers above the expression values. The first column should say "probeset\_id".
  * Paste transpose the column headers into a new spreadsheet. In Microsoft Excel, this is through the menu option Edit>Paste Special>Transpose. Make sure only the .CEL files names appear (not the probe set column header). With these CEL files names proceed to Setting Up the "Groups File section below".

## Setting Up the Groups File ##
| GSM86112-0hr1.CEL    | 1 | 0hr    |
|:---------------------|:--|:-------|
| GSM86114-0hr2.CEL    | 1 | 0hr    |
| GSM86116-0hr3.CEL    | 1 | 0hr    |
| GSM86118-6hr1.CEL    | 2 | 6hr    |
| GSM86120-6hr2.CEL    | 2 | 6hr    |
| GSM86122-6hr3.CEL    | 2 | 6hr    |
| GSM86124-12hr1.CEL   | 3 | 12hr   |
| GSM86126-12hr2.CEL   | 3 | 12hr   |
| GSM86128-12hr3.CEL   | 3 | 12hr   |
| GSM86136-4days1.CEL  | 4 | 4days  |
| GSM86138-4days2.CEL  | 4 | 4days  |
| GSM86140-4days3.CEL  | 4 | 4days  |
| GSM86142-7days1.CEL  | 5 | 7days  |
| GSM86144-7days2.CEL  | 5 | 7days  |
| GSM86146-7days3.CEL  | 5 | 7days  |
| GSM86154-14days1.CEL | 6 | 14days |
| GSM86156-14days2.CEL | 6 | 14days |
| GSM86158-14days3.CEL | 6 | 14days |
| GSM86300-24hr1.CEL   | 7 | 24hr   |
| GSM86302-24hr2.CEL   | 7 | 24hr   |
| GSM86304-24hr3.CEL   | 7 | 24hr   |

  * If the name of the expression file was "exp.mESC\_differentiation.txt" or "mESC\_differentiation.txt", this file should be saved as "groups.mESC\_differentiation.txt" or "groups." plus the name of the experiment, to the same folder as the expression file.
  * The comparison file is simply two columns, each with a single number. One row illustrates comparison between the groups designated in the groups file.  In this case, we wish to compare different time-points of differentiation to 0hr, where 0hr is the control group.

> | 2 | 1 |
|:--|:--|
> | 3 | 1 |
> | 4 | 1 |
> | 5 | 1 |
> | 6 | 1 |
> | 7 | 1 |

  * This file should be saved as "comps.mESC\_differentiation.txt", to the same folder.