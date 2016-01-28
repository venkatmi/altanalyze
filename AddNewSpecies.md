# Adding a New Species for Gene Expression Analysis #

Adding Affymetrix 3' array or gene array support for a new species couldn't be easier in AltAnalyze. The user just has to follow the three below easy steps:

  * **Step 1)** Open the file "Config/species.txt" in the AltAnalyze program directory and add the species two letter and long name (e.g., Dr and Danio rerio) to a new row. No need to add information to the third column (this column is reserved for arrays that allow for alternative splicing).

  * **Step 2)**  Open the file "Config/arrays.txt" and in the last column, with species two letter codes, add your species code to the existing list, separated by a pipe (e.g., Mm|Hs|Rn|Dr to Mm|Hs|Rn) to the row with the ID 3'array.

  * **Step 3)** Open the folder AltDatabase and affymetrix. Here add a new folder with the two letter species name (e.g., "AltDatabase/affymetrix/Dr"). To this save the Affymetrix CSV annotation file provided at the Affymetrix website (requires login). To see compatible arrays and select the CSV file for download, go to: http://www.affymetrix.com/support/technical/byproduct.affx?.

Once these steps are completed are you ready analyze your conventional microrarray data! By saving the CSV annotation for your array into the AltDatabase/affymetrix directory, AltAnalyze will be able to propperly annotate the genes and probesets on your array with the most up-to-date gene IDs and information, inlcuding the Gene Ontology and pathways.

Note: these additions only allow for the analysis of 3' Affymetrix and Gene 1.0 arrays and do not allow for more complex splicing analyses for a new array in AltAnalyze.  Adding such support can be complex, usually involving the help of AltAnalyze developers.