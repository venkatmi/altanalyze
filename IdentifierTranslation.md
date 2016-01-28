# Identifier Translation in AltAnalyze #

A common use case for biologists dealing with genomics datasets is conversion of one identifier type to another. For example, analysis of AltAnalyze obtained Ensembl or microarray IDs in the ToppGene analysis suite (http://toppgene.cchmc.org) for downstream functional analyses currently requires that these IDs be converted to gene symbols or EntrezGene IDs. To accomplish this, users can access the Identifier Translation menu, load a file containing the IDs to be translated (must be the first column of values and obtain a new file in which the first column of values matches the desired ID type.

# Performing Identifier Translations #

From the AltAnalyze main menu, select your species and platform. Instead of indicating the platform you can also select "Other ID" as the data type and the gene system you are analyzing as the "platform ID". Select **Continue** and then **Additional Analyses**. From here, select **Identifier Translations** to see the options for translating your input file.

![http://altanalyze.org/image/IDtranslation.png](http://altanalyze.org/image/IDtranslation.png)
## Annotation Database ##

These translations are accomplished through use of relationships obtained from Ensembl and HMDB (GO-Elite database > AltDatabase/EnsMart65/goelite/Hs/uid-gene). All original IDs and other column data will be present in the output file, along with the Ensembl or HMDB IDs used for translation. Where multiple Ensembl or HMDB IDs are related to the input ID, only one will be chosen (last listed).