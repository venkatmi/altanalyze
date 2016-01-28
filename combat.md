## Batch Effect Correction in AltAnalyze Using Combat ##

Batch effect removal is provided using the combat python package, which provides nearly identical results to that of the R version of combat (https://github.com/brentp/combat.py). When applying combat in AltAnalyze, the user must indicate which samples correspond to which batches. This is done in a batch annotation window within the GUI or can be manually created in Excel with the same format as the groups file, containing the prefix batch.

When AltAnalyze calls combat, it will produced an additional file with the prefix pheno. which it uses as input for combat (see http://www.bioconductor.org/packages/2.12/bioc/html/sva.html).