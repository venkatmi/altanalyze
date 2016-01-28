## Introduction ##

The mouse and human JAY arrays ([MJAY](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8940) and [HJAY](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8444)) are next generation Affymetrix exon and exon-junction microarrays that probe several hundred thousand thanscript features and over 20 thousand genes. Support for these arrays is introduced in AltAnalyze version 1.16. Gene database versions prior to EnsMart55 and executable versions of AltAnalyze prior to version 1.16 are not compatible with JAY array analyses.

Both exon and exon-junction analyses are executed in series on these arrays, produce separate result files for exon based alternative exon analyses (FIRMA, splicing-index, MiDAS) for single probesets and reciprocal exon-junction paired analyses (Linear Regression, ASPIRE). This allows the user to assess a change in exon inclusion based on combined junction and exon evidence. Methods used for reciprocal exon-junction analyses are the same as those used for the older generation [AltMouse](AltMouse.md) junction array ([Salomonis N. et al. PNAS 2010](http://www.pnas.org/content/107/23/10514)).

## Details ##

Publications describing AltAnalyze use with the MJAY and HJAY are arrays are anticipated in the near future, however, analyses with these arrays are also possible with other published methods (e.g., http://bioinformatics.oxfordjournals.org/content/26/2/268.abstract)