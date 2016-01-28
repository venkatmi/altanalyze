## Introduction ##

The human hGlue array (http://gluegrant1.stanford.edu/~DIC/GGHarray/) is a next generation Affymetrix exon and exon-exon junction microarray that probe several hundred thousand thanscript features and over 20 thousand genes. Support for this arrays is introduced in AltAnalyze version 2.0.6. Gene database versions prior to EnsMart62 and executable versions of AltAnalyze prior to version 2.0.6 are not compatible with JAY array analyses.

Both exon and exon-junction analyses are executed in series on these arrays, produce separate result files for exon based alternative exon analyses (FIRMA, splicing-index, MiDAS) for single probesets and reciprocal exon-junction paired analyses (Linear Regression, ASPIRE). This allows the user to assess a change in exon inclusion based on combined junction and exon evidence. Methods used for reciprocal exon-junction analyses are the same as those used for the older generation [AltMouse](AltMouse.md) junction array ([Salomonis N. et al. PNAS 2010](http://www.pnas.org/content/107/23/10514)) and the [JAY](JAY.md) array.

## Details ##

Publications describing AltAnalyze use with the hGlue are arrays are anticipated in the near future, however, analyses with these arrays are also possible with other published methods (e.g., http://www.pnas.org/content/108/9/3707.full?sid=39868761-b5e0-42e8-a5d7-9cf113fad5ac).