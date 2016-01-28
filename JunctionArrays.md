## Affymetrix Junction Arrays ##

AltAnalyze is capable of immediately analyzing Affymetrix files from several Affymetrix junction array platforms. These microarrays probe both known and predicted exon-exon junctions as well as known and predicted exon regions. Expression analysis is conducted similar to exon array analyses by processing the Affymetrix CEL files using RMA (via an automatic call to [APT](APT.md)) and examining probeset-level versus gene-level differences for the different array features. Currently, three array platforms are supported, the [JAY](JAY.md) array, [hGlue](hGlue.md) array and [AltMouse](AltMouse.md) array. Details on the analysis methods used for these platforms are available in the associated platform links. For users analyzing a different junction array, see the [following link](SplicingForAnyPlatform.md).