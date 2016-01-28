## Example Junction BED Format ##

Below is the column format used in the junction-level BED file recognized by AltAnalyze. When a .bed file containing these 12 fields is imported, AltAnalyze will assume the data are detected junctions from your experiment. This can be RNA-seq or even an unsupported microarray. Default expression is non-log read counts but can be log2 (select format as log in AltAnalyze). The notes field can include any or no information.

|chr|position1|position2|notes|expression|strand|start|stop|null|null|exon lengths|null|
|:--|:--------|:--------|:----|:---------|:-----|:----|:---|:---|:---|:-----------|:---|
|chr19|3266448  |3267268  |JUNC00000001|2         |-     |3266448|3267268|255,0,0|2   |11,27       |0,793|

_When the read is on the positive strand:_
  * 5' splice-site: position1+length1
  * 3' splice-site: position2-length2+1

_When the read is on the negative strand:_
  * 5' splice-site: position2-length2+1
  * 3' splice-site: position1+length1