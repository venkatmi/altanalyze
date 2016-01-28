## Example Exon BED Format ##

Below is the column format used in the exon-level BED file recognized by AltAnalyze. When a .bed file containing these 10 fields is imported, AltAnalyze will assume the data are detected exons from your experiment. This can be RNA-seq or even an unsupported microarray. Default expression is non-log read counts but can be log2 (select format as log in AltAnalyze). The notes field can include any or no information.

|chr|start|stop|notes|null|strand|expression|null|region-size|null|
|:--|:----|:---|:----|:---|:-----|:---------|:---|:----------|:---|
|chr15|76527163|76527449|ENSMUSE00000393034|0   |+     |22        |0   |286        |0   |