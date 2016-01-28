# Introduction #

RMA or [Robust MultiArray Average](http://www.ncbi.nlm.nih.gov/pubmed/12925520) is a method commonly used for microarray summarization and normalization. RMA is the default method using by AltAnalyze for Affymetrix microarray expression normalization. To do this, AltAnalyze calls an imbedded version of the Affymetrix software [APT](APT.md). The result of RMA are sample expression values with a near identical distribution of expression values. APT also produces detection above background (DABG) p-values, for all Affymetrix arrays, to determine if a particular probeset is reasonably expressed (e.g. Present or Absent).

# Details #
http://www.plexdb.org/modules/documentation/RMAexplained.pdf