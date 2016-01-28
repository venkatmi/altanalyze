### _If two probesets align to one exon and have opposite splicing changes, what does this mean?_ ###

_**Answer**_: Such a result may be meaningful. For example:

  * Two adjacent probesets may associate with an alternative 5' or 3' splice site, where one isoform is up and the other is down. In AltAnalyze, if this is a characterized splice event, one probeset will have the annoation "FivePrime?" or "ThreePrime?". In DomainGraph, you can also see if the probesets have different distributions in different transcripts. This can still be alternative splicing however, even if the annotation is missing if this is an uncharacterized event or in an EST not analyzed by AltAnalyze. We also recommend you check http://www.genome.ucsc.edu for that probeset.

  * Alternatively, this could be due to a hybridization issue, where a probeset targets multiple locations in the genome. You can BLAT the probeset sequence (obtain from NetAffx?) at http://www.genome.ucsc.edu to verify.