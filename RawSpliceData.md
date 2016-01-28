The RawSpliceData directory contains alternative exon scores or normalized intensity values for each sample analyzed in a study. These values are typically log2 values (e.g., the number one is equivalent to a fold of 2). Since scores are calculated by comparing two specified biological groups, the geometric difference of the mean of the values for the two reported groups should be equal to the reported score.

## Exon Arrays ##
For exon arrays, normalized probeset intensities are reported. These are the log2 values of the probeset sample expression divided by the calculated constitutive expression (aka gene expression) for that gene. Constitutive expression is calculated according to the [following methods](GeneExpCalculation.md).

## Junction Arrays ##
For junction arrays, the ratio of normalized intensities for the two compared probesets (e.g., reciprocal exon junctions) for each sample (inclusion probeset/exclusion probeset).