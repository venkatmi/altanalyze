The ExpressionInput directory contains files used for all downstream AltAnalyze expression analyses.

When beginning with Affymetrix CEL files, normalized expression values for all probesets are stored to this directory (exp.YourExperiment.txt) that result from APT ([RMA](RMA.md)).

For splicing-sensitive arrays, a file containing detection probability statistics ([DABG](DABG.md)) for all probesets is also saved (stats.YourExperiment.txt).

When groups and comparisons are designated for your experimental samples, the files groups.YourExperiment.txt and comps.YourExperiment.txt are also saved here.

In addition to these files, interim files created by AltAnalyze are also saved here and used by the program. Currently, these include probeset residuals (residuals.YourExperiment.txt) and constitutive gene summarized expression values (exp.YourExperiment-steady-state.txt). For RNASeq analyses, a file containing all original non-incremented read counts (counts.YourExperiment.txt) and gene-level total counts (counts.YourExperiment-steady-state.txt) will also be present.