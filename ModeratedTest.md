## Moderated Tests Implemented in AltAnalyze ##

AltAnalyze provides a moderated t-test option (unpaired, assuming equal variance), based on the limma empirical Bayes model. This model calculates an adjusted variance based on the variance of all genes or exons analyzed for each comparison. Unlike a standard limma analysis, only the variances for the two groups being compared are used to compute the adjusted variance (s2 prior and df prior) for each comparison as opposed to variance from all groups analyzed (when greater than 2). For these calculations (gamma functions), AltAnalyze uses the python mpmath library (http://code.google.com/p/mpmath).


## Conventional Tests Implemented in AltAnalyze ##

In addition to these moderated tests, multiple conventional statistical tests are provided for the various comparison analyses available (e.g. differential gene expression, alternative exon expression). For group comparisons involving greater than two groups, a standard f-test statistic is used. For pairwise group comparisons, users can choose from unpaired and paired t-test’s (assuming equal variance), rank sum, Mann Whitney and Kolmogorov Smirnov.  These tests are made available from the open-source statistical package [SalStat](http://salstat.sourceforge.net/).


## False Discovery Rate P-values ##

Adjusted p-values for these various tests are calculated using the Benjamini-Hochberg correction method.