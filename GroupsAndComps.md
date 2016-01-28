### _Why do I need to define groups and comparisons?_ ###

_**Answer**_: While some alternative exon analysis methods do not rely on defining groups or comparisons, AltAnalyze uses these to simplify user analyses. Indeed, most experimental studies require defining two or more groups for the purpose of calculating statistical measures (e.g, case and control). Groups and comparisons are designated by the user in the AltAnalyze user interface, by specifying which samples correspond to which groups and which groups should be compared to each other (pairwise comparisons).

_**Excluding Pairwise Comparisons**_

Since pairwise comparisons are sometimes minimally informative, users can also choose to compare all groups to each other. When this is done, the user does not need to designate comparisons, just groups.

_**Excluding Groups**_

In other studies, such as genome wide polymorphism analyses, the sample groups may actually change for each haplotype block (where SNPs are in linkage disequilibrium). Thus, if a user wishes to see which alternative exons correlate with a particular genotype, groups will need to be determined "on-the-fly" or no groups assigned at all. Some methods, such as FIRMA in R support such analyses, as they do not require group assignment. If analyses are run outside of AltAnalyze, the list of significant probesets (PSRs) can be annotated later to identify interesting functional associations using the "Annotate External Results" option. Alternatively, users can choose to not filter their data and arbitrarily assign sample groups. In such a case, the main goal is to export the normalized intensities and analyze this data outside of AltAnalyze. For more information on this option, see the [associated documentation](ReturnAll.md).