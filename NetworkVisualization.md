## Network Analysis and Visualization in AltAnalyze ##

In AltAnalyze version 2.0.8 and above, users can obtain automatically produced network visualization outputs from their analyses as well as perform custom network analyses using a variety of biological annotation resources. These network analyses consist of two functions:
  1. GO-Elite Enriched Gene-Term Associations
  1. NetPerspective hypothetical network generation

## GO-Elite Enriched Gene-Term Associations ##

GO-Elite enrichment analysis produces a set of gene-sets and ontology terms pruned based on various statistical parameters. Given that genes or metabolites can often be shared between distinct enriched terms, visualizing the relationships can be important to determine functional redundancy. GO-Elite will output the relationships as automatically generated networks, visualized with up and downregulated genes (red and blue respectively) produced through AltAnalyze.

Custom input and denominator lists can also be analyzed using the "Additional Analyses" menu option, in which case, nodes will be colored yellow when no quantitative values are included (third column of the input files).

## NetPerspective: AltAnalyze Biological Network Analysis Tool ##

NetPerspective is a new tool in AltAnalyze designed to visualize biological interactions between gene, proteins, RNAs, metabolites and drugs in a coherent and easy way. More information on this tool can be found [here](NetPerspective.md).