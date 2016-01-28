## Introduction ##

AltInteract is a currently undeveloped plugin for Cytoscape in the early planning stage. The purpose of AltInteract will be to identify and visualize the effects of alternative splicing on biological pathways, with an emphasis on domain-domain disruptions. This plugin will interpret predictions from AltAnalyze for genes (e.g., increased/decreased nonsense mediated decay, protein domain disruption, protein binding site removal, microRNA binding site removal, protein truncation) and extend these to interactions observed on pathways. When domain disruptions are identified, known domain-domain interactions between proteins for which protein interactions exist in a pathway will be examined to determine whether an interaction is valid given alternative exon changes.

As an additional functionality, this plugin will include the option of creating new interaction networks from alternatively expressed AltAnalyze genes via the plugin [GeneMANIA](http://www.genemania.org/plugin/). To obtain an overall of all disrupted interactions for all pathways, a separate project will be to obtain the same associations in the software [GO-Elite](GOElite.md). Future development is planned in association with AltAnalyze development and through Google Summer Code projects.

![http://www.altanalyze.org/image/ERbB-example.jpg](http://www.altanalyze.org/image/ERbB-example.jpg)

Mockup of AltInteract visualization upon the WikiPathways [ErbB signaling pathway](http://www.wikipathways.org/index.php/Pathway:WP673).