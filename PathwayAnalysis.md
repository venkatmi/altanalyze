## Performing Pathway Analysis Downstream of AltAnalyze ##

After running AltAnalyze to identify differentially or alternatively expressed genes, it is recommended that the users explore these results along biological pathways. Within AltAnalyze users can select the GO-Elite option to obtain over-representation results for a number of ontologies, pathway and gene-set databases. In addition, the input files created for [GO-Elite](GOElite.md) (GO-Elite/input directory) can be used to directly visualize gene expression or alternative exon results on WikiPathways, using the AltAnalyze interface **Additional Analyses** and **Pathway Visualization**. Since this interface has limited options for pathway visualization, the user may want to use more sophisticated pathway visualization programs. Two examples are listed below.

### GenMAPP-CS ###

GenMAPP-CS can be downloaded from http://www.genmapp.org/beta/genmappcs/. Compatibility and installation information can be found **[here](http://code.google.com/p/go-elite/wiki/Installation)**.

To load the GenMAPP input text file produced by AltAnalyze (ExpressionOutput directory), see the [following tutorials](http://opentutorials.cgl.ucsf.edu/index.php/Portal:GenMAPP-CS).

### PathVisio ###

PathVisio, like GenMAPP-CS, can be used to load pathways from WikiPathways. In addition, it can create new pathways and has a number of sophisticated options for streamlined data visualization.

The same input file used for GenMAPP-CS can also be loaded into PathVisio. This will require that the text file is imported and processed before creating criterion. For example directions, see their [tutorial here](http://www.pathvisio.org/wiki/PathVisioTutorials).