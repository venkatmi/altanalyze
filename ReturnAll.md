### Prevent Filtering of Probesets ###

If a user wishes to either export all interim normalized intensities, without filtering based on expression or get all alternative exons without filtering, they can change the default filtering options or use the "returnAll" option when running AltAnalyze by [command-line](CommandLineMode.md). To mimic these options in the user interface, set the corresponding below variables equal to the shown values.

Both options apply the following parameters:

```
--dabgp 1 --rawexp 1 --altp 1 --probetype full --altscore 1 --GEcutoff 10000
```

When done through the graphical user interface, manually change these defaults (e.g, set the minimum DABG p-value to < 1). When using the [command-line](CommandLineMode.md) option, add the flag --returnAll.

--returnAll: long variable name “return\_all”, default value for this variable is no. When set to yes, returns all un-filtered alternative exon results by setting all associated filtering parameters to the lowest stringency values. This is equivalent to providing the following flags:  --dabgp 1 --rawexp 1 --altp 1 --probetype full --altscore 1 --GEcutoff 10000. Since this option will output all alternative exon scores for all Ensembl annotated probe sets, the results file will be exceptionally large (>500,000 lines), unless the user has saved previously run alternative exon results (e.g., MADS) to the directory “AltDatabase/filtering” in the AltAnalyze program directory, with a name that matches the analyzed comparison. For example, if the user has a list of 2,000 MADS regulated probe sets for cortex versus cerebellum, then the MADS results should be saved to “AltDatabase/filtering” with the name “Cortex\_vs\_Cerebellum.txt” and in AltAnalyze the CEL file groups should be named Cortex and Cerebellum and the comparison should be Cortex versus Cerebellum. When the filename for a file in the “filtering” directory is contained within the comparison filename (ignoring “.txt”), only these probe sets will be selected when exporting the results. This analysis will produce a results file with all AltAnalyze statistics (default or custom) for just the selected probe sets, independent of the value of each statistic.