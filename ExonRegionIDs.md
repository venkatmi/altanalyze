## _What does the identifier E2-3 indicate and where does this annotation come from??_ ##

_**Answer**_: Along with a probeset identifier, alternative exon results include an exon annotation that indicates the relative position of a probeset in an exon for a gene (column: exons) as well as the position of the probeset in alternative regions of the exon (column: exon-region-ID).

These exon identifiers are defined for each new build of Ensembl, based on the overlap of exons from all Ensembl transcript in addition to mRNAs from the UCSC genome database. Based on the relative positions of known exons, exon block IDs and region IDs are defined. Below is a graphical representation of this process:

![http://www.altanalyze.org/help_files/image059.jpg](http://www.altanalyze.org/help_files/image059.jpg)

For example, if either an exon region ID or exon ID points to E2-3, the probeset is in the second exon within that gene and three indicates some sub-region within that exon.

An exon region ID is composed of two numbers (e.g E2-3), where the first number is the exon block and the second is the region within that block. The exon block indicates the exon position in the gene, 5' to 3', taking into account all annotated transcript exons. Thus, E2 indicates that this is the second exon in the gene after the first intron. The exon regions are individual regions within a specific exon block that indicate alternative 3' or 5' splice sites within the exon based on overlapping exon positions for different transcripts. In an exon with three regions, region 1 likely represents sequence proceeding an alternative 3' splice site, region 2 would be sequence common to both overlapping exons and region 3 would be sequence following an alternative 5' splice site. This is an idealized example, but looking at the splicing annotation column should indicate was is likely going on. Otherwise, check out the exon annotation track in the UCSC genome browser.

The exon ID is similar to the exon-region ID, except that the second number indicates the relative position of the probeset relative to any other probesets that align to that exon. Thus, there may be just on region within an exon, but multiple probesets aligning to that region.