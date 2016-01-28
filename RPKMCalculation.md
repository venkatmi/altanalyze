### How are RPKMs Calculated in AltAnalyze? ###

RPKMs or reads per kilobase per million is calculated similar to other applications. An example is outlined below:

`RPKM = 1,000,000*(Num-Reads)/(All-Reads*Seq-Len)`

For gene X for one sample, let's say you get an RPKM of 30.45. You sequenced a sample at a depth of 100 million reads, paired-end with 100nt reads. To get here, the following variables were calculated:

**Num-Reads/gene X/sample A = 53,504**

**Seq-Len gene X = 1.7kb**

**All-Reads/sample A = 1,033,588,886**

The reason that the All-Reads is so much more than the actual number of reads is AltAnalyze can actually chop up a single read into multiple read calls based on the exon "regions" it defines. When BedTools is called, it counts the number of reads that align to a given exon region. Since a single read can overlap with a lot of regions, you get over-calling at the gene and dataset-level (but not at the individual exon or junction level). This shouldn't cause problems for differential expression analysis since gene reads are biased in the same way.

RPKMs are also calculated for exons (exon-length) and junctions (assuming 60nt for any junction) in the same way. Prior to calculating any RPKM, the counts for the gene are incremented by a single count, to prevent non-zero RPKMs (no fold change calculated).