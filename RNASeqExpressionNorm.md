### _How is RNA-Seq data normalized in AltAnalyze?_ ###

_**Answer:**_ Typically, RNA-seq data imported into AltAnalyze will be in one of two formats:
# junction coordinates and read counts
# exon coordinate read counts

In both cases, the number of read counts is summarized at either the level of an annotated/novel exon junction, rather than importing the genomic position of each individual read ([details](RNASeq.md)). In AltAnalyze, these summarized features can be analyzed directly as read counts or normalized using a basic quantile normalization approach or RPKM (Reads Per Kilobase of sequence per Million mapped reads). For junctions, sequence length in this equation is set to 60nt. For exons, this length is the length of the full annotated exon (Ensembl/UCSC).

To calculate differential gene and alternative exon expression, some adjustments to these counts must be made, prior to normalization and analysis. For example, if a junction is undetected in condition 1 and highly expressed in condition 2, to obtain a fold change, both conditions must be incremented by a non-zero value to obtain a semi-realistic fold change. We currently increment each read count for an imported exon or junction by 1 (e.g., a non-expressed junction will have a read count of 1). An exon junction would have an RPKM of 0.955 (based on a single read), whereas an exon would have a value dependent on the exon length, typically close to zero. When expression is reported in the gene expression summary results (ExpressionOutput) or alternative expression results (AltResults/AlternativeOutput), the original non-incremented read counts will be reported if analyzing counts or quantile normalized counts. When analyzing RPKM, the incremented RPKM values will be reported.

For users wishing to view the non-incremented read counts for each exon and junction, along with the feature genomic coordinates, see the "counts." file in the ExpressionInput directory.