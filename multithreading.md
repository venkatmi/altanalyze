# Multithreaded Analysis in AltAnalyze #

AltAnalyze (version 2.0.8 and greater) has the built-in ability to use multiple cores for certain analyses. If multiple processors are not present, only a single process will be used. The places were this is currently implemented include:
  1. RNA-Seq genomic coordinate exon/junction annotation/assignment (NPU = NPA\*2)
  1. GO-Elite standard enrichment analysis (NPU = 4)
  1. GO-Elite parallel enrichment analysis (NPU = NPA\*2)
  1. GO-Elite based WikiPathways visualization (NPU = NPA\*1.5)

_NPU = number of processes used_

_NPA = number of processors available_

The above defaults were selected to optimize performance on most machines but is not currently customizable to include user-defined number of processes/processors to use. The two different GO-Elite enrichment options above are employed when analyzing all, both or one gene-set type (option 2) or multiple distinct gene-set types (option 3 - atypical). If using the python source code version of AltAnalyze these can be manually indicated in the code before running by finding:

**RNA-Seq coordinate annotation:**

`Line 947 RNASeq.py: pool_size = mlp.cpu_count() * 2`

**Parallel ORA processing:**

` Line 1172 GO-Elite.py: processes=mlp.cpu_count() * 2`

**ORA statistics:**

`Line 810 mappfinder.py: pool = mlp.Pool(processes=4) and 811: si = (len(pathway_db)/4)`

**WP Visualization:**

`Line 1236 GO-Elite.py: pool_size = mlp.cpu_count() * 2`

In future versions of AltAnalyze, we will provide user control over these options.