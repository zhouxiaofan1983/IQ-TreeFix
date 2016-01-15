IQ-TreeFix
==================================================

IQ-TreeFix is a modified version of TreeFix , a tool for statistically informed gene tree error correction using species tree. IQ-TreeFix uses IQ-TREE instead of RAxML (the default statistical module of TreeFix; still available in IQ-TreeFix) for likelihood calculation and statistical test. The motivations are: 1) IQ-TREE supports a much larger selection of evolutionary models for nucleotide and amino acid data, and also supports a number of codon models; 2) the version of RAxML integrated in TreeFix does not support multithreading; 3) the flexibility to calculate per-site likelihood values using IQ-TREE and conduct various statistical tests using CONSEL.

--------------------------------------------------

Please refer to the websites of TreeFix (http://compbio.mit.edu/treefix) and IQ-TREE (http://www.cibiv.at/software/iqtree) for their implementation details.
