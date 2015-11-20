DomStratStats
=============

Stratified statistics for protein domains.
Compute domain q-values and local FDRs, and tiered q-values (combining HMMER3 domain and sequence data).

About DomStratStats
===

The goal of this software is to bring better statistics to protein sequence analysis.  The standard approaches are based on p-values and E-values, but here we introduce q-values and lFDRs (local False Discovery Rates) for protein domains.  They key to make q-values and lFDRs work better is by stratifying, in this case analyzing each domain family separately, which best balances noise across domain families.

Practically, the "domain q-value" approach will be the most useful.  Theoretically, we proved in our work that lFDRs are optimal; however, we found in benchmarks that lFDR estimates do not perform as well as q-values because some HMMER3 p-values are bad (for certain repetitive domain families, in particular).  We found that q-values are much more robust than lFDRs given the p-values that HMMER3 gives.  Our code provides lFDRs anyway for research purposes, and with the hope that they can be better used in the future, when domain p-values improve.

Lastly, we included an even more experimental mode, "tiered q-values".  Our tiered analysis combines the HMMER sequence and domain p-values into a joint assessment of significance based on the FDR.  On the plus side, it produces even more predictions at a fixed combined FDR than the previous approaches.  However, the theoretical FDR is too far off of the empirical FDR, which means that you should exercise caution choosing which threshold to use.

Installation
===

You'll need Perl 5 (without additional Perl packages), [HMMER3](http://hmmer.janelia.org/), and gzip installed.

For the HMM database, you have two options.  For Pfam, you will need to download Pfam-A.hmm.gz and Pfam-A.hmm.dat.gz from [Pfam](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/). For Superfamily, register and download hmmlib_1.75.gz from [Superfamily](http://supfam.org/SUPERFAMILY/downloads.html).  Either way, use HMMER3's hmmpress to prepare database for searching.


Synopsis of scripts
===

Produce domain predictions with HMMER3 (provides hmmscan) with weak filters.  This is the slowest step.
```
perl -w 0runHmmscan.pl hmmscan Pfam-A.hmm PLAF7.fa PLAF7.p.txt 
```

(Optional) compress output, our scripts will read it either way 
```
gzip PLAF7.p.txt 
```

Remove overlaps, ranking by p-value, creating an intermediate output file 
```
perl -w 1noOvs.pl PLAF7.p.txt PLAF7.po.txt Pfam-A.hmm.dat 
```

Add domain q-value, local FDR, and FDR|lFDR columns (no thresholds set) 
```
perl -w 2domStratStats.pl PLAF7.fa PLAF7.po.txt PLAF7.pod.txt 
```

... and set a q-value threshold 
```
perl -w 2domStratStats.pl PLAF7.fa PLAF7.po.txt PLAF7.pod.q4e-4.txt 4e-4 
```

... or set a local FDR threshold 
```
perl -w 2domStratStats.pl PLAF7.fa PLAF7.po.txt PLAF7.pod.l2.5e-2.txt 1 2.5e-2 
```

... or set an FDR threshold via equal local FDR thresholds 
```
perl -w 2domStratStats.pl PLAF7.fa PLAF7.po.txt PLAF7.pod.fl4e-4.txt 1 1 4e-4 
```
Tiered q-values (for sequence and domain p-values), with mandatory threshold set 
```
perl -w 3tieredStratQ.pl PLAF7.fa PLAF7.po.txt PLAF7.pot.1e-4.txt 1e-4 
```

... or threshold set to sequence q-values only 
```
perl -w 3tieredStratQ.pl PLAF7.fa PLAF7.po.txt PLAF7.pot.1e-4,1.txt 1e-4 1 
```

Run entire pipeline on multiple organisms, pooling their stats.  Pass prefixes ORG instead of ORG.fa or ORG.fa.gz. Outputs follow patterns (see detailed examples) 
```
perl -w 4allManyOrgs.pl hmmscan Pfam-A.hmm Pfam-A.hmm.dat Pf Pv Pk Py Pb Pc Bb Ta Tp Tg Nc Et Ch Cp Cm
```

More details
===

All scripts give detailed usage instructions when executed without arguments.  I have more detailed recommendations, usage examples (code snippets and test sample inputs and outputs) at my personal website, viiia.org, in [English](http://viiia.org/domStratStats/?l=en-us) and [Spanish](http://viiia.org/domStratStats/).


Citation
===

2015-11-17. Alejandro Ochoa, John D Storey, Manuel Llinás, and Mona Singh. Beyond the E-value: stratified statistics for protein domain prediction. PLoS Comput Biol. 11 e1004509. [Article](http://dx.doi.org/10.1371/journal.pcbi.1004509) [arXiv](http://arxiv.org/abs/1409.6384) 2014-09-23.
