#!/bin/bash
### Documentation
# http://genomespot.blogspot.com/2014/09/data-analysis-step-8-pathway-analysis.html
#
#sh DESeq22GSEAPreranked.sh /Users/tiphainecmartin/Documents/Projects_MSSM/ParsonsLab/Elias/github/parsonslab/EZH2degrader/results/RNAseq/condition_degraderEZH2_results_ensembl_rankedGSEA.csv
#

DGE=$1
#RNK=`echo $DGE | sed 's/.csv/_inv.rnk/'`
#sed 1d $DGE \
#| sort -k7g \
#| awk -F "," '{ if ($4>0) printf "%s\t%4.3e\n", $10, 1/$7 ;
#else printf "%s\t%4.3e\n", $10, -1/$7 }' \
#| sort -k2gr > $RNK; 

RNK2=`echo $DGE | sed 's/.csv/_log.rnk/'`
sed 1d $DGE \
| sort -k7g \
| awk -F "," '{ if ($4>0) printf "%s\t%4.3e\n", $10, log($7) ;
else printf "%s\t%4.3e\n", $10, -log($7) }' \
| sort -k2gr > $RNK2
