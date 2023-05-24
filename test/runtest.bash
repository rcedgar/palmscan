#!/bin/bash -e

palmscan=../bin/palmscan2

rm -rf results
mkdir -p results

for name in \
  Coronavirus_RefSeq_complete_genomes.6f \
  PF00078_RT \
  PF00680_RdRP_1 \
  PF00978_RdRP_2 \
  PF00998_RdRP_3 \
  PF02123_RdRP_4
do
	$palmscan -search_pssms data/$name.fa \
	  -tsv results/$name.tsv \
	  -fev results/$name.fev \
	  -fasta results/$name.pp.fasta \
	  -core results/$name.core.fasta \
	  -report_pssms results/$name.report.txt
done
