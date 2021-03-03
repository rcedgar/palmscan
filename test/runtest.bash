#!/bin/bash -e

palmscan=../bin/palmscan

rm -rf results
mkdir -p results

name=Coronavirus_RefSeq_complete_genomes
$palmscan -search_pp data/$name.fna \
  -all \
  -rdrp \
  -rt \
  -ppout results/$name.pp.faa \
  -ppout_nt results/$name.pp.fna \
  -fevout results/$name.fev \
  -report results/$name.txt

name=PF00078_RT
$palmscan -search_pp data/$name.faa \
  -all \
  -rdrp \
  -rt \
  -ppout results/$name.pp.faa \
  -fevout results/$name.fev \
  -report results/$name.txt
  
name=PF00680_RdRP_1
$palmscan -search_pp data/$name.faa \
  -all \
  -rdrp \
  -rt \
  -ppout results/$name.pp.faa \
  -fevout results/$name.fev \
  -report results/$name.txt

name=PF00978_RdRP_2
$palmscan -search_pp data/$name.faa \
  -all \
  -rdrp \
  -rt \
  -ppout results/$name.pp.faa \
  -fevout results/$name.fev \
  -report results/$name.txt

name=PF00998_RdRP_3
$palmscan -search_pp data/$name.faa \
  -all \
  -rdrp \
  -rt \
  -ppout results/$name.pp.faa \
  -fevout results/$name.fev \
  -report results/$name.txt

name=PF02123_RdRP_4
$palmscan -search_pp data/$name.faa \
  -all \
  -rdrp \
  -rt \
  -ppout results/$name.pp.faa \
  -fevout results/$name.fev \
  -report results/$name.txt
