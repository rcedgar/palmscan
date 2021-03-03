![palmscan](http://drive5.com/images/palmscan_hdr.png)

### Palmscan algorithm

`Palmscan` is software to detect [viral polymerase](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4711277/) palmprint barcode sequences in longer sequences such as virus genomes and ORFs. Palmprints can be used to classify [RNA viruses](https://en.wikipedia.org/wiki/RNA_virus).

### Palmprint sequence

![PALMDB](http://drive5.com/images/palm_structure_figure.png)

### Repository layout

```
palmscan/
  src/               # Source code (C++)
  bin/               # Pre-built binaries for Linux and Windows
  motif_msas/        # Multiple alignments used to train PSSMs
  test/data          # Test data
  test/results       # Test results
  test/runtest.bash  # Script to run tests
```

### Software usage

Example command line:

```
palmscan -search_pp inputfile.fa -rt -rdrp -ppout pp.fa -report pp.txt -fevout pp.fev
```

Type `palmscan -help` for option details.