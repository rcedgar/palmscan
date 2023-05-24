
const char *help_txt[] =
	{
	"Command:",
	"    -search_pssms seqs.fasta    # input sequences, must be aa",
	"",
	"Optional arguments:",
	"    -report_pssms report.txt    # human-readable text report",
	"    -tsv hits.tsv               # hits in tab-separated text",
	"	-fev hits.fev               # hits in field-equals-value",
	"    -threads N                  # threads, default min(number of cores, 10)",
	"	-fasta pp.fa                # pamlmprints, FASTA",
	"	-core core.fa               # palmcore, FASTA",
	"	-minflanklen n              # minimum flank for palmcore (default 0)",
	"",
	"Typical usage:",
	"    palmscan -search_pssms seqs.fasta -tsv hits.tsv",
	};
