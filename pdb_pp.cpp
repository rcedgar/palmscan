#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "rdrpsearcher.h"
#include "outputfiles.h"
#include "abcxyz.h"

// Extracts palmprint + margin from PDB file.
// Motifs are identified by PSSMs.
// If there are multiple chains in the input PDB, the first one 
//   with a hit is reported.

static void WriteSegment(FILE *f, const vector<vector<string> > &ATOMs,
  uint Pos, uint L, uint OutPos)
	{
	if (f == 0)
		return;

	const uint N = SIZE(ATOMs);
	asserta(Pos + L <= N);
	for (uint i = 0; i < L; ++i)
		{
		const vector<string> &Lines = ATOMs[Pos+i];

		const uint M = SIZE(Lines);
		for (uint j = 0; j < M; ++j)
			{
			const string &Line = Lines[j];

			string OutputLine;
			PDBChain::SetResidueNrInATOMLine(Line, OutPos+i+1, OutputLine);
			fputs(OutputLine.c_str(), f);
			fputc('\n', f);
			}
		}
	}

void cmd_pdb_pp()
	{
	const string &FN = opt_pdb_pp;
	asserta(optset_output);
	FILE *fOut = CreateStdioFile(opt_output);

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);
	uint FlankSize = 0;
	if (optset_flanks)
		FlankSize = opt_flanks;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;
	RdRpModel Model;
	Model.FromModelFile(ModelFileName);
	//RdRpSearcher::InitOutput();

	const uint ChainCount = SIZE(Chains);
	for (uint i = 0; i < ChainCount; ++i)
		{
		PDBChain &Chain = *Chains[i];
		const string &QLabel = Chain.m_Label;
		const string &QSeq = Chain.m_Seq;
		const uint QL = SIZE(QSeq);

		RdRpSearcher RS;
		RS.Init(Model);
		RS.Search(QLabel, QSeq);
		RS.WriteOutput();
		if (!RS.IsHit())
			continue;

		fprintf(fOut, "TITLE %s\n", Chain.m_Label.c_str());

		uint APos = RS.GetMotifPos(0);
		uint BPos = RS.GetMotifPos(1);
		uint CPos = RS.GetMotifPos(2);
		if (CPos < APos)
			Die("Permuted motifs not supported");

		uint PPStart = APos;
		if (PPStart >= FlankSize)
			PPStart -= FlankSize;
		else
			PPStart = 0;
		uint PPEnd = CPos + CL - 1 + FlankSize;
		if (PPEnd >= QL)
			PPEnd = QL - 1;
		asserta(PPEnd > PPStart);
		uint PPL = PPEnd - PPStart + 1;

		WriteSegment(fOut, Chain.m_ATOMs, PPStart, PPL, 0);
		break;
		}
	CloseStdioFile(fOut);
	}
