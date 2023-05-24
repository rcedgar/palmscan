#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include <set>

void ReadLinesFromFile(const string &FileName, vector<string> &Lines);

static bool MatchLabelSubstr(const string &QueryLabel,
  const string &RefLabel)
	{
	return (QueryLabel.find(RefLabel) != string::npos);
	}

static bool MatchLabelSet(const string &QueryLabel,
  const set<string> &RefLabels, bool MatchSubstr)
	{
	if (!MatchSubstr)
		{
		bool Found = (RefLabels.find(QueryLabel) != RefLabels.end());
		return Found;
		}
	for (set<string>::const_iterator p = RefLabels.begin();
	  p != RefLabels.end(); ++p)
		{
		bool Found = MatchLabelSubstr(QueryLabel, *p);
		if (Found)
			return true;
		}
	return false;
	}

void cmd_chains_subsample()
	{
	const string &InputFileName = opt_chains_subsample;

	asserta(optset_sample_size);
	uint SampleSize = opt_sample_size;
	asserta(SampleSize > 0);

	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint N = min(SIZE(Chains), SampleSize);
	if (N < SampleSize)
		Warning("sample_size > input");
	random_shuffle(Chains.begin(), Chains.end());

	uint n = 0;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Writing output");

		PDBChain &Q = *Chains[i];
		Q.ToCal(g_fcal);
		}
	}
