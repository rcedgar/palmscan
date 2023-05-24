#include "myutils.h"
#include "calreader.h"
#include <set>

void SplitTrainTest(const vector<PDBChain *> &Chains, double TrainPct,
  vector<PDBChain *> &TrainChains, vector<PDBChain *> &TestChains)
	{
	const uint N = SIZE(Chains);
	if (N < 2)
		Die("%u chains, cannot split", N);

	uint TrainCount = uint((N*TrainPct)/100);
	if (TrainCount == 0)
		TrainCount = 1;
	uint TestCount = N - TrainCount;

	vector<uint> Indexes;
	for (uint i = 0; i < N; ++i)
		Indexes.push_back(i);

	random_shuffle(Indexes.begin(), Indexes.end());

	set<uint> DoneSet;
	for (uint i = 0; i < TrainCount; ++i)
		{
		uint Index = Indexes[i];
		asserta(DoneSet.find(Index) == DoneSet.end());
		DoneSet.insert(Index);
		PDBChain &Chain = *Chains[Index];
		TrainChains.push_back(&Chain);
		}

	for (uint i = TrainCount; i < N; ++i)
		{
		uint Index = Indexes[i];
		asserta(DoneSet.find(Index) == DoneSet.end());
		DoneSet.insert(Index);
		PDBChain &Chain = *Chains[Index];
		TestChains.push_back(&Chain);
		}
	asserta(SIZE(DoneSet) == N);
}
void cmd_split_train_test()
	{
	const string &InputFileName = opt_split_train_test;
	const string &TrainFileName = opt_train_cal;
	const string &TestFileName = opt_test_cal;

	double TrainPct = 50.0;
	if (optset_train_pct)
		TrainPct = opt_train_pct;
	if (TrainPct <= 0.0 || TrainPct >= 100.0)
		Die("Invalid -train_pct %.3g", TrainPct);

	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);

	FILE *fTest = CreateStdioFile(TestFileName);
	FILE *fTrain = CreateStdioFile(TrainFileName);

	vector<PDBChain *> TrainChains;
	vector<PDBChain *> TestChains;
	SplitTrainTest(Chains, TrainPct, TrainChains, TestChains);

	const uint TrainCount = SIZE(TrainChains);
	const uint TestCount = SIZE(TestChains);

	for (uint i = 0; i < TrainCount; ++i)
		{
		PDBChain &Chain = *TrainChains[i];
		TrainChains.push_back(&Chain);
		}

	for (uint i = 0; i < TestCount; ++i)
		{
		PDBChain &Chain = *TestChains[i];
		TestChains.push_back(&Chain);
		}

	CloseStdioFile(fTest);
	CloseStdioFile(fTrain);
	Progress("%u train, %u test\n", TrainCount, TestCount);
	}
