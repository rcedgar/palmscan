#include "myutils.h"
#include <set>
#include <map>

void GetConnectedComponents2(const vector<unsigned> &FromVec, const vector<unsigned> &ToVec, 
  vector<vector<unsigned> > &CCs,  bool ShowProgress);

bool IsEdge(const string &S, const string &T, unsigned MaxDiffs)
	{
	const unsigned L = SIZE(S);
	asserta(SIZE(T) == L);
	unsigned d = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		char s = S[i];
		char t = T[i];
		if (s != t)
			{
			++d;
			if (d > MaxDiffs)
				return false;
			}
		}
	asserta(d <= MaxDiffs);
	return true;
	}

unsigned SLCSet(const set<string> &Seqs, map<string, unsigned> &SeqToClusterIndex)
	{
	SeqToClusterIndex.clear();

	vector<string> SeqVec;
	const unsigned N = SIZE(Seqs);
	for (set<string>::const_iterator p = Seqs.begin(); p != Seqs.end(); ++p)
		{
		const string &Seq = *p;
		SeqVec.push_back(Seq);
		}

	vector<unsigned> FromVec;
	vector<unsigned> ToVec;
	for (unsigned i = 0; i < N; ++i)
		{
		FromVec.push_back(i);
		ToVec.push_back(i);

		const string &Seqi = SeqVec[i];
		for (unsigned j = i + 1; j < N; ++j)
			{
			const string &Seqj = SeqVec[j];
			if (IsEdge(Seqi, Seqj, opt_cluster_maxdiffs))
				{
				FromVec.push_back(i);
				ToVec.push_back(j);
				}
			}
		}

	vector<vector<unsigned> > CCs;
	GetConnectedComponents2(FromVec, ToVec, CCs, false);

	unsigned CCCount = SIZE(CCs);
	for (unsigned CCIndex = 0; CCIndex < CCCount; ++CCIndex)
		{
		const vector<unsigned> &CC = CCs[CCIndex];
		const unsigned K = SIZE(CC);
		for (unsigned i = 0; i < K; ++i)
			{
			unsigned SeqIndex = CC[i];
			const string &Seq = SeqVec[SeqIndex];
			asserta(SeqToClusterIndex.find(Seq) == SeqToClusterIndex.end());
			SeqToClusterIndex[Seq] = CCIndex;
			}
		}
	asserta(SIZE(SeqToClusterIndex) == N);
	return CCCount;
	}
