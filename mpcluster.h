#pragma once

#include <set>
#include "abcxyz.h"
#include "seqdb.h"
#include "motifprofile.h"

enum
	{
#define S(x)	SEG_##x,
#include "segs.h"
	};

// Motif profile clustering
class MPCluster
	{
public:
	const vector<MotifProfile *> *m_Input = 0;

	set<uint> m_PendingIndexes;

// Greedy cluster
	float m_MinScore = FLT_MAX;
	uint m_Sample1 = 64;
	uint m_Sample2 = 64;
	uint m_TopN = UINT_MAX;
	vector<uint> m_CentroidIndexes;
	vector<vector<uint> > m_CentroidIndexToMemberIndexes;
	vector<uint> m_ClusterSizeOrder;

// NN cluster
	uint m_InputMPCount = 0;
	vector<MotifProfile *> m_MPs;
	vector<uint> m_Parents;
	vector<uint> m_Lefts;
	vector<uint> m_Rights;
	vector<uint> m_Sizes;

public:
	static SeqDB m_FLSeqDB;
	static SeqDB m_MotifSeqDB;
	static vector<uint> m_PosAs;
	static vector<uint> m_PosBs;
	static vector<uint> m_PosCs;

	static uint m_LeftFlank_Core;	// see also RdRpSearcher
	static uint m_RightFlank_Core;	// see also RdRpSearcher


public:
	void Clear()
		{
		m_Input = 0;
		m_MinScore = FLT_MAX;
		m_PendingIndexes.clear();
		m_CentroidIndexes.clear();
		m_CentroidIndexToMemberIndexes.clear();
		m_ClusterSizeOrder.clear();
		m_MPs.clear();
		m_Parents.clear();
		m_Lefts.clear();
		m_Rights.clear();
		m_Sizes.clear();
		}

	void LogNN() const;
	void NNCluster(const vector<MotifProfile *> &Input,
	  float MinScore);
	void FindNN(uint &Index1, uint &Index2) const;
	void Join(uint Index1, uint Index2);
	MotifProfile &CreateProfileNN(uint Index1, uint Index2) const;
	void NNToNewickFile(const string &FileName) const;
	void NNNodeToNewickFile(FILE *f, uint NodeIndex) const;

	void GreedyCluster(const vector<MotifProfile *> &Input,
	  float MinScore);
	uint GetNextGreedyCentroid() const;
	void GetRandomPending(uint n, vector<uint> &v) const;
	uint GetBestCentroid(const vector<uint> &v) const;
	void GetMembers(uint CentroidIndex, vector<uint> &Members) const;
	float GetScore(const MotifProfile &MP1,
	  const MotifProfile &MP2) const;
	float GetScore_DotProduct(const MotifProfile &MP1,
	  const MotifProfile &MP2) const;
	float GetScoreNNPair(uint i1, uint i2) const;
	void LogPair(const MotifProfile &MP1,
	  const MotifProfile &MP2) const;
	MotifProfile &GetProfile(uint i) const
		{
		const vector<MotifProfile *> &v = *m_Input;
		asserta(i < SIZE(v));
		return *v[i];
		}
	void LogClusters() const;
	void LogCluster(uint i) const;
	void WriteOutput() const;
	void WriteCluster(uint ClusterOrderIndex) const;
	void WriteFasta(const MotifProfile &MP) const;

public:
	static void ReadMPs(const string &FileName,
	  vector<MotifProfile *> &MPs);
	static void ReadSeqsVec(const string &TsvFileName,
	  vector<vector<string> > &SeqsVec);
	static void LogLogos(const vector<MotifProfile *> &MPs);
	static void TrimLeft(string &Left);
	static void TrimRight(string &Right);
	};
