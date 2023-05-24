#pragma once

#include "pdbchain.h"
#include "cddata.h"
#include "cdtemplate.h"

class CDSearcher
	{
public:
	const PDBChain *m_Query = 0;
	const CDInfo *m_Info = 0;
	const CDData *m_Dists = 0;
	const CDData *m_StdDevs = 0;
	const CDTemplate *m_Template = 0;

public:
	void Clear()
		{
		m_Query = 0;
		m_Info = 0;
		m_Dists = 0;
		m_StdDevs = 0;
		}

	void InitSearch(const PDBChain &Query);
	void LogHits() const;
	void ClearSearch();
	uint GetSeqLength() const;
	void Init(const CDInfo &Info, const CDData &Dists,
	  const CDData &StdDevs);
	double GetScore(uint SeqPos1, uint SeqPos2,
	  uint Ix1, uint Ix2,
	  uint L1, uint L2) const;
	double GetScoreHit(const vector<uint> &MotifIndexes,
	  const vector<uint> &Hit) const;
	void GetRange(uint NextMotifIndex, 
	  uint PrevLo, uint PrevHi, uint &Lo, uint &Hi) const;
	void Search1(uint MotifIndex, uint Lo, uint Hi,
	  vector<uint> &Hits, vector<double> &Scores) const;
	void SearchMotifs(const vector<uint> &MotifIndexes,
	  vector<uint> &TopHit);
	void ExpandHitVec(const vector<vector<uint> > &HitVecRagged,
	  vector<vector<uint> > &HitVec);
	bool EnumIndexesNext(const vector<uint> &Sizes,
	  vector<uint> &Indexes) const;
	void AddMotif(
	  const vector<uint> &MotifIndexes,
	  const vector<uint> &Hits,
	  uint MotifIndex, uint &TopHit, double &TopScore) const;
	double SearchPalm(const PDBChain &Chain,
	  vector<uint> &Hit);
	};
