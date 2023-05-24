#pragma once

#include "ppcaligner.h"
#include "pdbchain.h"
#include "tshit.h"

class PPCSearcher
	{
public:
// Query is non-const so that SS can be set
//  only if there are first-pass hits.
	PDBChain *m_Query = 0;
	const PDBChain *m_Ref = 0;
	bool m_DoClassify = false;
	PPCAligner m_PA;
	const vector<PDBChain *> *m_Refs = 0;
	vector<TSHit> m_Hits;
	vector<uint> m_RefIndexes;
	vector<uint> m_HitOrder;
	vector<double> m_RMSDs;
	vector<double> m_Scores;
	vector<bool> m_RefIsRdRps;

	bool m_Classified = false;
	bool m_ClassifiedAsRdRp = false;
	const TSHit *m_TopHit = 0;
	const TSHit *m_TopOtherHit = 0;
	double m_TopScore = 0;
	double m_ScoreDiff = 0;

public:
	static uint m_MaxKeep;

public:
	void ClearSearch()
		{
		m_Hits.clear();
		m_RefIndexes.clear();
		m_HitOrder.clear();
		m_RMSDs.clear();
		m_Scores.clear();
		m_Classified = false;
		m_ClassifiedAsRdRp = false;
		m_TopHit = 0;
		m_TopOtherHit = 0;
		m_TopScore = 0;
		m_ScoreDiff = 0;
		}

	void InitClassify();
	void SetQuery(PDBChain &Q);
	void SetRefs(const vector<PDBChain *> &Refs);
	bool Search(PDBChain &Q);
	void WriteOutput() const;
	const TSHit *GetTopHit() const;
	bool TopHitIsRdRp() const;
	const uint GetTopClassifiedHitIndex(bool IsRdRp) const;
	void Classify();
	void WriteClassify(FILE *f) const;

public:
	static bool GetIsRdRpFromLabel(const string &Label);
	};
