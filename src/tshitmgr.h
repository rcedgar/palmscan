#pragma once

#include "trisearcher.h"

class TSHitMgr
	{
public:
	PDBChain *m_Query = 0;
	vector<TSHit> m_Hits;
	TSHit *m_TopHit = 0;

public:
	void SetQuery(PDBChain &Query);
	void Add(TSHit &TH);
	void WriteOutput();
	void SetTopHit();
	void WriteReport(FILE *f) const;
	void WritePPC(FILE* f) const;
	void GetPPC(PDBChain &ChainPPC) const;
	};
