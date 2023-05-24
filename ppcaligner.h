#pragma once

#include "pdbchain.h"
#include "tshit.h"
#include "xdpmem.h"

class PPCAligner
	{
	PDBChain *m_Q = 0;
	const PDBChain *m_R = 0;

	vector<double> m_RPt;
	vector<double> m_QPt;

	XDPMem m_Mem;
	float **m_DPScoreMx = 0;
	string m_Path;

public:
	PPCAligner()
		{
		Clear();

		float **AllocDPScoreMx();
		m_DPScoreMx = AllocDPScoreMx();
		}

	void Clear()
		{
		m_Q = 0;
		m_R = 0;
		m_RPt.clear();
		m_QPt.clear();
		m_RPt.resize(3);
		m_QPt.resize(3);
		m_Path.clear();
		}

	void SetQuery(PDBChain &Q);
	void SetRef(const PDBChain &R);
	double GetRMSD2Segment(uint QPos, uint RPos, uint n);
	double GetMotifRMSD();
	void Align(TSHit &Hit) const;
//	double GetTMScore(const TSHit &Hit, string &Path) const;
	double GetTMScore();
	};
