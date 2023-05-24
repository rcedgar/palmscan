#pragma once

#include "xdpmem.h"
#include "pathinfo.h"
#include "pdbchain.h"
#include "tshit.h"

class TSHit;

class TriSearcher
	{
public:

	PDBChain *m_Query = 0;
	const PDBChain *m_Ref = 0;
	vector<uint> m_PosAs;
	vector<uint> m_PosBs;
	vector<uint> m_PosCs;
	vector<double> m_TriRMSD2s;
	vector<double> m_MotifRMSD2s;
	vector<uint> m_HitOrder;

	XDPMem m_XDPMem;

	vector<double> m_TriForm_t;
	vector<vector<double> > m_TriForm_R;

	float **m_DPScoreMx = 0;

public:
	void Clear()
		{
		m_Query = 0;
		m_Ref = 0;
		m_PosAs.clear();
		m_PosBs.clear();
		m_PosCs.clear();
		m_TriRMSD2s.clear();
		m_MotifRMSD2s.clear();
		m_HitOrder.clear();
		}

	void LogMe(FILE *f = g_fLog) const;
	void Search(PDBChain &Query, const PDBChain &Ref);
	void SetHitOrder();
	double GetRMSDMotifs(uint QueryPosA, uint QueryPosB,
	  uint QueryPosC) const;
	double GetRMSD2Segment(uint QueryStartPos, uint RefStartPos, uint n,
	  const vector<double> &t, const vector<vector<double> > &R) const;
	double GetMotifsTM() const;
	double GetTMSum(uint QPos, uint RPos, uint n) const;
	bool GetTopHit(uint &QueryPosA, uint &QueryPosB, uint &QueryPosC) const;
	bool AlignPalm(uint QueryPosA, uint QueryPosB, uint QueryPosC,
	  string &Path);
	void AlignSeg(uint QPos, uint QLen, uint RPos, uint RLen,
	  const vector<double> &t, const vector<vector<double> > &R,
	  string &Path);
	void SetTriForm();
	void AllocDPScoreMx();
	void WriteOutput();
	void WriteTsv();
	void WriteReport();
	void WriteAln(FILE *f);
	void GetHit(uint Index, TSHit &TH) const;
	bool GetTopHit(TSHit &TH) const;
	};
