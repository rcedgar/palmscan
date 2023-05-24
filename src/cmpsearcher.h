#pragma once

#include "pdbchain.h"
#include "cmp.h"

class DSHit;

class CMPSearcher
	{
public:
	const PDBChain *m_Query = 0;
	string m_Seq;
	const vector<vector<double> > *m_DistMx = 0;
	const vector<vector<double> > *m_StdDevs = 0;

	vector<uint> m_Ads;
	vector<uint> m_Bgs;
	vector<uint> m_Cds;
	vector<double> m_Scores;

public:
	void ClearSearch()
		{
		m_Query = 0;
		m_Seq.clear();
		m_Ads.clear();
		m_Bgs.clear();
		m_Cds.clear();
		m_Scores.clear();
		}

	void SetProfRef(const CMP &Prof, uint RefIndex);
	void SetProf(const CMP &Prof);
	bool GoodScore1(double Score) const;
	bool GoodScore3(double Score) const;
	void Search(const PDBChain &Query);
	void Search_ABC(const PDBChain &Query);
	void Search_CAB(const PDBChain &Query);
	double SearchRefs(const PDBChain &Query, const CMP &Prof,
	  uint &PosA, uint &PosB, uint &PosC, string &BestRefLabel);
	double GetPSSMStarts(uint &PosA, uint &PosB, uint &PosC) const;

	void SearchAd(uint AdLo, uint AdHi, vector<uint> &PosVec);
	void SearchBg(uint BgLo, uint BgHi, vector<uint> &PosVec);
	void SearchCd(uint CdLo, uint CdHi, vector<uint> &PosVec);

	bool MatchAd(uint AdPos) const;
	bool MatchBg(uint BgPos) const;
	bool MatchCd(uint CdPos) const;

	void GetBgLoHi(uint Ad, uint &BgLo, uint &BgHi) const;
	void GetCdLoHi(uint Bg, uint &CdLo, uint &CdHi) const;

	void GetAdLoHi_Permuted(uint Cd, uint &AdLo, uint &AdHi) const;
	void GetBgLoHi_Permuted(uint Ad, uint &BgLo, uint &BgHi) const;

	double GetScore(uint Ad, uint Bg, uint Cd) const;
	void CheckHit(uint Ad, uint Bg, uint Cd, double Score);

	double GetScore3(const PDBChain &Chain,
	  uint PosA, uint PosB, uint PosC) const;
	double GetScoreA(const PDBChain &Chain, uint PosA) const;
	double GetScoreB(const PDBChain &Chain, uint PosB) const;
	double GetScoreC(const PDBChain &Chain, uint PosC) const;
	double GetScoreAB(const PDBChain &Chain, uint PosA, uint PosB) const;
	double GetScoreBC(const PDBChain &Chain, uint PosB, uint PosC) const;
	double GetScoreAC(const PDBChain &Chain, uint PosA, uint PosC) const;
	double GetScore(const PDBChain &Chain,
	  uint SeqPos, uint Ix, uint L) const;
	double GetScore2(const PDBChain &Chain,
	  uint SeqPos1, uint SeqPos2,
	  uint Ix1, uint Ix2,
	  uint L1, uint L2) const;
	};
