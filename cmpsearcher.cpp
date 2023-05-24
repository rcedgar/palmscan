#include "myutils.h"
#include "cmpsearcher.h"

/***
From D:\src\py\palmscan3d_params.py
-----------------------------------
static const uint min_aadist_AdBg = 21;
static const uint max_aadist_AdBg = 136;

static const uint min_aadist_BgCd = 28;
static const uint max_aadist_BgCd = 140;

static const uint min_aadist_AdCd = 85;
static const uint max_aadist_AdCd = 173;
***/

// Manually extended to allow outliers
static const uint min_aadist_AdBg = 10;
static const uint max_aadist_AdBg = 150;

static const uint min_aadist_BgCd = 10;
static const uint max_aadist_BgCd = 150;

//static const uint min_aadist_AdCd = 70;
//static const uint max_aadist_AdCd = 200;

static const double MINSCORE1 = 0.4;
static const double MINSCORE3 = 0.5;

static bool g_Trace = false;
static bool g_Trace2 = false;

void CMPSearcher::SetProfRef(const CMP &Prof, uint RefIndex)
	{
	asserta(RefIndex < SIZE(Prof.m_DistMxVec));
	m_DistMx = &Prof.m_DistMxVec[RefIndex];
	m_StdDevs = &Prof.m_StdDevs;
	}

void CMPSearcher::SetProf(const CMP &Prof)
	{
	m_DistMx = &Prof.m_RefMeans;
	m_StdDevs = &Prof.m_StdDevs;
	}

bool CMPSearcher::GoodScore1(double Score) const
	{
	if (Score >= MINSCORE1)
		return true;
	return false;
	}

bool CMPSearcher::GoodScore3(double Score) const
	{
	if (Score >= MINSCORE3)
		return true;
	return false;
	}

bool CMPSearcher::MatchAd(uint Pos) const
	{
	if (Pos < 3)
		return false;
	double Score = GetScoreA(*m_Query, Pos-3);
	if (GoodScore1(Score))
		return true;
	return false;
	}

bool CMPSearcher::MatchBg(uint Pos) const
	{
	if (Pos < 1)
		return false;
	double Score = GetScoreB(*m_Query, Pos-1);
	if (GoodScore1(Score))
		return true;
	return false;
	}

bool CMPSearcher::MatchCd(uint Pos) const
	{
	if (Pos < 3)
		return false;
	double Score = GetScoreC(*m_Query, Pos-3);
	if (GoodScore1(Score))
		return true;
	return false;
	}

void CMPSearcher::SearchAd(uint AdLo, uint AdHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = 3; Pos + AL+BL+CL < L; ++Pos)
		{
		if (m_Seq[Pos] != 'D')
			continue;
		bool Ok = MatchAd(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void CMPSearcher::SearchBg(uint BgLo, uint BgHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = BgLo; Pos <= BgHi && Pos + BL+CL < L; ++Pos)
		{
		if (m_Seq[Pos] != 'G')
			continue;
		bool Ok = MatchBg(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void CMPSearcher::SearchCd(uint CdLo, uint CdHi,
  vector<uint> &PosVec)
	{
	PosVec.clear();
	const uint L = SIZE(m_Seq);
	for (uint Pos = CdLo; Pos <= CdHi && Pos + 4 < L; ++Pos)
		{
		if (m_Seq[Pos] != 'D')
			continue;
		bool Ok = MatchCd(Pos);
		if (Ok)
			PosVec.push_back(Pos);
		}
	}

void CMPSearcher::GetBgLoHi(uint Ad, uint &BgLo, uint &BgHi) const
	{
	BgLo = Ad + min_aadist_AdBg;
	BgHi = Ad + max_aadist_AdBg;
	}

void CMPSearcher::GetCdLoHi(uint Bg, uint &CdLo, uint &CdHi) const
	{
	CdLo = Bg + min_aadist_BgCd;
	CdHi = Bg + max_aadist_BgCd;
	}

void CMPSearcher::Search(const PDBChain &Query)
	{
	ClearSearch();
	m_Query = &Query;
	m_Seq = Query.m_Seq;

	if (!opt_permuted)
		Search_ABC(Query);

	if (!opt_notpermuted)
		Search_CAB(Query);
	}

void CMPSearcher::Search_ABC(const PDBChain &Query)
	{
	const uint QL = SIZE(Query.m_Seq);

	vector<uint> AdVec;
	vector<uint> BgVec;
	vector<uint> CdVec;

	SearchAd(0, QL, AdVec);
	const uint nad = SIZE(AdVec);
	for (uint iad = 0; iad < nad; ++iad)
		{
		uint Ad = AdVec[iad];
		uint BgLo, BgHi;
		GetBgLoHi(Ad, BgLo, BgHi);

		SearchBg(BgLo, BgHi, BgVec);
		const uint nbg = SIZE(BgVec);
		for (uint ibg = 0; ibg < nbg; ++ibg)
			{
			uint Bg = BgVec[ibg];
			uint CdLo, CdHi;
			GetCdLoHi(Bg, CdLo, CdHi);

			SearchCd(CdLo, CdHi, CdVec);
			const uint ncd = SIZE(CdVec);
			for (uint icd = 0; icd < ncd; ++icd)
				{
				uint Cd = CdVec[icd];
				double Score = GetScore(Ad, Bg, Cd);
				if (g_Trace)
					Log("[Score3] %.4f  (%u, %u, %u)\n",
					  Score, Ad, Bg, Cd);
				CheckHit(Ad, Bg, Cd, Score);
				}
			}
		}
	}

double CMPSearcher::GetPSSMStarts(uint &PosA, uint &PosB, uint &PosC) const
	{
	PosA = UINT_MAX;
	PosB = UINT_MAX;
	PosC = UINT_MAX;
	const uint N = SIZE(m_Scores);
	if (N == 0)
		return 0;
	double TopScore = m_Scores[0];
	PosA = m_Ads[0];
	PosB = m_Bgs[0];
	PosC = m_Cds[0];

	for (uint i = 1; i < N; ++i)
		{
		bool Better = (m_Scores[i] > TopScore);
		if (Better)
			{
			TopScore = m_Scores[i];
			PosA = m_Ads[i];
			PosB = m_Bgs[i];
			PosC = m_Cds[i];
			}
		}

	asserta(PosA >= 3);
	asserta(PosB >= 1);
	asserta(PosC >= 3);
	PosA -= 3;
	PosB -= 1;
	PosC -= 3;

	const string &Q = m_Query->m_Seq;
	asserta(Q[PosA+3] == 'D');
	asserta(Q[PosB+1] == 'G');
	asserta(Q[PosC+3] == 'D');

	return TopScore;
	}

double CMPSearcher::GetScore(uint Ad, uint Bg, uint Cd) const
	{
	const uint QL = SIZE(m_Query->m_Seq);
	asserta(Ad >= 3 && Ad < QL);
	asserta(Bg >= 2 && Bg < QL);
	asserta(Cd >= 3 && Cd < QL);
	double Score = GetScore3(*m_Query, Ad-3, Bg-1, Cd-3);
	return Score;
	}

void CMPSearcher::CheckHit(uint Ad, uint Bg, uint Cd, double Score)
	{
	bool Ok = (Ad < Bg && Bg < Cd) || (Cd < Ad && Ad < Bg);
	if (!Ok)
		return;
	if (!GoodScore3(Score))
		return;

	m_Ads.push_back(Ad);
	m_Bgs.push_back(Bg);
	m_Cds.push_back(Cd);
	m_Scores.push_back(Score);
	}

double CMPSearcher::GetScoreA(const PDBChain &Chain, uint PosA) const
	{
	asserta(Chain.m_Seq[PosA+3] == 'D');
	double Score = GetScore(Chain, PosA, AIX, AL);
	return Score;
	}

double CMPSearcher::GetScoreB(const PDBChain &Chain, uint PosB) const
	{
	asserta(Chain.m_Seq[PosB+1] == 'G');
	double Score = GetScore(Chain, PosB, BIX, BL);
	return Score;
	}

double CMPSearcher::GetScoreC(const PDBChain &Chain, uint PosC) const
	{
	asserta(Chain.m_Seq[PosC+3] == 'D');
	double Score = GetScore(Chain, PosC, CIX, CL);
	return Score;
	}

double CMPSearcher::GetScore(const PDBChain &Chain, uint SeqPos,
  uint Ix, uint L) const
	{
	double Score = GetScore2(Chain, SeqPos, SeqPos, Ix, Ix, L, L);
	return Score;
	}

double CMPSearcher::GetScore2(const PDBChain &Chain,
  uint SeqPos1, uint SeqPos2,
  uint Ix1, uint Ix2,
  uint L1, uint L2) const
	{
	bool Diag = (Ix1 == Ix2);
	double XS = (Diag ? 1.5 : 2);
	double Sum = 0;
	uint n = 0;
	asserta(m_DistMx != 0);
	asserta(m_StdDevs != 0);
	const vector<vector<double> > &DistMx = *m_DistMx;
	const vector<vector<double> > &StdDevs = *m_StdDevs;
	for (uint i = 0; i < L1; ++i)
		{
		uint jhi = (Diag ? i : L2);
		for (uint j = 0; j < jhi; ++j)
			{
			double Observed_d = Chain.GetDist(SeqPos1+i, SeqPos2+j);
			double Mu = DistMx[Ix1+i][Ix2+j];
			double Sigma = StdDevs[Ix1+i][Ix2+j];
			double y = GetNormal(Mu, XS*Sigma, Observed_d);
			double Max = GetNormal(Mu, XS*Sigma, Mu);
			double Ratio = y/Max;
			Sum += Ratio;
			++n;

			if (g_Trace2)
				{
				char MotifChari = CMP::GetMotifChar(Ix1);
				char MotifCharj = CMP::GetMotifChar(Ix2);

				char ci = Chain.m_Seq[SeqPos1+i];
				char cj = Chain.m_Seq[SeqPos2+j];

				Log("%c[%3u]%c", MotifChari, SeqPos1+i, ci);
				Log("  %c[%3u]%c", MotifCharj, SeqPos2+j, cj);
				Log("  d %6.2f", Observed_d);
				Log("  mu %6.2f", Mu);
				Log("  sd %6.2f", Sigma);
				Log("  r %6.4f", Ratio);
				Log("\n");
				}
			}
		}
	asserta(n > 0);
	double Score = Sum/n;
	if (g_Trace)
		{
		char MotifChari = CMP::GetMotifChar(Ix1);
		char MotifCharj = CMP::GetMotifChar(Ix2);

		Log("%c[%3u]", MotifChari, SeqPos1);
		Log("  %c[%3u]", MotifCharj, SeqPos2);
		Log("  %6.4f\n", Score);
		}
	return Score;
	}

double CMPSearcher::GetScoreAB(const PDBChain &Chain, uint PosA, uint PosB) const
	{
	double ScoreAB = GetScore2(Chain, PosA, PosB, AIX, BIX, AL, BL);
	return ScoreAB;
	}

double CMPSearcher::GetScoreBC(const PDBChain &Chain, uint PosB, uint PosC) const
	{
	double ScoreBC = GetScore2(Chain, PosB, PosC, BIX, CIX, BL, CL);
	return ScoreBC;
	}

double CMPSearcher::GetScoreAC(const PDBChain &Chain, uint PosA, uint PosC) const
	{
	double ScoreAC = GetScore2(Chain, PosA, PosC, AIX, CIX, AL, CL);
	return ScoreAC;
	}

double CMPSearcher::GetScore3(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC) const
	{
	double Score = 0;

	Score += GetScoreA(Chain, PosA);
	Score += GetScoreB(Chain, PosB);
	Score += GetScoreC(Chain, PosC);

	Score += GetScoreAB(Chain, PosA, PosB);
	Score += GetScoreBC(Chain, PosB, PosC);
	Score += GetScoreAC(Chain, PosA, PosC);

	Score /= 6;
	return Score;
	}

double CMPSearcher::SearchRefs(const PDBChain &Query, const CMP &Prof,
  uint &BestPosA, uint &BestPosB, uint &BestPosC, string &BestRefLabel)
	{
	const uint RefCount = SIZE(Prof.m_DistMxVec);
	double BestScore = 0;
	BestPosA = UINT_MAX;
	BestPosB = UINT_MAX;
	BestPosC = UINT_MAX;
	BestRefLabel.clear();
	for (uint RefIndex = 0; RefIndex < RefCount; ++RefIndex)
		{
		SetProfRef(Prof, RefIndex);
		Search(Query);
		uint PosA, PosB, PosC;
		double Score = GetPSSMStarts(PosA, PosB, PosC);
		if (Score > BestScore)
			{
			BestPosA = PosA;
			BestPosB = PosB;
			BestPosC = PosC;
			BestRefLabel = Prof.m_RefLabels[RefIndex];
			BestScore = Score;
			}
		}
	return BestScore;
	}
