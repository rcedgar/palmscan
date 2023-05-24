#include "myutils.h"
#include "cdsearcher.h"

#define TRACE	0

double GetNormal(double Mu, double Sigma, double x);

void CDSearcher::LogHits() const
	{
	Log("\n");
	Log(">%s\n", m_Query->m_Label.c_str());
	Die("TODO");
	}

void CDSearcher::ClearSearch()
	{
	m_Query = 0;
	}

void CDSearcher::Init(const CDInfo &Info, const CDData &Dists,
  const CDData &StdDevs)
	{
	Clear();
	m_Info = &Info;
	m_Dists = &Dists;
	m_StdDevs = &StdDevs;
	ClearSearch();
	}

double CDSearcher::GetScore(
  uint SeqPos1, uint SeqPos2,
  uint Ix1, uint Ix2,
  uint L1, uint L2) const
	{
	bool Diag = (Ix1 == Ix2);
	double XS = (Diag ? 1.5 : 2);
	double Sum = 0;
	uint n = 0;
	asserta(m_Dists != 0);
	asserta(m_StdDevs != 0);
	const CDData &Dists = *m_Dists;
	const CDData &StdDevs = *m_StdDevs;
	for (uint i = 0; i < L1; ++i)
		{
		uint jhi = (Diag ? i : L2);
		for (uint j = 0; j < jhi; ++j)
			{
			double Observed_d = m_Query->GetDist(SeqPos1+i, SeqPos2+j);
			double Mu = Dists.GetByIx(Ix1 + i, Ix2 + j);
			double Sigma = StdDevs.GetByIx(Ix1 + i, Ix2 + j);
			assert(Sigma > 0);
			double y = GetNormal(Mu, XS*Sigma, Observed_d);
			double Max = GetNormal(Mu, XS*Sigma, Mu);
			double Ratio = y/Max;
			Sum += Ratio;
			++n;
			}
		}
	asserta(n > 0);
	double Score = Sum/n;
	return Score;
	}

/***
Minimum start position of motif is:
	{Position of lowest hit to previous motif} +
	  {minimum number of AAs from previous motif}

Maximum start position of motif is:
	{Position of highest hit to previous motif} +
	  {maximum number of AAs from previous motif}
***/
void CDSearcher::GetRange(uint MotifIndex, 
  uint PrevLo, uint PrevHi, uint &Lo, uint &Hi) const
	{
	Lo = UINT_MAX;
	Hi = UINT_MAX;

	uint QL = m_Query->GetSeqLength();
	uint MinAAEnd = m_Template->GetMinAAsSeqEnd(MotifIndex);
	uint Hi_end = QL - MinAAEnd;
	if (MotifIndex == 0 || PrevLo == UINT_MAX)
		{
		assert(PrevHi = UINT_MAX);
		if (MinAAEnd < QL)
			{
			Lo = 0;
			Hi = Hi_end;
			}
		return;
		}

	assert(PrevLo <= PrevHi);
	uint MinAA = m_Template->GetMinAAsNext(MotifIndex-1);
	uint MaxAA = m_Template->GetMaxAAsNext(MotifIndex-1);

	Lo = PrevLo + MinAA;
	Hi = PrevHi + MaxAA;
	if (Hi >= Hi_end)
		Hi = Hi_end;
	if (Hi < Lo)
		{
		Lo = UINT_MAX;
		Hi = UINT_MAX;
		}
	}

void CDSearcher::Search1(uint MotifIndex, uint Lo, uint Hi,
  vector<uint> &Hits, vector<double> &Scores) const
	{
	const char AnchorAA = m_Template->m_AnchorAAs[MotifIndex];
	const uint AnchorAAOffset = m_Template->m_AnchorAAOffsets[MotifIndex];
	double MinScore = m_Template->m_MinScores[MotifIndex];
#if TRACE
	Log("Search1(Mf=%u, Lo=%u, Hi=%u) aa=%c(%u)\n",
	  MotifIndex, Lo, Hi, AnchorAA, AnchorAAOffset);
#endif
	const string &Seq = m_Query->m_Seq;
	const uint QL = SIZE(Seq);
	asserta(Lo <= Hi && Hi < QL);
	const uint Ix = m_Info->GetIx(MotifIndex, 0);
	const uint ML = m_Info->GetMotifLength(MotifIndex);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		if (AnchorAA != 'x' && Seq[Pos+AnchorAAOffset] != AnchorAA)
			continue;

		double Score = GetScore(Pos, Pos, Ix, Ix, ML, ML);
		if (Score >= MinScore)
			{
#if TRACE
			Log("Mf=%u Pos=%u, score=%.4f hits=%u\n",
			  MotifIndex, Pos, Score, SIZE(Hits));
#endif
			if (!Hits.empty())
				{
				uint LastHit = Hits.back();
				double LastScore = Scores.back();
				if (Pos - LastHit < 10)
					{
					if (Score > LastScore)
						{
						Hits.back() = Pos;
						Scores.back() = Score;
						}
					continue;
					}
				}
			Hits.push_back(Pos);
			Scores.push_back(Score);
			}
		}
	}

void CDSearcher::InitSearch(const PDBChain &Query)
	{
	ClearSearch();
	m_Query = &Query;
	}

uint CDSearcher::GetSeqLength() const
	{
	return m_Query->GetSeqLength();
	}

void CDSearcher::SearchMotifs(const vector<uint> &MotifIndexes,
  vector<uint> &TopHit)
	{
#if TRACE
	Log("\n");
	Log("SearchMotifs(%s)\n", Q.m_Label.c_str());
#endif
	TopHit.clear();
	const uint QL = GetSeqLength();
	const uint NM = SIZE(MotifIndexes);
	asserta(NM > 0);

	vector<vector<uint> > HitVecRagged;
	uint PrevLo = UINT_MAX;
	uint PrevHi = UINT_MAX;
	for (uint m = 0; m < NM; ++m)
		{
		const uint MotifIndex = MotifIndexes[m];

		uint Lo, Hi;
		GetRange(MotifIndex, PrevLo, PrevHi, Lo, Hi);

#if TRACE
		{
		const char *Name = m_Info->GetMotifName(MotifIndex);
		Log("Motif %s, range %u - %u\n", Name, Lo, Hi);
		}
#endif
		if (Lo == UINT_MAX)
			return;

		vector<uint> Hits;
		vector<double> Scores;
		Search1(MotifIndex, Lo, Hi, Hits, Scores);
		const uint HitCount = SIZE(Hits);
#if TRACE
		{
		Log("  %u hits ", SIZE(Hits));
		for (uint i = 0; i < SIZE(Hits); ++i)
			Log(" %u", Hits[i]);
		Log("\n");
		}
#endif
		if (HitCount == 0)
			return;

		HitVecRagged.push_back(Hits);
		PrevLo = Hits.front();
		PrevHi = Hits.back();
		}

	vector<vector<uint> > HitVec;
	ExpandHitVec(HitVecRagged, HitVec);

#if TRACE
	{
	Log("\n");
	Log("Combined hits (%u) >%s\n",
	  SIZE(HitVec), m_Query->m_Label.c_str());
	for (uint i = 0; i < SIZE(HitVec); ++i)
		{
		Log("[%u]  ", SIZE(HitVec[i]));
		for (uint m = 0; m < SIZE(HitVec[i]); ++m)
			Log("  %4u", HitVec[i][m]);
		double Score = GetScoreHit(MotifIndexes, HitVec[i]);
		Log("  %.4f\n", Score);
		}
	}
#endif

	double TopScore = 0;
	for (uint i = 0; i < SIZE(HitVec); ++i)
		{
		double Score = GetScoreHit(MotifIndexes, HitVec[i]);
		if (Score > TopScore)
			{
			TopHit = HitVec[i];
			TopScore = Score;
			}
		}

#if TRACE
	{
	Log("\n");
	if (TopHit.empty())
		{
		Log("No hit >%s\n", m_Query->m_Label.c_str());
		return;
		}

	Log("Top hit >%s %.4f\n", m_Query->m_Label.c_str(), TopScore);
	asserta(SIZE(TopHit) == NM);
	for (uint m = 0; m < NM; ++m)
		{
		Log(" %u", TopHit[m]);
		}
	Log("\n");
	}
#endif
	}

double CDSearcher::GetScoreHit(const vector<uint> &MotifIndexes,
  const vector<uint> &Hit) const
	{
	const uint NM = SIZE(MotifIndexes);
	asserta(SIZE(Hit) == NM);

	double Sum = 0;
	uint PairCount = 0;
	for (uint mi = 0; mi < NM; ++mi)
		{
		uint Motifi = MotifIndexes[mi];
		uint Posi = Hit[mi];
		uint Ixi = m_Info->GetIx(Motifi, 0);
		uint Li = m_Info->GetMotifLength(Motifi);
		for (uint mj = 0; mj < NM; ++mj)
			{
			uint Motifj = MotifIndexes[mj];
			uint Posj = Hit[mj];
			uint Ixj = m_Info->GetIx(Motifj, 0);
			uint Lj = m_Info->GetMotifLength(Motifj);
			double Score = GetScore(Posi, Posj, Ixi, Ixj, Li, Lj);
			Sum += Score;
			++PairCount;
			}
		}
	return Sum/PairCount;
	}

bool CDSearcher::EnumIndexesNext(const vector<uint> &Sizes,
  vector<uint> &Indexes) const
	{
	const uint N = SIZE(Sizes);
	if (Indexes.empty())
		{
		for (uint i = 0; i < N; ++i)
			{
			if (Sizes[i] == 0)
				return false;
			Indexes.push_back(0);
			}
		return true;
		}

	asserta(SIZE(Indexes) == N);
	for (uint i = 0; i < N; ++i)
		{
		++Indexes[i];
		if (Indexes[i] < Sizes[i])
			return true;
		Indexes[i] = 0;
		}

	return false;
	}

void CDSearcher::ExpandHitVec(const vector<vector<uint> > &HitVecRagged,
  vector<vector<uint> > &HitVec)
	{
	HitVec.clear();
	const uint M = SIZE(HitVecRagged);
	if (M == 0)
		return;

	vector<uint> Sizes;
	for (uint m = 0; m < M; ++m)
		{
		uint Size = SIZE(HitVecRagged[m]);
		Sizes.push_back(Size);
		}

	vector<uint> Indexes;
	while (EnumIndexesNext(Sizes, Indexes))
		{
		vector<uint> Hits;
		asserta(SIZE(Indexes) == M);
		for (uint m = 0; m < M; ++m)
			{
			uint Ix = Indexes[m];
			Hits.push_back(HitVecRagged[m][Ix]);
			}
		HitVec.push_back(Hits);
		}
	}

void CDSearcher::AddMotif(
  const vector<uint> &MotifIndexes,
  const vector<uint> &Hits,
  uint MotifIndex, uint &TopHit, double &TopScore) const
	{
	TopHit = UINT_MAX;
	TopScore = 0;

	const uint QL = GetSeqLength();
	uint NM = SIZE(MotifIndexes);
	asserta(NM > 0);

	uint Lo = UINT_MAX;
	uint Hi = UINT_MAX;
	vector<uint> MotifIndexesX;
	vector<uint> HitsX;
	uint MI = UINT_MAX;
	const uint ML = m_Info->GetMotifLength(MotifIndex);
	if (MotifIndex + 1 == MotifIndexes[0])
		{
		uint MinAAsNext = m_Template->GetMinAAsNext(MotifIndex);
		if (MinAAsNext > Hits[0])
			return;
		Hi = Hits[0] - MinAAsNext;
		uint MinAAsStart = m_Template->GetMinAAsSeqStart(MotifIndex);
		Lo = MinAAsStart;
		if (Lo > Hi)
			return;
		MotifIndexesX.push_back(MotifIndex);
		HitsX.push_back(UINT_MAX);
		for (uint i = 0; i < NM; ++i)
			{
			MotifIndexesX.push_back(MotifIndexes[i]);
			HitsX.push_back(Hits[i]);
			}
		MI = 0;
		}
	else if (MotifIndexes[NM-1] + 1 == MotifIndex)
		{
		uint MinAAsNext = m_Template->GetMinAAsNext(MotifIndexes[NM-1]);
		Lo = Hits[NM-1] + MinAAsNext;
		if (Lo >= QL)
			return;
		uint MaxAAsNext = m_Template->GetMaxAAsNext(MotifIndexes[NM-1]);
		Hi = Hits[NM-1] + MaxAAsNext;
		if (Hi + ML >= QL)
			Hi = QL - ML - 1;
		MotifIndexesX = MotifIndexes;
		HitsX = Hits;
		MotifIndexesX.push_back(MotifIndex);
		HitsX.push_back(UINT_MAX);
		MI = NM;
		}
	else
		asserta(false);
	if (Lo + ML >= QL)
		return;

	const char AnchorAA = m_Template->m_AnchorAAs[MotifIndex];
	const uint AnchorAAOffset = m_Template->m_AnchorAAOffsets[MotifIndex];
	double MinScore = m_Template->m_MinScores[MotifIndex];
#if TRACE
	Log("AddMotif(Mf=%u, Lo=%u, Hi=%u) aa=%c(%u)\n",
	  MotifIndex, Lo, Hi, AnchorAA, AnchorAAOffset);
#endif
	const string &Seq = m_Query->m_Seq;
	const uint Ix = m_Info->GetIx(MotifIndex, 0);
	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
		if (AnchorAA != 'x' && Seq[Pos+AnchorAAOffset] != AnchorAA)
			continue;

		double Score1 = GetScore(Pos, Pos, Ix, Ix, ML, ML);
#if TRACE
		Log("Mf=%u Pos=%u, score1=%.4f\n",
			MotifIndex, Pos, Score1);
#endif
		if (Score1 >= MinScore)
			{
			HitsX[MI] = Pos;
			double Score = GetScoreHit(MotifIndexesX, HitsX);
#if TRACE
			Log("Mf=%u Pos=%u, score1=%.4f score=%.4f\n",
				MotifIndex, Pos, Score1, Score);
#endif
			if (Score > TopScore)
				{
				TopScore = Score;
				TopHit = Pos;
				}
			}
		}
	}

double CDSearcher::SearchPalm(const PDBChain &Query,
  vector<uint> &Hit)
	{
	InitSearch(Query);

	vector<uint> MotifsABC;
	MotifsABC.push_back(1);
	MotifsABC.push_back(2);
	MotifsABC.push_back(3);

	vector<uint> TopHitABC;
	SearchMotifs(MotifsABC, TopHitABC);
	if (TopHitABC.empty())
		return 0;

	uint HitF2 = UINT_MAX;
	double ScoreF2 = 0;
	AddMotif(MotifsABC, TopHitABC, 0, HitF2, ScoreF2);
	if (HitF2 == UINT_MAX)
		return 0;

	vector<uint> MotifsF2ABC;
	MotifsF2ABC.push_back(0);
	MotifsF2ABC.push_back(1);
	MotifsF2ABC.push_back(2);
	MotifsF2ABC.push_back(3);

	vector<uint> TopHitF2ABC;
	TopHitF2ABC.push_back(HitF2);
	TopHitF2ABC.push_back(TopHitABC[0]);
	TopHitF2ABC.push_back(TopHitABC[1]);
	TopHitF2ABC.push_back(TopHitABC[2]);

	uint HitD = UINT_MAX;
	double ScoreD = 0;
	AddMotif(MotifsF2ABC, TopHitF2ABC, 4, HitD, ScoreD);
	if (HitD == UINT_MAX)
		return 0;

	vector<uint> MotifsF2ABCD;
	MotifsF2ABCD.push_back(0);
	MotifsF2ABCD.push_back(1);
	MotifsF2ABCD.push_back(2);
	MotifsF2ABCD.push_back(3);
	MotifsF2ABCD.push_back(4);

	vector<uint> TopHitF2ABCD;
	TopHitF2ABCD.push_back(TopHitF2ABC[0]);
	TopHitF2ABCD.push_back(TopHitF2ABC[1]);
	TopHitF2ABCD.push_back(TopHitF2ABC[2]);
	TopHitF2ABCD.push_back(TopHitF2ABC[3]);
	TopHitF2ABCD.push_back(HitD);

	uint HitE = UINT_MAX;
	double ScoreE = 0;
	AddMotif(MotifsF2ABCD, TopHitF2ABCD, 5, HitE, ScoreE);
	if (HitE == UINT_MAX)
		return 0;

	Hit = TopHitF2ABCD;
	Hit.push_back(HitE);
	return ScoreE;
	}
