#include "myutils.h"
#include "rdrpsearcher.h"
#include "sort.h"
#include "abcxyz.h"

void RdRpSearcher::SearchMotif(uint GroupIndex, uint MotifIndex,
  uint QLo, uint QHi, float MinScore, vector<RPHit> &Hits)
	{
	if (m_Trace)
		{
		string GroupName;
		GetGroupName(GroupIndex, GroupName);
		Log("SearchMotif(%u=%s, %c, QLo=%u, QHi=%u, MinScore=%.1f)\n",
		  GroupIndex, GroupName.c_str(), MotLet, QLo, QHi, MinScore);
		}
	Hits.clear();

	asserta(QLo < QHi);
	const uint QL = SIZE(m_QuerySeq);
	asserta(QHi < QL);

	asserta(m_Model != 0);
	const PSSM &P = m_Model->GetPSSM(GroupIndex, MotifIndex);
	const unsigned PL = P.GetColCount();
	const char *Q = m_QuerySeq.c_str();
	for (unsigned QPos = QLo; QPos + PL - 1 <= QHi; ++QPos)
		{
		float Score = P.GetScore(Q + QPos);
		if (Score > MinScore)
			{
			Hits.resize(Hits.size()+1);
			RPHit &NewHit = Hits.back();
			NewHit.Init(P, QPos, Score);
			if (m_Trace)
				{
				uint QLo2, QHi2;
				string QRow, PRow, ARow;
				GetAln(NewHit, QRow, PRow, ARow, QLo2, QHi2);

				char MotifLetter = "ABC"[MotifIndex];
				string GroupName;
				GetGroupName(GroupIndex, GroupName);

				Log("\n");
				Log("%s.%c", GroupName.c_str(), MotifLetter);
				Log(" (%.1f)", NewHit.m_Score);
				Log(" %u-%u", QLo2, QHi2);
				Log("\n");

				Log("  %s\n", QRow.c_str());
				Log("  %s\n", ARow.c_str());
				Log("  %s\n", PRow.c_str());
				}
			}
		}
	if (m_Trace)
		Log("  %u hits\n", SIZE(Hits));
	}

void RdRpSearcher::SearchMotif_TopHitOnly(uint GroupIndex, uint MotifIndex,
  uint QLo, uint QHi, float MinScore, RPHit &Hit)
	{
	if (m_Trace)
		{
		string GroupName;
		GetGroupName(GroupIndex, GroupName);
		Log("SearchMotif_TopHitOnly(%u=%s, %c, QLo=%u, QHi=%u, MinScore=%.1f)\n",
		  GroupIndex, GroupName.c_str(), MotLet, QLo, QHi, MinScore);
		}

	Hit.SetNoHit();

	asserta(QLo < QHi);
	const uint QL = SIZE(m_QuerySeq);
	asserta(QHi < QL);

	asserta(m_Model != 0);
	const PSSM &P = m_Model->GetPSSM(GroupIndex, MotifIndex);
	const unsigned PL = P.GetColCount();
	const char *Q = m_QuerySeq.c_str();
	for (unsigned QPos = QLo; QPos + PL < QHi; ++QPos)
		{
		float Score = P.GetScore(Q + QPos);
		if (Score >= MinScore && Score > Hit.m_Score)
			Hit.Init(P, QPos, Score);
		}
	if (m_Trace)
		{
		if (Hit.m_Score <= 0)
			Log("  SearchMotif_TopHitOnly: no hit\n");
		else
			{
			uint QLo2, QHi2;
			string QRow, PRow, ARow;
			GetAln(Hit, QRow, PRow, ARow, QLo2, QHi2);

			char MotifLetter = "ABC"[MotifIndex];
			string GroupName;
			GetGroupName(GroupIndex, GroupName);

			Log("\n");
			Log("%s.%c", GroupName.c_str(), MotifLetter);
			Log(" (%.1f)", Hit.m_Score);
			Log(" %u-%u", QLo, QHi);
			Log("\n");

			Log("  %s\n", QRow.c_str());
			Log("  %s\n", ARow.c_str());
			Log("  %s\n", PRow.c_str());
			}
		}
	}

void RdRpSearcher::CFilterGroup(uint GroupIndex)
	{
	uint QL = SIZE(m_QuerySeq);
	if (QL == 0)
		return;
	vector<RPHit> Hits;
	SearchMotif(GroupIndex, MOTIF_C, 0, QL-1, m_MinCScore, Hits);
	MergeCFilterHits(Hits);
	}

void RdRpSearcher::SearchGroup(uint GroupIndex)
	{
	uint QL = SIZE(m_QuerySeq);
	if (QL < AL+BL+CL+10)
		return;
	vector<RPHit> Hits_C;
	SearchMotif(GroupIndex, MOTIF_C, 0, QL-1, m_MinCScore, Hits_C);
	uint HitCount_C = SIZE(Hits_C);
	if (HitCount_C == 0)
		return;

	for (uint i = 0; i < HitCount_C; ++i)
		{
		RPHit Hit_A, Hit_B;
		SearchAB(GroupIndex, Hits_C[i]);
		}
	}

void RdRpSearcher::SearchAB(uint GroupIndex, const RPHit &Hit_C)
	{
	SearchAB_NotPermuted(GroupIndex, Hit_C);
	SearchAB_Permuted(GroupIndex, Hit_C);
	}

/***
            Insert           Insert
           <----->     <----------------->
   Alo   Ahi    Blo   Bhi               Clo   Chi
	|     |      |     |                 |     |
    AAAAAAA      BBBBBBB                 CCCCCCC
    <------------------------------------------>
                       PPlen
***/
bool RdRpSearcher::SearchAB_NotPermuted(uint GroupIndex,
  const RPHit &Hit_C)
	{
	int CLo = int(Hit_C.m_QPos);
	int CHi = int(Hit_C.GetQHi());

	if (m_Trace)
		{
		Log("\n");
		Log("SearchAB_NotPermuted()\n");
		}

	uint QL = SIZE(m_QuerySeq);
	uint PAL = GetPSSMLength(GroupIndex, MOTIF_A);

	int MinALo = CHi - m_MaxPPLength + 1;
	if (MinALo < 0)
		MinALo = 0;
	int MaxAHi = CHi - m_MinPPLength + PAL;
	if (MaxAHi < 0)
		{
		if (m_Trace)
			Log("  MaxAHi<0\n");
		return false;
		}
	if (MaxAHi >= int(QL))
		MaxAHi = int(QL) - 1;
	if (MinALo >= MaxAHi)
		{
		if (m_Trace)
			Log("  MinALo-%d <MaxAHi=%d\n", MinALo, MaxAHi);
		return false;
		}

	RPHit Hit_A;
	SearchMotif_TopHitOnly(GroupIndex, MOTIF_A, uint(MinALo), uint(MaxAHi),
	  m_MinPSSMScore, Hit_A);
	if (Hit_A.m_Score <= 0)
		return false;

	int AHi = int(Hit_A.GetQHi());

	int MinBLo = AHi + m_MinInsert;
	if (MinBLo >= int(QL))
		MinBLo = int(QL) - 1;
	int MaxBHi = CLo - m_MinInsert;
	if (MaxBHi >= int(QL))
		MaxBHi = int(QL) - 1;
	if (MinBLo >= MaxBHi)
		{
		if (m_Trace)
			Log("  MinBLo-%d >= MaxBHi=%d\n", MinBLo, MaxBHi);
		return false;
		}

	RPHit Hit_B;
	SearchMotif_TopHitOnly(GroupIndex, MOTIF_B, MinBLo, MaxBHi,
	  m_MinPSSMScore, Hit_B);
	if (Hit_B.m_Score <= 0)
		return false;

	OnPalmHit(GroupIndex, Hit_A, Hit_B, Hit_C, false);
	return true;
	}

/***
            Insert           Insert
           <----->     <----------------->
   Clo   Chi    Alo   Ahi               Blo   Bhi
	|     |      |     |                 |     |
    CCCCCCC      AAAAAAA                 BBBBBBB
    <------------------------------------------>
                       PPlen
***/
bool RdRpSearcher::SearchAB_Permuted(uint GroupIndex,
  const RPHit &Hit_C)
	{
	if (m_Trace)
		{
		Log("\n");
		Log("SearchAB_Permuted()\n");
		}

	int CLo = int(Hit_C.m_QPos);
	int CHi = int(Hit_C.GetQHi());

	uint QL = SIZE(m_QuerySeq);
	uint PBL = GetPSSMLength(GroupIndex, MOTIF_B);

	int MinBLo = CLo + m_MinPPLength - 1;
	if (MinBLo >= int(QL))
		MinBLo = int(QL) - 1;
	int MaxBHi = CLo + m_MaxPPLength + PBL;
	if (MaxBHi < 0)
		{
		if (m_Trace)
			Log(" MaxBHi<0\n");
		return false;
		}
	if (MaxBHi >= int(QL))
		MaxBHi = int(QL) - 1;
	if (MinBLo >= MaxBHi)
		{
		if (m_Trace)
			Log(" MinBLo=%d >= MaxBHi=%d\n", MinBLo, MaxBHi);
		return false;
		}

	RPHit Hit_B;
	SearchMotif_TopHitOnly(GroupIndex, MOTIF_B, uint(MinBLo), uint(MaxBHi),
	  m_MinPSSMScore, Hit_B);
	if (Hit_B.m_Score <= 0)
		return false;

	int BLo = int(Hit_B.m_QPos);

	int MinALo = CHi + m_MinInsert;
	if (MinALo >= int(QL))
		{
		if (m_Trace)
			Log(" MinALo >= QL\n");
		return false;
		}
	int MaxAHi = BLo - m_MinInsert;
	if (MaxAHi < 0)
		{
		if (m_Trace)
			Log("  MaxAHi < 0\n");
		return false;
		}
	if (MinALo >= MaxAHi)
		{
		if (m_Trace)
			Log("  MinALo=%d >= MaxAHi=%d\n", MinALo, MaxAHi);
		return false;
		}

	RPHit Hit_A;
	SearchMotif_TopHitOnly(GroupIndex, MOTIF_A, MinALo, MaxAHi,
	  m_MinPSSMScore, Hit_A);
	if (Hit_A.m_Score <= 0)
		return false;

	OnPalmHit(GroupIndex, Hit_A, Hit_B, Hit_C, true);
	return true;
	}

void RdRpSearcher::OnPalmHit(uint GroupIndex, const RPHit &Hit_A,
  const RPHit &Hit_B, const RPHit &Hit_C, bool Permuted)
	{
	if (m_Trace)
		Log("+++ OnPalmHit()\n");
	float TotalScore = Hit_A.m_Score + Hit_B.m_Score + Hit_C.m_Score;
	if (TotalScore > m_TopPalmHit.m_Score)
		{
		if (GroupIndex != m_SecondPalmHit.m_GroupIndex)
			m_SecondPalmHit = m_TopPalmHit;

		m_TopPalmHit.m_GroupIndex = GroupIndex;
		m_TopPalmHit.m_A = Hit_A;
		m_TopPalmHit.m_B = Hit_B;
		m_TopPalmHit.m_C = Hit_C;
		m_TopPalmHit.m_Score = TotalScore;
		m_TopPalmHit.m_Permuted = Permuted;
		return;
		}

	if (TotalScore > m_SecondPalmHit.m_Score
	  && GroupIndex != m_TopPalmHit.m_GroupIndex)
		{
		m_SecondPalmHit.m_GroupIndex = GroupIndex;
		m_SecondPalmHit.m_A = Hit_A;
		m_SecondPalmHit.m_B = Hit_B;
		m_SecondPalmHit.m_C = Hit_C;
		m_SecondPalmHit.m_Score = TotalScore;
		m_SecondPalmHit.m_Permuted = Permuted;
		}
	}
