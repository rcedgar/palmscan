#include "myutils.h"
#include "rdrpmodel.h"

uint RdRpModel::GetNextSeedPos(const string &QuerySeq, uint SeedPos)
	{
	const uint QL = SIZE(QuerySeq);
	const char *Q = QuerySeq.c_str();
	while (SeedPos + 2 < QL)
		{
		char c1 = Q[SeedPos];
		char c2 = Q[SeedPos+1];
		char c3 = Q[SeedPos+2];
#define	eq(x, y, z)	if (c1 == x && c2 == y && c3 == z) return SeedPos;
		eq('G', 'D', 'D')
		eq('S', 'D', 'D')
		eq('G', 'D', 'N')
#undef eq
		++SeedPos;
		}
	return UINT_MAX;
	}

static void CorrectResultSeeded(RPResult &Res, uint LoPos)
	{
	Res.m_StartPos += LoPos;
	}

void RdRpModel::SearchAA_Seeded(const string &QueryLabel, const string &QuerySeq)
	{
	RPResult BestResult;
	BestResult.Clear();
	const uint QL = SIZE(QuerySeq);
	SetResult_NoAaHit(QueryLabel, QL, m_Result);
	SetResult_NoAaHit(QueryLabel, QL, BestResult);
	uint SeedPos = 0;
	for (;;)
		{
		SeedPos = GetNextSeedPos(QuerySeq, SeedPos);
		if (SeedPos == UINT_MAX)
			break;

		int LoPos = int(SeedPos) - 350;
		int HiPos = int(SeedPos) + 159;
		if (LoPos < 0)
			LoPos = 0;
		if (HiPos >= int(QL))
			HiPos = QL- 1;
		int ChunkL = HiPos - LoPos + 1;
		if (ChunkL < int(MIN_SEG))
			break;
		string Chunk;
		for (int i = LoPos; i <= HiPos; ++i)
			Chunk += QuerySeq[i];

		SearchAA(QueryLabel, Chunk);
		if (m_Result.m_FinalScore > BestResult.m_FinalScore)
			{
			BestResult = m_Result;
			asserta(LoPos >= 0);
			CorrectResultSeeded(BestResult, uint(LoPos));
			}

		SeedPos += 50;
		}
	m_Result = BestResult;
	}
