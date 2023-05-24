#include "myutils.h"
#include "ppcsearcher.h"
#include "searchparams.h"
#include "sort.h"
#include "outputfiles.h"
#include "abcxyz.h"

uint PPCSearcher::m_MaxKeep = 10;

void PPCSearcher::SetQuery(PDBChain &Q)
	{
	asserta(Q.CheckMotifCoords());

	m_Query = &Q;
	m_PA.SetQuery(Q);
	}

void PPCSearcher::SetRefs(const vector<PDBChain *> &Refs)
	{
	const uint N = SIZE(Refs);
	for (uint i = 0; i < N; ++i)
		asserta(Refs[i]->CheckPPCMotifCoords());
	m_Refs = &Refs;
	}

bool PPCSearcher::Search(PDBChain &Q)
	{
	ClearSearch();
	SetQuery(Q);
	asserta(m_Refs != 0);

	const vector<PDBChain *> &Rs = *m_Refs;
	const uint RefN = SIZE(Rs);

	vector<uint> RefIndexesPass1;
	for (uint RefIndex = 0; RefIndex < RefN; ++RefIndex)
		{
		const PDBChain &R = *Rs[RefIndex];

		if (opt_self && Q.m_Label == R.m_Label)
			continue;

		m_PA.SetRef(R);
		double MotifRMSD = m_PA.GetMotifRMSD();
		if (MotifRMSD <= MaxMotifRMSD)
			{
			m_RMSDs.push_back(MotifRMSD);
			RefIndexesPass1.push_back(RefIndex);
			}
		}
	uint HitCount = SIZE(m_RMSDs);
	if (HitCount == 0)
		return false;

	Q.SetSS();

	vector<uint> RMSDOrder(HitCount);
	QuickSortOrder(m_RMSDs.data(), HitCount, RMSDOrder.data());

	TSHit Hit;
	Hit.m_Query = &Q;
	Hit.m_QPosA = Q.m_MotifPosVec[A];
	Hit.m_QPosB = Q.m_MotifPosVec[B];
	Hit.m_QPosC = Q.m_MotifPosVec[C];

	uint N = min(m_MaxKeep, HitCount);
	m_Hits.clear();
	for (uint k = 0; k < N; ++k)
		{
		uint i = RMSDOrder[k];
		double RMSD = m_RMSDs[i];
		uint RefIndex = RefIndexesPass1[i];
		const PDBChain &R = *Rs[RefIndex];

		Hit.m_Ref = &R;
		Hit.m_MotifRMSD2 = RMSD*RMSD;
		Hit.SetSketch();
		m_Scores.push_back(RMSD);
		m_Hits.push_back(Hit);
		m_RefIndexes.push_back(RefIndex);
		}

	m_HitOrder.resize(N);
	QuickSortOrder(m_Scores.data(), N, m_HitOrder.data());

	uint TopHitIndex = m_HitOrder[0];
	asserta(TopHitIndex < SIZE(m_Hits));
	m_TopHit = &m_Hits[TopHitIndex];

	if (m_DoClassify)
		Classify();

	return true;
	}

const TSHit *PPCSearcher::GetTopHit() const
	{
	if (m_HitOrder.empty())
		return 0;
	uint TopHitIndex = m_HitOrder[0];
	asserta(TopHitIndex < SIZE(m_Hits));
	const TSHit *TopHit = &m_Hits[TopHitIndex];
	return TopHit;
	}

void PPCSearcher::WriteOutput() const
	{
#pragma omp critical
	{
	const uint HitCount = SIZE(m_Hits);
	asserta(SIZE(m_HitOrder) == HitCount);
	asserta(SIZE(m_Scores) == HitCount);
	asserta(SIZE(m_RefIndexes) == HitCount);

	for (uint k = 0; k < HitCount; ++k)
		{
		uint i = m_HitOrder[k];
		const TSHit &Hit = m_Hits[i];
		Hit.WriteTsv(g_ftsv);
		}

	const TSHit *TopHit = GetTopHit();
	if (TopHit != 0)
		TopHit->WriteReport(g_freport_3d);

	if (m_DoClassify)
		WriteClassify(g_fclassify_tsv);
	}
	}

void PPCSearcher::Classify()
	{
	const uint HitCount = SIZE(m_Hits);
	asserta(SIZE(m_HitOrder) == HitCount);
	asserta(SIZE(m_Scores) == HitCount);
	asserta(SIZE(m_RefIndexes) == HitCount);

	m_Classified = true;
	m_ClassifiedAsRdRp = false;
	m_TopScore = 0;
	m_ScoreDiff = 0;
	if (HitCount == 0)
		return;

	m_TopScore = m_TopHit->m_Score;
	m_ScoreDiff = m_TopScore;
	bool TopRdRp = TopHitIsRdRp();
	asserta(m_TopHit != 0);

	m_TopOtherHit = 0;
	uint TopOtherHitIndex = GetTopClassifiedHitIndex(!TopRdRp);
	if (TopOtherHitIndex != UINT_MAX)
		{
		asserta(TopOtherHitIndex < SIZE(m_Hits));
		m_TopOtherHit = &m_Hits[TopOtherHitIndex];
		}

	if (m_TopOtherHit != 0)
		m_ScoreDiff -= m_TopOtherHit->m_Score;
	m_ClassifiedAsRdRp = (TopRdRp && m_TopScore >= MinScoreClassify);
	}

void PPCSearcher::InitClassify()
	{
	m_DoClassify = true;
	m_RefIsRdRps.clear();
	const vector<PDBChain *> &Rs = *m_Refs;
	const uint N = SIZE(Rs);
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &R = *Rs[i];
		bool IsRdRp = GetIsRdRpFromLabel(R.m_Label);
		m_RefIsRdRps.push_back(IsRdRp);
		}
	}

bool PPCSearcher::GetIsRdRpFromLabel(const string &Label)
	{
	bool IsRdRp = StartsWith(Label, "rdrp.");
	return IsRdRp;
	}

const uint PPCSearcher::GetTopClassifiedHitIndex(bool IsRdRp) const
	{
	asserta(!m_RefIsRdRps.empty());
	const uint HitCount = SIZE(m_Hits);
	asserta(SIZE(m_HitOrder) == HitCount);
	asserta(SIZE(m_Scores) == HitCount);
	asserta(SIZE(m_RefIndexes) == HitCount);

	for (uint k = 0; k < HitCount; ++k)
		{
		uint HitIndex = m_HitOrder[k];
		double Score = m_Scores[HitIndex];
		uint RefIndex = m_RefIndexes[HitIndex];
		asserta(RefIndex < SIZE(m_RefIsRdRps));
		bool RefIsRdRp = m_RefIsRdRps[RefIndex];
		if (RefIsRdRp == IsRdRp)
			return HitIndex;
		}
	return UINT_MAX;
	}

bool PPCSearcher::TopHitIsRdRp() const
	{
	asserta(!m_RefIsRdRps.empty());
	asserta(!m_HitOrder.empty());
	uint HitIndex = m_HitOrder[0];
	asserta(HitIndex < SIZE(m_RefIndexes));
	uint RefIndex = m_RefIndexes[HitIndex];
	asserta(RefIndex < SIZE(m_RefIsRdRps));
	bool IsRdRp = m_RefIsRdRps[RefIndex];
	return IsRdRp;
	}

void PPCSearcher::WriteClassify(FILE *f) const
	{
	if (f == 0)
		return;

	if (m_Hits.empty())
		return;

	asserta(m_TopHit != 0);
	fprintf(f, "%.3f", m_TopScore);
	fprintf(f, "\t%s", m_ClassifiedAsRdRp ? "rdrp" : "other");
	fprintf(f, "\t%+.3f", m_ScoreDiff);
	fprintf(f, "\t%s", m_Query->m_Label.c_str());
	fprintf(f, "\t%s", m_TopHit->m_Ref->m_Label.c_str());
	if (m_TopOtherHit == 0)
		fprintf(f, "\t.\t.");
	else
		{
		fprintf(f, "\t%.3f\t%s",
		  m_TopOtherHit->m_Score,
		  m_TopOtherHit->m_Ref->m_Label.c_str());
		}
	fputc('\n', f);
	}
