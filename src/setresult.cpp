#include "myutils.h"
#include "rdrpmodel.h"
#include "heuristics.h"

static bool MatchMotif1(char q, char t)
	{
	if (q == t || t == 'x' || t == '.')
		return true;
	return false;
	}

static bool MatchMotif(const string &Motif, const string &Pattern)
	{
	const uint L = SIZE(Motif);
	asserta(SIZE(Pattern) == L);
	for (uint i = 0; i < L; ++i)
		if (!MatchMotif1(Motif[i], Pattern[i]))
			return false;
	return true;
	}

void RdRpModel::SetResult()
	{
	CheckThread();
	m_Result.Clear();
	m_Result.m_QueryLabel = m_QueryLabel;
	m_Result.m_Comments.clear();
	m_Result.m_Gene = "unclassified";
	m_Result.m_HiConf = false;
	m_Result.m_FinalScore = 0;

	GetAlnRows(m_Result.m_Aln);
	//GetFullAln(m_Result.m_FullAln);
	GetMotifsSeq(m_Result.m_MotifsSeq);
	GetMotifsSeq2(m_Result.m_MotifsSeq2);
	GetSuperMotif(m_Result.m_MotifsSeq, m_Result.m_SuperMotif);
	GetTrimmedSeq(m_Result.m_TrimmedSeq);

	bool AllowHigh = true;
	uint XCount = 0;
	for (uint i = 0; i < SIZE(m_Result.m_TrimmedSeq); ++i)
		{
		char c = m_Result.m_TrimmedSeq[i];
		if (c == '*')
			{
			m_Result.m_Comments += "stop-codon.";
			m_Result.m_FinalScore = 0.0;
			return;
			}

		if (c == 'X')
			++XCount;
		}
	if (XCount > m_MaxX)
		{
		m_Result.m_Comments += "too-many-Xs.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	m_Result.m_QL = SIZE(m_QuerySeq);
	if (m_Result.m_QL < MIN_SEG)
		{
		m_Result.m_Comments += "query-too-short.";
		m_Result.m_FinalScore = 0.0;
		return;
		}
	if (m_Result.m_QL < SHORT_QUERY)
		m_Result.m_Comments += "short-query.";

	uint TopGroup;
	float TopScore;
	GetTopGroup(TopGroup, TopScore);
	if (TopGroup == UINT_MAX)
		{
		m_Result.m_Comments += "motifs-not-found.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	GetGroupName(TopGroup, m_Result.m_TopGroupName);
	m_Result.m_PSSMTotalScore = TopScore;
	m_Result.m_FinalScore = TopScore;
	if (TopScore >= MIN_HICONF)
		m_Result.m_Comments += "high-PSSM-score.";

	if (StartsWith(m_Result.m_TopGroupName, "RT"))
		m_Result.m_Gene = "RT";
	else
		m_Result.m_Gene = "RdRP";

	const RPHit *ptrHit1 = 0;
	const RPHit *ptrHit2 = 0;
	const RPHit *ptrHit3 = 0;
	GetOrderedHits(m_Result.m_XXX, ptrHit1, ptrHit2, ptrHit3);
	if (m_Result.m_XXX.empty())
		{
		AllowHigh = false;
		m_Result.m_Comments += "motifs-out-of-order.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	bool Permuted = (m_Result.m_XXX != "ABC");
	if (Permuted && m_Result.m_Gene == "RT")
		{
		m_Result.m_Comments += "reject-permuted-RT.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	float CScore = GetCScore();
	if (CScore < MIN_C_SCORE)
		{
		AllowHigh = false;
		m_Result.m_Comments += "reject-low-C-score.";
		m_Result.m_FinalScore = 0;
		return;
		}

	float MinScore = GetMinScore();
	m_Result.m_PSSMMinScore = MinScore;
	if (Permuted && MinScore < MIN_PSSM_SCORE_PERMUTED)
		{
		AllowHigh = false;
		m_Result.m_Comments += "permuted-low-PSSM-score.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	if (MinScore < MIN_PSSM_SCORE)
		{
		AllowHigh = false;
		m_Result.m_Comments += "penalty-low-PSSM-score.";
		m_Result.m_FinalScore -= MIN_PSSM_PENALTY;
		}

	asserta(ptrHit1 != 0 && ptrHit2 != 0 && ptrHit3 != 0);
	m_Result.m_StartPos = ptrHit1->m_QPos;

	m_Result.m_Dist12 = GetDist(ptrHit1, ptrHit2);
	m_Result.m_Dist23 = GetDist(ptrHit2, ptrHit3);
	if (m_Result.m_Dist12 <= 0 || m_Result.m_Dist23 <= 0)
		{
		AllowHigh = false;
		m_Result.m_Comments += "reject-overlapping-motifs.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	uint Span = GetSpan(ptrHit1, ptrHit3);
	m_Result.m_SegLength = Span;
	if (Span < MIN_SEG || Span > MAX_SEG)
		{
		AllowHigh = false;
		m_Result.m_Comments += "bad-segment-length.";
		m_Result.m_FinalScore = 0.0;
		return;
		}

	if (m_Result.m_Gene == "RdRP")
		{
		const string &Super = m_Result.m_SuperMotif;
		if (m_Result.m_XXX == "ABC" && !Super.empty())
			{
			const string &GDD = Super.substr(3, 3);
			if (GDD != "GDD" && GDD != "SDD" && GDD != "GDN")
				{
				m_Result.m_Comments += "penalty-missing-GDD/SDD/GDN.";
				m_Result.m_FinalScore -= PENALTY_NOT_GDD_SDD_GDN;
				}
			else if (MatchMotif(Super, "DDGGDD"))
				{
				m_Result.m_Comments += "reward-DDGGDD.";
				m_Result.m_FinalScore += REWARD_DDGGDD;
				}
			else if (MatchMotif(Super, "DDGSDD"))
				{
				m_Result.m_Comments += "reward-DDGSDD.";
				m_Result.m_FinalScore += REWARD_DDGSDD;
				}
			else if (MatchMotif(Super, "DNMSDD"))
				{
				m_Result.m_Comments += "reward-DNMSDD.";
				m_Result.m_FinalScore += REWARD_DNMSDD;
				}
			else if (MatchMotif(Super, "DNxGDN"))
				{
				m_Result.m_Comments += "reward-DNxGDN.";
				m_Result.m_FinalScore += REWARD_DNxGDN;
				}
			}

		if (Span < MIN_SEG1 || Span > MAX_SEG1)
			{
			m_Result.m_Comments += "segment-length-penalty1.";
			m_Result.m_FinalScore -= SEG_LENGTH_PENALTY1;
			}
		else if (Span < MIN_SEG2 || Span > MAX_SEG2)
			{
			m_Result.m_Comments += "segment-length-penalty2.";
			m_Result.m_FinalScore -= SEG_LENGTH_PENALTY2;
			}
		else
			m_Result.m_Comments += "good-segment-length.";

		if (m_Result.m_XXX == "CAB")
			{
			m_Result.m_Comments += "permuted-penalty.";
			m_Result.m_FinalScore -= PERMUTED_PENALTY;
			}
		else
			{
			int AB = m_Result.m_Dist12;
			int BC = m_Result.m_Dist23;
			if (AB < MIN_AB1 || AB > MAX_AB1)
				{
				m_Result.m_Comments += "AB-dist-penalty.";
				m_Result.m_FinalScore -= AB_PENALTY1;
				}
			if (BC < MIN_BC1 || BC > MAX_BC1)
				{
				m_Result.m_Comments += "BC-dist-penalty.";
				m_Result.m_FinalScore -= BC_PENALTY1;
				}
			}
		} // end if RdRP

	if (m_Result.m_FinalScore < 0)
		{
		m_Result.m_FinalScore = 0;
		m_Result.m_Gene = "unclassified";
		m_Result.m_Comments += "negative-score.";
		return;
		}

	if (!AllowHigh && m_Result.m_FinalScore >= MIN_HICONF)
		m_Result.m_FinalScore = MIN_LOCONF;
	m_Result.m_HiConf = (m_Result.m_FinalScore >= MIN_HICONF);
	}

void RPResult::ToFEVStr(string &s) const
	{
	Ps(s, "score=%.1f", m_FinalScore);
	Psa(s, "\tquery=%s", m_QueryLabel.c_str());
	Psa(s, "\tgene=%s", m_Gene.c_str());
	if (m_Gene == "RdRP" || m_Gene == "RT")
		Psa(s, "\tconfidence=%s", m_HiConf ? "high" : "low");
	Psa(s, "\tcomments=%s", m_Comments.c_str());
	if (!m_XXX.empty())
		Psa(s, "\torder=%s", m_XXX.c_str());
	if (m_QLnt != 0 && m_QLnt != UINT_MAX)
		Psa(s, "\tqlen=%u", m_QLnt);
	else if (m_QL != 0 && m_QL != UINT_MAX)
		Psa(s, "\tqlen=%u", m_QL);
	if (m_StartPosNt != UINT_MAX)
		{
		Psa(s, "\txlat=yes");
		Psa(s, "\tpp_start=%u", m_StartPosNt + 1);
		Psa(s, "\tpp_end=%u", m_StartPosNt + 3*m_SegLength);
		Psa(s, "\tpp_length=%u", 3*m_SegLength);
		}
	else if (m_StartPos != UINT_MAX)
		{
		Psa(s, "\tpp_start=%u", m_StartPos + 1);
		Psa(s, "\tpp_end=%u", m_StartPos + m_SegLength);
		Psa(s, "\tpp_length=%u", m_SegLength);
		}
	if (m_Frame != 0 && m_Frame != -999)
		Psa(s, "\tframe=%+d", m_Frame);
	if (m_PSSMTotalScore > -10)
		Psa(s, "\tpssm_total_score=%.1f", m_PSSMTotalScore);
	if (m_PSSMMinScore > -10)
		Psa(s, "\tpssm_min_score=%.1f", m_PSSMMinScore);
	if (m_XXX == "ABC" || m_XXX == "CAB")
		{
		Psa(s, "\torder=%s", m_XXX.c_str());
		Psa(s, "\tv1_length=%d", m_Dist12);
		Psa(s, "\tv2_length=%d", m_Dist23);
		}
	if (!m_TopGroupName.empty())
		Psa(s, "\tgroup=%s", m_TopGroupName.c_str());
	if (!m_SuperMotif.empty())
		Psa(s, "\tsuper=%s", m_SuperMotif.c_str());
	if (!m_MotifsSeq.empty())
		Psa(s, "\tmotifs=%s", m_MotifsSeq.c_str());
	}

void RdRpModel::SetResult_NoNtHit(const string &QueryLabel, uint QLnt,
  RPResult &Result)
	{
	Result.Clear();
	Result.m_QueryLabel = QueryLabel;
	Result.m_QL = 0;
	Result.m_QLnt = QLnt;
	Result.m_Comments.clear();
	Result.m_Gene = "unclassified";
	Result.m_HiConf = false;
	Result.m_Comments = "six-frame-none-found.";
	Result.m_FinalScore = 0;
	}

void RdRpModel::SetResult_NoAaHit(const string &QueryLabel, uint QLaa,
  RPResult &Result)
	{
	Result.Clear();
	Result.m_QueryLabel = QueryLabel;
	Result.m_QL = QLaa;
	Result.m_QLnt = UINT_MAX;
	Result.m_Comments.clear();
	Result.m_Gene = "unclassified";
	Result.m_HiConf = false;
	Result.m_Comments = "aa-none-found.";
	Result.m_FinalScore = 0;
	}
