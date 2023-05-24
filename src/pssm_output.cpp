#include "myutils.h"
#include "rdrpsearcher.h"
#include "outputfiles.h"
#include "abcxyz.h"

extern vector<string> g_ExcludeNames;

void RdRpSearcher::WriteOutput() const
	{
	if (!g_ExcludeNames.empty())
		{
		if (m_TopPalmHit.m_Score > 0)
			{
			uint GroupIndex = m_TopPalmHit.m_GroupIndex;
			string GroupName;
			GetGroupName(GroupIndex, GroupName);
			for (uint i = 0; i < SIZE(g_ExcludeNames); ++i)
				{
				if (GroupName == g_ExcludeNames[i])
					return;
				}
			}
		}

#pragma omp critical
	{
	WriteReport(g_freport_pssms);
	WriteTsv(g_ftsv);
	WriteFev(g_ffev);
	WritePalmprintFasta(g_ffasta_abc, g_ffasta_cab);
	WritePalmprintFasta(g_ffasta, g_ffasta);
	WriteFlankFasta(g_fleft_abc, g_fleft_cab, true);
	WriteFlankFasta(g_fright_abc, g_fright_cab, false);
	WriteMotifs(g_fmotifs, g_fmotifs);
	WriteMotifs(g_fmotifs_abc, g_fmotifs_cab);
	WriteCore(g_fcore, g_fcore);
	bool Ok = WriteCore(g_fcore_abc, g_fcore_cab);
	if (!Ok)
		SeqToFasta(g_fcore_notmatched, m_QueryLabel.c_str(),
			m_QuerySeq.c_str(), SIZE(m_QuerySeq));
	}
	}

void RdRpSearcher::WritePalmprintFasta(FILE *fABC, FILE *fCAB) const
	{
	if (m_TopPalmHit.m_Score <= 0)
		return;

	uint PosA = GetMotifPos(0);
	uint PosB = GetMotifPos(1);
	uint PosC = GetMotifPos(2);
	if (PosA == UINT_MAX || PosB == UINT_MAX || PosC == UINT_MAX)
		return;

	bool Perm = m_TopPalmHit.m_Permuted;

	uint PPLo = UINT_MAX;
	uint PPHi = 0;
	if (Perm)
		{
		PPLo = PosC;
		PPHi = PosB + BL - 1;
		}
	else
		{
		PPLo = PosA;
		PPHi = PosC + CL - 1;
		}
	asserta(PPHi > PPLo);
	const uint PPL = PPHi - PPLo + 1;

	string SeqA;
	string SeqB;
	string SeqC;
	GetMotifSeq(0, SeqA);
	GetMotifSeq(1, SeqB);
	GetMotifSeq(2, SeqC);
	if (SeqA == "" || SeqB == "" || SeqC == "")
		return;

	uint GroupIndex = m_TopPalmHit.m_GroupIndex;
	string GroupName;
	GetGroupName(GroupIndex, GroupName);

	string Label = m_QueryLabel;

	Psa(Label, " A:%u:%s", PosA - PPLo + 1, SeqA.c_str());
	Psa(Label, " B:%u:%s", PosB - PPLo + 1, SeqB.c_str());
	Psa(Label, " C:%u:%s", PosC - PPLo + 1, SeqC.c_str());
	Psa(Label, " PSSM:%s", GroupName.c_str());

	const char *PPSeq = m_QuerySeq.c_str() + PPLo;
	SeqToFasta(Perm ? fCAB : fABC, Label.c_str(), PPSeq, PPL);
	}

void RdRpSearcher::WriteReport(FILE *f) const
	{
	if (f == 0)
		return;

	if (m_TopPalmHit.m_Score <= 0)
		return;

	uint GroupIndex = m_TopPalmHit.m_GroupIndex;
	bool Permuted = m_TopPalmHit.m_Permuted;
	float Score = m_TopPalmHit.m_Score;

	string GroupName;
	GetGroupName(GroupIndex, GroupName);

	vector<string> Rows;
	GetAlnRows(Rows);

	fprintf(f, "\n");
	fprintf(f, ">%s\n", m_QueryLabel.c_str());
	for (uint i = 0; i < SIZE(Rows); ++i)
		fprintf(f, "%s\n", Rows[i].c_str());

	fprintf(f, "Score %.1f", Score);
	if (Permuted)
		fprintf(f, " permuted");
	fprintf(f, ", %s", GroupName.c_str());

	float SecondScore = m_SecondPalmHit.m_Score;
	if (SecondScore > 0)
		{
		string SecondGroupName;
		uint SecondGroupIndex = m_SecondPalmHit.m_GroupIndex;
		GetGroupName(SecondGroupIndex, SecondGroupName);
		float Diff = Score - SecondScore;
		fprintf(f, " (%s +%.1f)", SecondGroupName.c_str(), Diff);
		}
	fprintf(f, "\n");
	}

void RdRpSearcher::WriteMotifs(FILE* fABC, FILE *fCAB) const
	{
	if (m_TopPalmHit.m_Score <= 0)
		return;

	string SeqA;
	string SeqB;
	string SeqC;
	GetMotifSeq(0, SeqA);
	GetMotifSeq(1, SeqB);
	GetMotifSeq(2, SeqC);
	if (SeqA == "" || SeqB == "" || SeqC == "")
		return;
	FILE *f = (m_TopPalmHit.m_Permuted ? fCAB : fABC);
	if (f == 0)
		return;

	string xxx = "";
	if (optset_xxx)
		xxx = opt_xxx;

	string Motifs = SeqA + xxx + SeqB + xxx + SeqC;
	SeqToFasta(f, m_QueryLabel.c_str(), Motifs.c_str(), SIZE(Motifs));
	}

void RdRpSearcher::WriteFlankFasta(FILE *fABC, FILE *fCAB, bool Left) const
	{
	if (m_TopPalmHit.m_Score < m_MinScore_Core)
		return;

	uint PosA = GetMotifPos(0);
	uint PosB = GetMotifPos(1);
	uint PosC = GetMotifPos(2);
	if (PosA == UINT_MAX || PosB == UINT_MAX || PosC == UINT_MAX)
		return;

	uint QL = SIZE(m_QuerySeq);
	bool Perm = m_TopPalmHit.m_Permuted;
	uint FlankLo = UINT_MAX;
	uint FlankHi = UINT_MAX;
	uint FlankLen = UINT_MAX;
	if (Left)
		{
		uint FlankHi1 = (Perm ? PosC : PosA);
		if (FlankHi1 <= 1)
			return;
		FlankHi = FlankHi1 - 1;
		FlankLo = 0;
		FlankLen = FlankHi - FlankLo + 1;
		FlankLen = min(FlankLen, m_LeftFlank_Core);
		FlankLo = FlankHi + 1 - FlankLen;
		}
	else
		{
		FlankLo = (Perm ? PosB + BL : PosC + CL);
		if (FlankLo >= QL - 1)
			return;
		FlankLen = QL - FlankLo;
		FlankLen = min(FlankLen, m_RightFlank_Core);
		FlankHi = FlankLo + FlankLen - 1;
		asserta(FlankHi < QL);
		}

	if (FlankLen < m_MinFlankLen)
		return;

	uint GroupIndex = m_TopPalmHit.m_GroupIndex;
	string GroupName;
	GetGroupName(GroupIndex, GroupName);

	string Label = m_QueryLabel;
	Label += " group=";
	Label += GroupName;

	const char *Start = m_QuerySeq.c_str() + FlankLo;
	SeqToFasta(Perm ? fCAB : fABC, Label.c_str(), Start, FlankLen);
	}

bool RdRpSearcher::WriteCore(FILE *fABC, FILE *fCAB) const
	{
	if (m_TopPalmHit.m_Score < m_MinScore_Core)
		return false;

	uint PosA = GetMotifPos(0);
	uint PosB = GetMotifPos(1);
	uint PosC = GetMotifPos(2);
	if (PosA == UINT_MAX || PosB == UINT_MAX || PosC == UINT_MAX)
		return false;

	uint QL = SIZE(m_QuerySeq);
	uint LF = UINT_MAX;
	uint RF = UINT_MAX;
	bool Perm = m_TopPalmHit.m_Permuted;
	if (Perm)
		{
		LF = PosC;
		asserta(PosA < QL);
		RF = QL - PosA - 1;
		}
	else
		{
		LF = PosA;
		asserta(PosC < QL);
		RF = QL - PosC - 1;
		}
	LF = min(LF, m_LeftFlank_Core);
	RF = min(RF, m_RightFlank_Core);
	if (LF < m_MinFlankLen)
		return false;
	if (RF < m_MinFlankLen)
		return false;

	uint CoreLo = UINT_MAX;
	uint CoreHi = UINT_MAX;
	if (Perm)
		{
		asserta(LF <= PosC);
		CoreLo = PosC - LF;
		CoreHi = PosA + RF;
		}
	else
		{
		asserta(LF <= PosA);
		CoreLo = PosA - LF;
		CoreHi = PosC + RF;
		}
	asserta(CoreHi < QL);
	uint CoreL = CoreHi - CoreLo + 1;

	string SeqA;
	string SeqB;
	string SeqC;
	GetMotifSeq(0, SeqA);
	GetMotifSeq(1, SeqB);
	GetMotifSeq(2, SeqC);

	uint GroupIndex = m_TopPalmHit.m_GroupIndex;
	string GroupName;
	GetGroupName(GroupIndex, GroupName);

	string NewLabel = m_QueryLabel;
	Psa(NewLabel, " A:%u:%s", PosA+1-CoreLo, SeqA.c_str());
	Psa(NewLabel, " B:%u:%s", PosB+1-CoreLo, SeqB.c_str());
	Psa(NewLabel, " C:%u:%s", PosC+1-CoreLo, SeqC.c_str());
	Psa(NewLabel, " PSSM:%s", GroupName.c_str());

	const char *CoreStart = m_QuerySeq.c_str() + CoreLo;
	SeqToFasta(Perm ? fCAB : fABC,
	  NewLabel.c_str(), CoreStart, CoreL);
	return true;
	}

void RdRpSearcher::GetCoreSeq(string &CoreSeq) const
	{
	CoreSeq.clear();
	uint PosA = GetMotifPos(0);
	uint PosB = GetMotifPos(1);
	uint PosC = GetMotifPos(2);
	if (PosA == UINT_MAX || PosB == UINT_MAX || PosC == UINT_MAX)
		return;

	uint QL = SIZE(m_QuerySeq);
	uint LF = UINT_MAX;
	uint RF = UINT_MAX;
	bool Perm = m_TopPalmHit.m_Permuted;
	if (Perm)
		{
		LF = PosC;
		asserta(PosA < QL);
		RF = QL - PosA - 1;
		}
	else
		{
		LF = PosA;
		asserta(PosC < QL);
		RF = QL - PosC - 1;
		}
	LF = min(LF, m_LeftFlank_Core);
	RF = min(RF, m_RightFlank_Core);
	if (LF < m_MinFlankLen)
		return;
	if (RF < m_MinFlankLen)
		return;

	uint CoreLo = UINT_MAX;
	uint CoreHi = UINT_MAX;
	if (Perm)
		{
		asserta(LF <= PosC);
		CoreLo = PosC - LF;
		CoreHi = PosA + RF;
		}
	else
		{
		asserta(LF <= PosA);
		CoreLo = PosA - LF;
		CoreHi = PosC + RF;
		}
	asserta(CoreHi < QL);
	uint CoreL = CoreHi - CoreLo + 1;

	const char *CoreStart = m_QuerySeq.c_str() + CoreLo;
	for (uint i = 0; i < CoreL; ++i)
		CoreSeq += CoreStart[i];
	}

void RdRpSearcher::WriteTsv(FILE *f) const
	{
	if (f == 0)
		return;

	static bool HdrDone = false;
	if (!HdrDone)
		{
		HdrDone = true;

		fprintf(f, "Label");
		fprintf(f, "\tScore");
		fprintf(f, "\tGroup");
		fprintf(f, "\tGroup2");
		fprintf(f, "\tDiff2");
		fprintf(f, "\tABC");
		fprintf(f, "\tQL");
		fprintf(f, "\tLo");
		fprintf(f, "\tHi");
		fprintf(f, "\tPPL");
		fprintf(f, "\tSuff");
		fprintf(f, "\tPosA");
		fprintf(f, "\tSeqA");
		fprintf(f, "\tPosB");
		fprintf(f, "\tSeqB");
		fprintf(f, "\tPosC");
		fprintf(f, "\tSeqC");
		fprintf(f, "\n");
		}

	if (m_TopPalmHit.m_Score <= 0)
		return;

	string QueryLabel = m_QueryLabel;
	float Score = 0;
	string GroupName = ".";
	string SecondGroupName = ".";
	float Diff2 = 0;
	const char *ABC = ".";
	uint QL = SIZE(m_QuerySeq);
	uint Lo = 0;
	uint Hi = 0;
	uint PPL = 0;
	uint Suff = 0;

	if (m_TopPalmHit.m_Score > 0)
		{
		Score = m_TopPalmHit.m_Score;
		uint GroupIndex = m_TopPalmHit.m_GroupIndex;
		GetGroupName(GroupIndex, GroupName);
		if (m_TopPalmHit.m_Permuted)
			ABC = "CAB";
		else
			ABC = "ABC";
		GetSpan(Lo, Hi);
		++Lo;
		++Hi;
		assert(Hi <= QL);
		PPL = Hi - Lo + 1;
		Suff = QL - Hi;
		}

	if (m_SecondPalmHit.m_Score > 0)
		{
		float SecondScore = m_SecondPalmHit.m_Score;
		uint SecondGroupIndex = m_SecondPalmHit.m_GroupIndex;
		GetGroupName(SecondGroupIndex, SecondGroupName);
		Diff2 = Score - SecondScore;
		}

	uint PosA = GetMotifPos(0);
	uint PosB = GetMotifPos(1);
	uint PosC = GetMotifPos(2);

	string SeqA = ".";
	string SeqB = ".";
	string SeqC = ".";

	if (PosA == UINT_MAX)
		PosA = 0;
	else
		{
		GetMotifSeq(0, SeqA);
		++PosA;
		}
	if (PosB == UINT_MAX)
		PosB = 0;
	else
		{
		GetMotifSeq(1, SeqB);
		++PosB;
		}
	if (PosC == UINT_MAX)
		PosC = 0;
	else
		{
		GetMotifSeq(2, SeqC);
		++PosC;
		}

	fprintf(f, "%s", QueryLabel.c_str());
	fprintf(f, "\t%.1f", Score);
	fprintf(f, "\t%s", GroupName.c_str());
	fprintf(f, "\t%s", SecondGroupName.c_str());
	fprintf(f, "\t%+.1f", Diff2);
	fprintf(f, "\t%s", ABC);
	fprintf(f, "\t%u", QL);
	fprintf(f, "\t%u", Lo);
	fprintf(f, "\t%u", Hi);
	fprintf(f, "\t%u", PPL);
	fprintf(f, "\t%u", Suff);
	fprintf(f, "\t%u", PosA);
	fprintf(f, "\t%s", SeqA.c_str());
	fprintf(f, "\t%u", PosB);
	fprintf(f, "\t%s", SeqB.c_str());
	fprintf(f, "\t%u", PosC);
	fprintf(f, "\t%s", SeqC.c_str());
	fprintf(f, "\n");
	}

void RdRpSearcher::WriteFev(FILE *f) const
	{
	if (f == 0)
		return;

	if (m_TopPalmHit.m_Score <= 0)
		return;

	string QueryLabel = m_QueryLabel;
	float Score = 0;
	string GroupName = ".";
	string SecondGroupName = ".";
	float Diff2 = 0;
	const char *ABC = ".";
	uint QL = SIZE(m_QuerySeq);
	uint Lo = 0;
	uint Hi = 0;
	uint PPL = 0;
	uint Suff = 0;

	if (m_TopPalmHit.m_Score > 0)
		{
		Score = m_TopPalmHit.m_Score;
		uint GroupIndex = m_TopPalmHit.m_GroupIndex;
		GetGroupName(GroupIndex, GroupName);
		if (m_TopPalmHit.m_Permuted)
			ABC = "CAB";
		else
			ABC = "ABC";
		GetSpan(Lo, Hi);
		++Lo;
		++Hi;
		assert(Hi <= QL);
		PPL = Hi - Lo + 1;
		Suff = QL - Hi;
		}

	if (m_SecondPalmHit.m_Score > 0)
		{
		float SecondScore = m_SecondPalmHit.m_Score;
		uint SecondGroupIndex = m_SecondPalmHit.m_GroupIndex;
		GetGroupName(SecondGroupIndex, SecondGroupName);
		Diff2 = Score - SecondScore;
		}

	uint PosA = GetMotifPos(0);
	uint PosB = GetMotifPos(1);
	uint PosC = GetMotifPos(2);

	string SeqA;
	string SeqB;
	string SeqC;
	GetMotifSeq(0, SeqA);
	GetMotifSeq(1, SeqB);
	GetMotifSeq(2, SeqC);

	fprintf(f, "%s", QueryLabel.c_str());
	fprintf(f, "\tpssm_score=%.1f", Score);
	fprintf(f, "\tpssm_group=%s", GroupName.c_str());
	fprintf(f, "\tpssm_group2=%s", SecondGroupName.c_str());
	fprintf(f, "\tpssm_diff2=%+.1f", Diff2);
	fprintf(f, "\tpssm_ABC=%s", ABC);
	if (SeqA != "")
		{
		fprintf(f, "\tpssm_seqA=%s", SeqA.c_str());
		fprintf(f, "\tpssm_posA=%u", PosA+1);
		}
	if (SeqB != "")
		{
		fprintf(f, "\tpssm_seqB=%s", SeqB.c_str());
		fprintf(f, "\tpssm_posB=%u", PosB+1);
		}
	if (SeqC != "")
		{
		fprintf(f, "\tpssm_seqC=%s", SeqC.c_str());
		fprintf(f, "\tpssm_posC=%u", PosC+1);
		}
	fprintf(f, "\n");
	}
