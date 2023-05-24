#include "myutils.h"
#include "rdrpsearcher.h"
#include "sort.h"

void RdRpSearcher::Init(const RdRpModel &Model)
	{
	Clear();

	m_Model = &Model;
	}

const PSSM &RdRpSearcher::GetPSSM(uint GroupIndex, uint MotifIndex) const
	{
	asserta(m_Model != 0);
	return m_Model->GetPSSM(GroupIndex, MotifIndex);
	}

uint RdRpSearcher::GetPSSMLength(uint GroupIndex, uint MotifIndex) const
	{
	asserta(m_Model != 0);
	uint L = m_Model->GetPSSMLength(GroupIndex, MotifIndex);
	return L;
	}

void RdRpSearcher::GetMotifSeq(uint MotifIndex, string &Seq) const
	{
	Seq.clear();
	const RPHit *Hit = &m_TopPalmHit.GetHit(MotifIndex);
	if (Hit == 0)
		return;
	uint QPos = Hit->m_QPos;
	uint QLo = QPos + 1;
	uint QHi = Hit->GetQHi() + 1;
	const PSSM *P = Hit->m_PSSM;
	if (P == 0)
		return;
	uint PL = P->m_ColCount;
	uint PPos = 0;
	uint QL = SIZE(m_QuerySeq);
	for (uint PPos = 0; PPos < PL; ++PPos)
		{
		asserta(QPos+PPos < QL);
		char q = m_QuerySeq[QPos+PPos];
		Seq += q;
		}
	}

void RdRpSearcher::GetAln(const RPHit &Hit,
  string &QRow, string &PRow, string &AnnotRow,
  uint &QLo, uint &QHi) const
	{
	QRow.clear();
	PRow.clear();
	AnnotRow.clear();

	const uint QL = SIZE(m_QuerySeq);
	if (Hit.m_Score <= 0)
		return;
	asserta(Hit.m_PSSM != 0);
	const PSSM &P = *Hit.m_PSSM;
	uint PL = P.m_ColCount;

	uint QPos = Hit.m_QPos;
	QLo = QPos + 1;
	QHi = Hit.GetQHi() + 1;

	uint PPos = 0;
	for (uint PPos = 0; PPos < PL; ++PPos)
		{
		char q = m_QuerySeq[QPos+PPos];
		float Score = P.GetScore1(PPos, q);
		char p = P.GetConsChar(PPos);

		if (Score >= 1)
			AnnotRow += '|';
		else if (Score >= 0.5)
			AnnotRow += '+';
		else if (Score > 0)
			AnnotRow += ".";
		else
			AnnotRow += " ";

		QRow += q;
		PRow += p;
		}
	}

void RdRpSearcher::Search(const string &QueryLabel, const string &QuerySeq)
	{
	ClearSearch();

	m_QueryLabel = QueryLabel;
	m_QuerySeq = QuerySeq;

	const uint GroupCount = GetGroupCount();
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		SearchGroup(GroupIndex);
	}

void RdRpSearcher::CFilter(const string &QueryLabel, const string &QuerySeq)
	{
	ClearSearch();

	m_QueryLabel = QueryLabel;
	m_QuerySeq = QuerySeq;

	const uint GroupCount = GetGroupCount();
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		CFilterGroup(GroupIndex);
	}

void RdRpSearcher::GetMotifPositions(uint &APos, uint &BPos, uint &CPos) const
	{
	Die("TODO");
	}

void RdRpSearcher::GetTrimmedSeq(string &Seq) const
	{
	Die("TODO");
	//Seq.clear();

	//uint QLo, QHi;
	//GetSpan(QLo, QHi);
	//if (QLo == UINT_MAX)
	//	return;

	//const uint QL = SIZE(m_QuerySeq);
	//asserta(QLo < QHi);
	//asserta(QHi < QL);

	//for (uint Pos = QLo; Pos <= QHi; ++Pos)
	//	Seq += m_QuerySeq[Pos];
	}

uint RdRpSearcher::GetMotifPos(uint MotifIndex) const
	{
	const RPHit *Hit = &m_TopPalmHit.GetHit(MotifIndex);
	if (Hit == 0 || Hit->m_Score == 0)
		return UINT_MAX;
	uint Pos = Hit->m_QPos;
	return Pos;
	}

void RdRpSearcher::GetMotifsSeq(const string &Sep, string &s) const
	{
	s.clear();

	string A, B, C;
	GetMotifSeq(0, A);
	GetMotifSeq(1, B);
	GetMotifSeq(2, C);
	if (A.empty() || B.empty() || C.empty())
		return;

	s = A + Sep + B + Sep + C;
	}

static bool IsNegarna(const string &GroupName)
	{
#define x(Name)	if (GroupName == #Name) return true;
	x(Articulavirales)
	x(Bunyavirales)
	x(Goujianvirales)
	x(Jingchuvirales)
	x(Mononegavirales)
	x(Muvirales)
	x(Serpentovirales)
#undef x
	return false;
	}

void RdRpSearcher::GetSpan(uint &QLo, uint &QHi) const
	{
	m_TopPalmHit.GetSpan(QLo, QHi);
	}

void RdRpSearcher::GetAlnRows(vector<string> &Rows) const
	{
	Rows.clear();
	if (m_TopPalmHit.m_Score <= 0)
		return;

	uint GroupIndex = m_TopPalmHit.m_GroupIndex;
	bool Permuted = m_TopPalmHit.m_Permuted;

	static const uint ABCIndexes[3] = { MOTIF_A, MOTIF_B, MOTIF_C };
	static const uint CABIndexes[3] = { MOTIF_C, MOTIF_A, MOTIF_B };

	Rows.resize(4);

	string &TopLine = Rows[0];
	string &QLine = Rows[1];
	string &ALine = Rows[2];
	string &PLine = Rows[3];

	const uint *MotifIndexes = (Permuted ? CABIndexes : ABCIndexes);
	const RPHit *PrevHit = 0;
	for (uint k = 0; k < 3; ++k)
		{
		uint MotifIndex = MotifIndexes[k];
		char MotifLetter = "ABC"[MotifIndex];

		uint QLo, QHi;
		string QRow, PRow, ARow;
		const RPHit *Hit = &m_TopPalmHit.GetHit(MotifIndex);
		GetAln(*Hit, QRow, PRow, ARow, QLo, QHi);

		uint n = SIZE(QRow);
		string s;
		string t;
		Psa(t, "%c:%u-%u", MotifLetter, QLo, QHi);
		Psa(t, "(%.1f)", Hit->m_Score);
		uint m = SIZE(t);

		if (PrevHit != 0)
			{
			uint IL = Hit->m_QPos - PrevHit->m_QPos + PrevHit->m_PSSM->m_ColCount;
			Psa(s, " <%u> ", IL);
			QLine += s;
			for (uint k = 0; k < SIZE(s); ++k)
				{
				TopLine += ' ';
				PLine += ' ';
				ALine += ' ';
				}
			}
		else
			{
			TopLine += "   ";
			QLine += "   ";
			PLine += "   ";
			ALine += "   ";
			}

		TopLine += t;
		for (uint k = m; k < n; ++k)
			TopLine += ' ';

		Psa(QLine, "%s", QRow.c_str());
		Psa(PLine, "%s", PRow.c_str());
		Psa(ALine, "%s", ARow.c_str());

		if (k < 2)
			{
			for (uint kk = n; kk < m; ++kk)
				{
				QLine += ' ';
				PLine += ' ';
				ALine += ' ';
				}
			}

		PrevHit = Hit;
		}

	uint QLo;
	uint QHi;
	GetSpan(QLo, QHi);
	Psa(QLine, "  [%u]", QHi - QLo + 1);
	}

extern vector<string> g_IncludeNames;
extern vector<string> g_ExcludeNames;

bool ExcludeName(const string &Name)
	{
	for (uint i = 0; i < SIZE(g_ExcludeNames); ++i)
		{
		if (Name == g_ExcludeNames[i])
			return true;
		}
	return false;
	}

bool IncludeName(const string &Name)
	{
	if (g_IncludeNames.empty())
		return true;

	for (uint i = 0; i < SIZE(g_IncludeNames); ++i)
		{
		if (Name == g_IncludeNames[i])
			return true;
		}
	return false;
	}

bool RdRpSearcher::IsHit() const
	{
	if (m_TopPalmHit.m_Score <= opt_minscore)
		return false;

	if (!g_ExcludeNames.empty())
		{
		if (m_TopPalmHit.m_Score > 0)
			{
			uint GroupIndex = m_TopPalmHit.m_GroupIndex;
			string GroupName;
			GetGroupName(GroupIndex, GroupName);
			if (ExcludeName(GroupName))
				return false;
			if (!IncludeName(GroupName))
				return false;
			}
		}
	return true;
	}
