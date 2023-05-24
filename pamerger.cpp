#include "myutils.h"
#include "pamerger.h"
#include "cmp.h"

void PAMerger::LogFields()
	{
#define F(x)	Log("STR  " #x "\n");
#include "pastrs.h"

#define F(x)	Log("INT  " #x "\n");
#include "paints.h"

#define F(x)	Log("DBL  " #x "\n");
#include "padbls.h"
	}

void PAMerger::ToFasta(FILE *f) const
	{
	if (f == 0)
		return;

	if (aaseq == "")
		return;
	if (m_FullTrimLo == UINT_MAX)
		return;

	const uint L = SIZE(aaseq);
	asserta(m_FullTrimLo > 0);
	asserta(m_FullTrimHi <= L);

	string OutSeq;
	for (uint Pos = m_FullTrimLo; Pos <= m_FullTrimHi; ++Pos)
		OutSeq += aaseq[Pos-1];
	asserta(SIZE(OutSeq) > 0);

	string OutLabel = m_Label;

	if (m_MotifA != "")
		Psa(OutLabel, " A:%u:%s", m_MotifPosA - m_FullTrimLo + 1, m_MotifA.c_str());
	if (m_MotifB != "")
		Psa(OutLabel, " B:%u:%s", m_MotifPosB - m_FullTrimLo + 1, m_MotifB.c_str());
	if (m_MotifC != "")
		Psa(OutLabel, " C:%u:%s", m_MotifPosC - m_FullTrimLo + 1, m_MotifC.c_str());

	SeqToFasta(f, OutLabel, OutSeq.c_str());
	}

void PAMerger::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	if (m_Label == "")
		return;

	double RdRpScore = GetRdRpScore();

	static bool HdrDone = false;
	if (!HdrDone)
		{
		fprintf(f, "score");
		fprintf(f, "\tlabel");
		fprintf(f, "\tlength");
		fprintf(f, "\tmotifA");
		fprintf(f, "\tmotifB");
		fprintf(f, "\tmotifC");
		fprintf(f, "\tposA");
		fprintf(f, "\tposB");
		fprintf(f, "\tposC");
		fprintf(f, "\textlo");
		fprintf(f, "\texthi");
		fprintf(f, "\tpplo");
		fprintf(f, "\tpphi");
		fprintf(f, "\n");

		HdrDone = true;
		}

	fprintf(f, "%.1f", RdRpScore);
	fprintf(f, "\t%s", m_Label.c_str());
	fprintf(f, "\t%s", m_MotifA == "" ? "." : m_MotifA.c_str());
	fprintf(f, "\t%s", m_MotifB == "" ? "." : m_MotifB.c_str());
	fprintf(f, "\t%s", m_MotifC == "" ? "." : m_MotifC.c_str());
	fprintf(f, "\t%u", m_MotifPosA == UINT_MAX ? 0 : m_MotifPosA);
	fprintf(f, "\t%u", m_MotifPosB == UINT_MAX ? 0 : m_MotifPosB);
	fprintf(f, "\t%u", m_MotifPosC == UINT_MAX ? 0 : m_MotifPosC);
	fprintf(f, "\t%u", m_FullTrimLo == UINT_MAX ? 0 : m_FullTrimLo);
	fprintf(f, "\t%u", m_FullTrimHi == UINT_MAX ? 0 : m_FullTrimHi);
	fprintf(f, "\t%u", m_PPLo == UINT_MAX ? 0 : m_PPLo);
	fprintf(f, "\t%u", m_PPHi == UINT_MAX ? 0 : m_PPHi);
	fprintf(f, "\n");
	}

void PAMerger::ToFev(FILE *f) const
	{
	if (f == 0)
		return;
	if (m_Label == "")
		return;

	if (opt_compact)
		fputs(m_Label.c_str(), f);
	else
		{
		asserta(m_Line != 0);
		fputs(m_Line->c_str(), f);
		}

	if (m_MotifA != "")
		fprintf(f, "\tseqA=%s", m_MotifA.c_str());

	if (m_MotifB != "")
		fprintf(f, "\tseqB=%s", m_MotifB.c_str());

	if (m_MotifC != "")
		fprintf(f, "\tseqC=%s", m_MotifC.c_str());

	if (m_MotifPosA != UINT_MAX)
		fprintf(f, "\tposA=%u", m_MotifPosA);

	if (m_MotifPosB != UINT_MAX)
		fprintf(f, "\tposB=%u", m_MotifPosB);

	if (m_MotifPosC != UINT_MAX)
		fprintf(f, "\tposC=%u", m_MotifPosC);

	string DGD;
	GetDGD(DGD);
	if (DGD == "DGD")
		fprintf(f, "\tdgd=yes");

	string GDD;
	GetGDD(GDD);
	if (GDD != "XXX")
		fprintf(f, "\tgdd=%s", GDD.c_str());

	char Gate = GetGate();
	if (Gate != 'X')
		fprintf(f, "\tgate=%c", Gate);

	double P_gate = CMP::GetRdRpProb_Gate(Gate);
	if (P_gate >= 0.8 || P_gate <= 0.2)
		fprintf(f, "\tgate_prob=%.2f", P_gate);

	double P_gdd = CMP::GetRdRpProb_GDD(GDD);
	if (P_gdd >= 0.8 || P_gdd <= 0.2)
		fprintf(f, "\tgdd_prob=%.2f", P_gdd);

	double P = (P_gate + P_gdd)/2;
	if (P > 0.8 || P < 0.2)
		fprintf(f, "\tgate_gdd_prob=%.2f", P);

	if (m_FullTrimLo != UINT_MAX)
		{
		asserta(m_FullTrimHi != UINT_MAX);
		fprintf(f, "\text_lo=%u", m_FullTrimLo);
		fprintf(f, "\text_hi=%u", m_FullTrimHi);
		fprintf(f, "\text_method=%s", m_FullTrimMethod.c_str());
		}

	if (m_PPLo != UINT_MAX)
		{
		asserta(m_PPHi != UINT_MAX);
		fprintf(f, "\tpp_lo=%u", m_PPLo);
		fprintf(f, "\tpp_hi=%u", m_PPHi);
		}

	double RdRpScore = GetRdRpScore();
	if (RdRpScore != DBL_MAX)
		fprintf(f, "\trdrp=%.1f", RdRpScore);

	fputc('\n', f);
	}

void PAMerger::LogMe() const
	{
	Log("\n");
	Log(">%s\n", m_Label.c_str());
#define F(x)	if (x != "") Log("%+16.16s  " #x "\n", x.c_str());
#include "pastrs.h"

#define F(x)	if (x != UINT_MAX) Log("%16u  " #x "\n", x);
#include "paints.h"

#define F(x)	if (x != DBL_MAX) Log("%16.3g  " #x "\n", x);
#include "padbls.h"

	if (m_MotifA != "")
		Log("%16.16s  A\n", m_MotifA.c_str());

	if (m_MotifB != "")
		Log("%16.16s  B\n", m_MotifB.c_str());

	if (m_MotifC != "")
		Log("%16.16s  C\n", m_MotifC.c_str());

	if (m_MotifPosA != UINT_MAX)
		Log("%16u  PosA\n", m_MotifPosA);

	if (m_MotifPosB != UINT_MAX)
		Log("%16u  PosB\n", m_MotifPosB);

	if (m_MotifPosC != UINT_MAX)
		Log("%16u  PosC\n", m_MotifPosC);

	if (m_FullTrimLo != UINT_MAX)
		Log("%16u  FullTrimLo\n", m_FullTrimLo);

	if (m_FullTrimHi != UINT_MAX)
		Log("%16u  FullTrimHi\n", m_FullTrimHi);

	if (m_FullTrimMethod != "")
		Log("%16s  FullTrimMethod\n", m_FullTrimMethod);

	if (m_PPLo != UINT_MAX)
		Log("%16u  PPLo\n", m_PPLo);

	if (m_PPHi != UINT_MAX)
		Log("%16u  PPHi\n", m_PPHi);
	}

void PAMerger::FromLine(const string &Line)
	{
	Reset();

	m_Line = &Line;
	vector<string> Fields;

	Split(Line, Fields, '\t');
	const uint N = SIZE(Fields);
	m_Label = Fields[0];
	for (uint i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		size_t n = Field.find('=');
		if (n == string::npos)
			continue;
		const string &Name = Field.substr(0, n);
		const string &Value = Field.substr(n+1);
		m_NameToValue[Name] = Value;
		}

#define F(x)	SetStr(#x, x);
#include "pastrs.h"

#define F(x)	SetInt(#x, x);
#include "paints.h"

#define F(x)	SetDbl(#x, x);
#include "padbls.h"

	SetMotifs();
	SetFullTrim();
	SetPP();
	}

void PAMerger::SetInt(const string &Name, uint &Value)
	{
	map<string, string>::const_iterator p =
	  m_NameToValue.find(Name);
	if (p == m_NameToValue.end())
		Value = UINT_MAX;
	else
		Value = StrToUint(p->second);
	}

void PAMerger::SetStr(const string &Name, string &Value)
	{
	map<string, string>::const_iterator p =
	  m_NameToValue.find(Name);
	if (p == m_NameToValue.end())
		Value = "";
	else
		Value = p->second;
	}

void PAMerger::SetDbl(const string &Name, double &Value)
	{
	map<string, string>::const_iterator p =
	  m_NameToValue.find(Name);
	if (p == m_NameToValue.end())
		Value = DBL_MAX;
	else
		Value = StrToFloat(p->second);
	}

bool PAMerger::PSSMGroupIsRT(const string &Group)
	{
	if (Group.size() < 2)
		return false;
	return Group[0] == 'R' && Group[1] == 'T';
	}

bool PAMerger::MotifHMMIsRT(const string &HMM)
	{
	if (HMM.size() < 2)
		return false;
	return HMM.substr(0, 2) == "PN";
	}

bool PAMerger::MotifHMMPermuted() const
	{
	if (motif_hmm == "" || motif_hmm_evalue > 1e-3)
		return false;

	if (motif_hmm_posC == UINT_MAX || motif_hmm_posA == UINT_MAX)
		return false;

	return motif_hmm_posC < motif_hmm_posA;
	}

char PAMerger::GetGate() const
	{
	if (m_MotifA.size() < 8)
		return 'X';
	return m_MotifA[3+5];
	}

const char *PAMerger::GetGDD(string &GDD) const
	{
	GDD = "XXX";
	if (m_MotifC.size() >= 6)
		{
		GDD[0] = m_MotifC[2];
		GDD[1] = m_MotifC[3];
		GDD[2] = m_MotifC[4];
		}
	return GDD.c_str();
	}

const char *PAMerger::GetDGD(string &DGD) const
	{
	DGD = "XXX";
	if (m_MotifA.size() >= 8)
		DGD[0] = m_MotifA[3];
	if (m_MotifB.size() >= 10)
		DGD[1] = m_MotifB[1];
	if (m_MotifC.size() >= 6)
		DGD[2] = m_MotifC[3];
	return DGD.c_str();
	}

void PAMerger::SetMotifs()
	{
	SetMotif(motif_hmm_seqA, pssm_seqA, dmnd_seqA, m_MotifA);
	SetMotif(motif_hmm_seqB, pssm_seqB, dmnd_seqB, m_MotifB);
	SetMotif(motif_hmm_seqC, pssm_seqC, dmnd_seqC, m_MotifC);

	SetMotifPos(m_MotifA, m_MotifPosA);
	SetMotifPos(m_MotifB, m_MotifPosB);
	SetMotifPos(m_MotifC, m_MotifPosC);
	}

void PAMerger::SetMotif(const string &Seq1, const string &Seq2,
  const string &Seq3, string &Motif)
	{
	vector<string> Seqs;
	Seqs.push_back(Seq1);
	Seqs.push_back(Seq2);
	Seqs.push_back(Seq3);
	SetMotif(Seqs, Motif);
	}

void PAMerger::SetMotifPos(const string &MotifSeq, uint &Pos) const
	{
	Pos = UINT_MAX;
	if (SIZE(MotifSeq) < 4)
		return;
	uint n = (uint) aaseq.find(MotifSeq);
	uint r = (uint) aaseq.rfind(MotifSeq);
	if (n == r)
		Pos = n + 1;
	}

void PAMerger::SetMotif(const vector<string> &Seqs,
  string &Motif) const
	{
	const uint N = SIZE(Seqs);
	map<string, uint> SeqToCount;
	uint MaxCount = 0;
	Motif = "";
	for (uint i = 0; i < N; ++i)
		{
		const string &Seq = Seqs[i];
		if (Seq == "")
			continue;
		if (SeqToCount.find(Seq) == SeqToCount.end())
			SeqToCount[Seq] = 0;
		uint n = SeqToCount[Seq] + 1;
		SeqToCount[Seq] = n;
		if (n > MaxCount)
			{
			MaxCount = n;
			Motif = Seq;
			}
		}
	}

static uint GetMin2(uint x, uint y)
	{
	if (x != UINT_MAX && y != UINT_MAX)
		return min(x, y);
	else if (x != UINT_MAX && y == UINT_MAX)
		return x;
	else if (x == UINT_MAX && y != UINT_MAX)
		return y;
	asserta(x == UINT_MAX && y == UINT_MAX);
	return UINT_MAX;
	}

static uint GetMin3(uint x, uint y, uint z)
	{
	return GetMin2(x, GetMin2(y, z));
	}

static uint GetMax2(uint x, uint y)
	{
	if (x != UINT_MAX && y != UINT_MAX)
		return max(x, y);
	else if (x != UINT_MAX && y == UINT_MAX)
		return x;
	else if (x == UINT_MAX && y != UINT_MAX)
		return y;
	asserta(x == UINT_MAX && y == UINT_MAX);
	return UINT_MAX;
	}

static uint GetMax3(uint x, uint y, uint z)
	{
	return GetMax2(x, GetMax2(y, z));
	}

uint PAMerger::GetLoMotifPos() const
	{
	return GetMin3(m_MotifPosA, m_MotifPosB, m_MotifPosC);
	}

uint PAMerger::GetHiMotifPos() const
	{
	return GetMax3(m_MotifPosA, m_MotifPosB, m_MotifPosC);
	}

uint PAMerger::GetMinHitLo() const
	{
	return GetMin3(hmm_rdrp_plus_lo, hmm_rdrp_minus_lo, dmnd_lo);
	}

uint PAMerger::GetMaxHitHi() const
	{
	return GetMin3(hmm_rdrp_plus_hi, hmm_rdrp_minus_hi, dmnd_hi);
	}

void PAMerger::SetFullTrim()
	{
	m_FullTrimLo = UINT_MAX;
	m_FullTrimHi = UINT_MAX;

	uint L = SIZE(aaseq);
	if (L == 0)
		return;
	if (m_MotifPosA != UINT_MAX && m_MotifPosC != UINT_MAX
	  && m_MotifPosA < m_MotifPosC)
		{
		if (m_MotifPosA <= 150)
			m_FullTrimLo = 1;
		else
			m_FullTrimLo = m_MotifPosA - 150;

		m_FullTrimHi = m_MotifPosC + 150;
		if (m_FullTrimHi > L)
			m_FullTrimHi = L;
		m_FullTrimMethod = "150AC150";
		return;
		}

	if (m_MotifPosC != UINT_MAX && m_MotifPosB != UINT_MAX
	  && m_MotifPosC < m_MotifPosB)
		{
		if (m_MotifPosC <= 150)
			m_FullTrimLo = 1;
		else
			m_FullTrimLo = m_MotifPosC - 150;

		m_FullTrimHi = m_MotifPosB + 150;
		if (m_FullTrimHi > L)
			m_FullTrimHi = L;
		m_FullTrimMethod = "150CB150";
		return;
		}

	if (m_MotifPosB != UINT_MAX)
		{
		if (m_MotifPosB <= 200)
			m_FullTrimLo = 1;
		else
			m_FullTrimLo = m_MotifPosB - 200;

		m_FullTrimHi = m_MotifPosB + 200;
		if (m_FullTrimHi > L)
			m_FullTrimHi = L;
		m_FullTrimMethod = "200B200";
		return;
		}

	if (m_MotifPosC != UINT_MAX)
		{
		if (m_MotifPosC <= 250)
			m_FullTrimLo = 1;
		else
			m_FullTrimLo = m_MotifPosC - 250;

		m_FullTrimHi = m_MotifPosC + 150;
		if (m_FullTrimHi > L)
			m_FullTrimHi = L;
		m_FullTrimMethod = "250C150";
		return;
		}

	if (m_MotifPosA != UINT_MAX)
		{
		if (m_MotifPosA <= 150)
			m_FullTrimLo = 1;
		else
			m_FullTrimLo = m_MotifPosA - 150;

		m_FullTrimHi = m_MotifPosA + 250;
		if (m_FullTrimHi > L)
			m_FullTrimHi = L;
		m_FullTrimMethod = "150A250";
		return;
		}

	uint Lo = GetMinHitLo();
	uint Hi = GetMaxHitHi();
	if (Lo != UINT_MAX && Hi != UINT_MAX)
		{
		m_FullTrimLo = Lo;
		m_FullTrimHi = Hi;
		m_FullTrimMethod = "aln";
		}
	}

void PAMerger::SetPP()
	{
	m_PPLo = UINT_MAX;
	m_PPHi = UINT_MAX;

	if (m_MotifPosA != UINT_MAX && m_MotifPosB != UINT_MAX
	  && m_MotifPosA < m_MotifPosC)
		{
		m_PPLo = m_MotifPosA;
		m_PPHi = m_MotifPosC + 7;
		}
	else if (m_MotifPosC != UINT_MAX && m_MotifPosB != UINT_MAX
	  && m_MotifPosC < m_MotifPosB)
		{
		m_PPLo = m_MotifPosC;
		m_PPHi = m_MotifPosB + 13;
		}
	}

double PAMerger::GetLE(double E)
	{
	if (E == DBL_MAX)
		return 0;
	if (E > 1e-3)
		return 0;
	if (E < 1e-100)
		return 100;
	double LE = -log10(E);
	return LE;
	}

bool PAMerger::HMMMinusWins() const
	{
	if (motif_hmm != "" && MotifHMMIsRT(motif_hmm))
		return true;
	double Ep = hmm_rdrp_plus_evalue;
	double Em = hmm_rdrp_minus_evalue;
	if (Em != DBL_MAX)
		{
		if (Ep == DBL_MAX)
			return true;
		if (Em < Ep)
			return true;
		}
	return false;
	}

double PAMerger::GetRdRpScore() const
	{
	if (HMMMinusWins())
		return 0;
	double LEp = GetBestPlusLE();
	double LEm = GetBestMinusLE();

	if (LEp != 0 && LEm == 0)
		{
		if (LEp >= 10)
			return 100;
		return 50 + LEp*5.0;
		}
	else if (LEp == 0 && LEm != 0)
		return 0;
	else if (LEp != 0 && LEm != 0)
		{
		if (LEm - LEp > 1)
			return 0;
		else if (LEm - LEp > 0)
			return 50;
		else if (LEp - LEm > 10)
			return 100;
		double d = LEp - LEm;
		asserta(d >= 0 && d <= 10);
		return 50 + 5*d;
		}
	else if (LEp == 0 && LEm == 0)
		return 0;
	asserta(false);
	return 0;
	}

double PAMerger::GetPSSME() const
	{
	if (pssm_score == DBL_MAX || pssm_score < 5)
		return 1;
	double LE = pssm_score/5.0;
	double E = pow(10, -LE);
	return E;
	}

double PAMerger::GetBestPlusE() const
	{
	if (motif_hmm != "" && MotifHMMIsRT(motif_hmm))
		return 1;

	double E = hmm_rdrp_plus_evalue;
	if (motif_hmm != "" && !MotifHMMIsRT(motif_hmm))
		E = min(E, motif_hmm_evalue);
	if (E != DBL_MAX)
		return E;

	if (PSSMGroupIsRT(pssm_group))
		return 1;

	double PSSME = GetPSSME();
	double E_dmnd = dmnd_evalue;
	double MinE = min(PSSME, E_dmnd);

	char Gate = GetGate();
	double P_gate = CMP::GetRdRpProb_Gate(Gate);
	if (P_gate <= 0.2)
		MinE *= 1e5;

	string GDD;
	GetGDD(GDD);
	double P_gdd = CMP::GetRdRpProb_GDD(GDD);
	if (P_gdd <= 0.2)
		MinE *= 1e5;

	return MinE;
	}

double PAMerger::GetBestMinusE() const
	{
	double E = hmm_rdrp_minus_evalue;
	if (motif_hmm != "" && MotifHMMIsRT(motif_hmm))
		E = min(E, motif_hmm_evalue);
	if (E != DBL_MAX)
		return E;

	if (!PSSMGroupIsRT(pssm_group))
		return 1;

	double PSSME = GetPSSME();
	return PSSME;
	}

double PAMerger::GetBestPlusLE() const
	{
	double E = GetBestPlusE();
	if (E > 0.1)
		return 0;
	if (E < 1e-100)
		return 100;
	double LE = -log10(E);
	return LE;
	}

double PAMerger::GetBestMinusLE() const
	{
	double E = GetBestMinusE();
	if (E > 0.1)
		return 0;
	if (E < 1e-100)
		return 100;
	double LE = -log10(E);
	return LE;
	}
