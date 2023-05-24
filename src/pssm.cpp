#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "alpha.h"
#include "sort.h"

void PSSM::SetIsNucleo(const vector<string> &Seqs)
	{
	unsigned CharCount = 0;
	unsigned NucleoCharCount = 0;
	unsigned SeqCount = SIZE(Seqs);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Seq = Seqs[SeqIndex];
		const unsigned L = SIZE(Seq);
		for (unsigned i = 0; i < L; ++i)
			{
			char c = Seq[i];
			++CharCount;
			if (g_IsNucleoChar[c])
				++NucleoCharCount;
			}
		}
	m_IsNucleo = (GetPct(NucleoCharCount, CharCount) > 80.0);
	}

void PSSM::FromSeqs(const vector<string> &Seqs)
	{
	//m_Labels = Labels;
	//m_Seqs = Seqs;
	unsigned SeqCount = SIZE(Seqs);
	asserta(SeqCount > 0);
	m_ColCount = SIZE(Seqs[0]);
	for (unsigned SeqIndex = 1; SeqCount < SeqCount; ++SeqIndex)
		asserta(SIZE(Seqs[SeqIndex]) == m_ColCount);
//	SetIsNucleo(Seqs);
	SetCounts(Seqs);
	SetFreqs();
	SetPseudoFreqs();
	AddPseudo(0.1f);
	SetScores();
	CalcConsSeq(m_ConsSeq);
	}

void PSSM::FromFasta(const string &FileName)
	{
	SeqDB DB;
	DB.FromFasta(FileName);
	FromSeqs(DB.m_Seqs);
	}

void PSSM::SetCounts(const vector<string> &Seqs)
	{
	unsigned ColCount = GetColCount();
	m_Counts.clear();
	m_Counts.resize(ColCount);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		SetCountsCol(Seqs, ColIndex);
	}

void PSSM::SetCountsCol(const vector<string> &Seqs, unsigned ColIndex)
	{
	unsigned AlphaSize = (m_IsNucleo ? 4 : 20);
	const byte *CharToLetter = (m_IsNucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	unsigned SeqCount = SIZE(Seqs);

	vector<unsigned> &Counts = m_Counts[ColIndex];
	Counts.clear();
	Counts.resize(AlphaSize);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		asserta(ColIndex < SIZE(Seqs[SeqIndex]));
		char c = Seqs[SeqIndex][ColIndex];
		unsigned Letter = CharToLetter[c];
		if (Letter >= AlphaSize)
			continue;
		++(Counts[Letter]);
		}
	}

void PSSM::SetFreqs()
	{
	unsigned ColCount = GetColCount();
	m_RawFreqs.clear();
	m_RawFreqs.resize(ColCount);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		SetFreqsCol(ColIndex);
	}

void PSSM::SetFreqsCol(unsigned ColIndex)
	{
	unsigned AlphaSize = (m_IsNucleo ? 4 : 20);

	vector<float> &Freqs = m_RawFreqs[ColIndex];
	Freqs.clear();
	Freqs.resize(AlphaSize);

	vector<unsigned> &Counts = m_Counts[ColIndex];
	unsigned N = 0;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		{
		unsigned n = Counts[Letter];
		N += n;
		}
	if (N == 0)
		{
		// Special case for all Xs
		Warning("PSSM::SetFreqsCol(%u): no valid letters", ColIndex);
		for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
			Counts[Letter] = 1;
		N = AlphaSize;
		}

	float Sum = 0.0;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		{
		unsigned n = Counts[Letter];
		float f = N == 0 ? 0.0f : float(n)/N;
		Sum += f;
		Freqs[Letter] = f;
		}
	asserta(Sum > 0.99f && Sum < 1.01f);
	}

void PSSM::SetPseudoFreqs()
	{
	unsigned AlphaSize = GetAlphaSize();
	m_PseudoFreqs.clear();
	unsigned ColCount = GetColCount();
	m_PseudoFreqs.resize(ColCount);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		m_PseudoFreqs[ColIndex].resize(AlphaSize, 0.0f);
		SetPseudoFreqsCol(ColIndex);
		}
	}

void PSSM::AddPseudo(float w)
	{
	m_PseudoWeight = w;
	unsigned ColCount = GetColCount();
	unsigned AlphaSize = GetAlphaSize();
	m_Freqs.clear();
	m_Freqs.resize(ColCount);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		m_Freqs[ColIndex].resize(AlphaSize);
		AddPseudoCol(ColIndex, w);
		}
	}

void PSSM::AddPseudoCol(unsigned ColIndex, float w)
	{
	unsigned AlphaSize = GetAlphaSize();
	const vector<float> &RawFreqs = m_RawFreqs[ColIndex];
	const vector<float> &PseudoFreqs = m_PseudoFreqs[ColIndex];
	vector<float> &Freqs = m_Freqs[ColIndex];
	asserta(SIZE(RawFreqs) == AlphaSize);
	asserta(SIZE(PseudoFreqs) == AlphaSize);
	asserta(SIZE(Freqs) == AlphaSize);

// Weighted sum
	float Sum = 0.0f;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		{
		float RawFreq = RawFreqs[Letter];
		float PseudoFreq = PseudoFreqs[Letter];
		float Freq = RawFreq + w*PseudoFreq;
		Freqs[Letter] = Freq;
		Sum += Freq;
		}

// Normalize
	float Sum2 = 0.0f;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		{
		float Freq = Freqs[Letter]/Sum;
		Freqs[Letter] = Freq;
		Sum2 += Freq;
		}
	asserta(Sum2 > 0.99f && Sum2 < 1.01f);
	}

void PSSM::SetPseudoFreqsCol(unsigned ColIndex)
	{
	unsigned AlphaSize = GetAlphaSize();
	const vector<float> &Freqs = m_RawFreqs[ColIndex];
	vector<float> &PseudoFreqs = m_PseudoFreqs[ColIndex];
	asserta(SIZE(Freqs) == AlphaSize);
	asserta(SIZE(PseudoFreqs) == AlphaSize);
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		{
		float f = Freqs[Letter];
		float SumJ = 0.0f;
		for (unsigned Letter2 = 0; Letter2 < AlphaSize; ++Letter2)
			{
			float JP = GetJP(Letter, Letter2);
			SumJ += JP;
			float Incf = f*JP;
			PseudoFreqs[Letter2] += Incf;
			}
		asserta(SumJ > 0.99f && SumJ < 1.01f);
		}

	float Sum = 0.0f;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Sum += PseudoFreqs[Letter];
	asserta(Sum > 0.99f && Sum < 1.01f);
	}

const vector<float> &PSSM::GetFreqs(unsigned ColIndex) const
	{
	asserta(ColIndex < SIZE(m_RawFreqs));
	return m_RawFreqs[ColIndex];
	}

void PSSM::LogFreqs(const vector<vector<float> > &FreqsVec) const
	{
	const unsigned AlphaSize = GetAlphaSize();
	const char *LetterToChar = GetLetterToChar();
	const unsigned ColCount = GetColCount();

	Log("Col");
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Log(" %5c", LetterToChar[Letter]);
	Log("\n");

	Log("---");
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Log(" -----");
	Log("\n");

	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		Log("%3u", ColIndex);
		const vector<float> &Freqs = FreqsVec[ColIndex];
		for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
			{
			float f = Freqs[Letter];
			if (f == 0)
				Log("     .");
			else
				Log(" %5.3f", f);
			}
		string Logo;
		Log("  %s\n", GetLogo(Freqs, Logo).c_str());
		}
	}

void PSSM::LogMe() const
	{
	unsigned ColCount = GetColCount();
	unsigned AlphaSize = GetAlphaSize();
	const char *LetterToChar = GetLetterToChar();
	//string ConsSeq;
	//CalcConsSeq(ConsSeq);

	Log("%u cols  %s\n", ColCount, m_ConsSeq.c_str());

	if (!m_RawFreqs.empty())
		{
		Log("\n");
		Log("Raw freqs:\n");
		LogFreqs(m_RawFreqs);

		Log("\n");
		Log("Pseudo freqs, weight %.3f:\n", m_PseudoWeight);
		LogFreqs(m_PseudoFreqs);

		Log("\n");
		Log("Freqs:\n");
		LogFreqs(m_Freqs);
		}

	if (m_Scores.empty())
		return;

	Log("\n");
	Log("Scores:\n");
	Log("Col");
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Log(" %5c", LetterToChar[Letter]);
	Log(" %5.5s", "Wild");
	Log("\n");

	Log("---");
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		Log(" -----");
	Log(" -----");
	Log("\n");
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		Log("%3u", ColIndex);
		const vector<float> &Scores = m_Scores[ColIndex];
		for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
			{
			float s = m_Scores[ColIndex][Letter];
			Log(" %5.2f", s);
			}
		Log(" %5.2f", Scores[AlphaSize]);
		if (m_Freqs.empty())
			Log(" ");
		else
			{
			string Logo;
			Log("  %s ", GetLogo(m_Freqs[ColIndex], Logo).c_str());
			}

		vector<unsigned> Order;
		SortDescending(Scores, Order);

		for (unsigned i = 0; i < AlphaSize; ++i)
			{
			unsigned Letter = Order[i];
			float Score = Scores[Letter];
			if (Score <= 0.0f)
				break;
			Log(" %c(%.2f)", LetterToChar[Letter], Score);
			}
		Log("\n");
		}
	}

const char *PSSM::CalcConsSeq(string &Seq) const
	{
	Seq.clear();
	uint ColCount = GetColCount();
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		char c = CalcConsChar(ColIndex);
		Seq += c;
		}
	return Seq.c_str();
	}

char PSSM::CalcConsChar(uint ColIndex) const
	{
	const vector<float> &Freqs = GetFreqs(ColIndex);
	const uint N = SIZE(Freqs);
	asserta(N == 20);
	float TopFreq = 0;
	uint TopLetter = 0;
	for (uint i = 0; i < N; ++i)
		{
		if (Freqs[i] > TopFreq)
			{
			TopFreq = Freqs[i];
			TopLetter = i;
			}
		}
	const char *LetterToChar = GetLetterToChar();
	char c = LetterToChar[TopLetter];
	if (TopFreq < 0.5)
		c = tolower(c);
	return c;
	}

char PSSM::GetConsChar(uint ColIndex) const
	{
	asserta(ColIndex < SIZE(m_ConsSeq));
	char c = m_ConsSeq[ColIndex];
	return c;
	}

const string &PSSM::GetLogo(const vector<float> &Freqs, string &Logo, unsigned n) const
	{
	Logo.clear();
	vector<unsigned> Order;
	SortAscending(Freqs, Order);
	unsigned AlphaSize = GetAlphaSize();
	asserta(SIZE(Order) == AlphaSize);
	const char *LetterToChar = GetLetterToChar();
	char c = '.';
	for (unsigned i = 0; i < AlphaSize; ++i)
		{
		unsigned Letter = Order[i];
		float f = Freqs[Letter];
		c = LetterToChar[Letter];
		unsigned w = unsigned(f*n);
		for (unsigned k = 0; k < w; ++k)
			Logo.push_back(c);
		}
	while (SIZE(Logo) < n)
		Logo.push_back(c);
	return Logo;
	}

void PSSM::SetScores()
	{
	unsigned AlphaSize = GetAlphaSize();
	m_Scores.clear();
	unsigned ColCount = GetColCount();
	m_Scores.resize(ColCount);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		m_Scores[ColIndex].resize(AlphaSize+1);
		SetScoresCol(ColIndex);
		}
//	SetBackgroundScore();
	}

void PSSM::CalcProbs(uint ColIndex, vector<double> &Probs) const
	{
	Probs.clear();
	extern double g_AAFreqs[20];
	asserta(ColIndex < SIZE(m_Scores));
	asserta(SIZE(m_Scores[ColIndex]) == 21);
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		// double Score = log(f) - log((float) g_AAFreqs[Letter]);
		double Score = m_Scores[ColIndex][Letter];
		double LogBack = log(g_AAFreqs[Letter]);
		double LogFreq = Score + LogBack;
		double Prob = exp(LogFreq);
		Probs.push_back(Prob);
		}
	}

void PSSM::CalcProbs(uint ColIndex, vector<float> &Probs) const
	{
	Probs.clear();
	vector<double> dProbs;
	CalcProbs(ColIndex, dProbs);
	asserta(SIZE(dProbs) == 20);
	for (uint i = 0; i < 20; ++i)
		{
		float P = (float) dProbs[i];
		Probs.push_back(P);
		}
	}

void PSSM::SetScoresCol(unsigned ColIndex)
	{
	extern double g_AAFreqs[20];
	unsigned AlphaSize = GetAlphaSize();
	float WildScore = 0.0f;
	for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
		{
		asserta(ColIndex < SIZE(m_Freqs));
		asserta(ColIndex < SIZE(m_Scores));
		asserta(Letter < SIZE(m_Freqs[ColIndex]));
		asserta(Letter < SIZE(m_Scores[ColIndex]));
		float f = m_Freqs[ColIndex][Letter];
		if (m_IsNucleo)
			{
			double Score = log(f) - log((float) 0.25);
			WildScore += float(Score*0.25);
			m_Scores[ColIndex][Letter] = float(Score);
			}
		else
			{
			double Score = log(f) - log((float) g_AAFreqs[Letter]);
			WildScore += float(Score*g_AAFreqs[Letter]);
			m_Scores[ColIndex][Letter] = float(Score);
			}
		}
	asserta(AlphaSize < SIZE(m_Scores[ColIndex]));
	if (WildScore > 0.0)
		WildScore = 0.0;
	m_Scores[ColIndex][AlphaSize] = WildScore;
	}

float PSSM::GetScore1(uint ColIndex, char c) const
	{
	const byte *CharToLetter = GetCharToLetter();
	if (c == '*')
		return (float) opt_stop_score;
	unsigned AlphaSize = GetAlphaSize();
	unsigned Letter = CharToLetter[c];
	if (Letter == INVALID_LETTER)
		return m_Scores[ColIndex][AlphaSize];
	asserta(Letter < AlphaSize);
	return m_Scores[ColIndex][Letter];
	}

float PSSM::GetScore(const char *Seq) const
	{
	unsigned AlphaSize = GetAlphaSize();
	const byte *CharToLetter = GetCharToLetter();
	unsigned ColCount = GetColCount();
	float Score = 0.0f;
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		char c = Seq[ColIndex];
		if (c == '*')
			Score += (float) opt_stop_score;
		else
			{
			unsigned Letter = CharToLetter[c];
			if (Letter == INVALID_LETTER)
				Score += m_Scores[ColIndex][AlphaSize];
			else
				{
				asserta(Letter < AlphaSize);
				Score += m_Scores[ColIndex][Letter];
				}
			}
		}
	return Score;
	}

const char *PSSM::GetRandomSeq(string &s) const
	{
	s.clear();
	unsigned AlphaSize = GetAlphaSize();
	const char *LetterToChar = GetLetterToChar();
	unsigned ColCount = GetColCount();
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		byte Letter = GetRandomLetter();
		char c = LetterToChar[Letter];
		s.push_back(c);
		}
	return s.c_str();
	}

byte PSSM::GetRandomLetter() const
	{
	unsigned AlphaSize = GetAlphaSize();
	return rand()%AlphaSize;
	}

void PSSM::ToFile(const string &FileName) const
	{
	FILE *f = CreateStdioFile(FileName);
	ToFile(f);
	CloseStdioFile(f);
	}

void PSSM::ToFile(FILE *f) const
	{
	unsigned ColCount = GetColCount();
	unsigned AlphaSize = GetAlphaSize();
	const char *LetterToChar = GetLetterToChar();

	fprintf(f, "PSSM %s %u\n", m_IsNucleo ? "nt" : "aa", ColCount);
	fprintf(f, "consseq	%s\n", m_ConsSeq.c_str());
	if (m_IsNucleo)
		fprintf(f, "Col     A     C     G     T     N\n");
	else
		fprintf(f, "Col     A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y     X\n");

	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		fprintf(f, "%3u", ColIndex + 1);
		const vector<float> &Scores = m_Scores[ColIndex];
		for (unsigned Letter = 0; Letter < AlphaSize; ++Letter)
			{
			float s = m_Scores[ColIndex][Letter];
			fprintf(f, " %5.2f", s);
			}
		fprintf(f, " %5.2f", Scores[AlphaSize]);
		if (!m_Freqs.empty())
			{
			string Logo;
			fprintf(f, "  %s ", GetLogo(m_Freqs[ColIndex], Logo).c_str());
			}

		vector<unsigned> Order;
		SortDescending(Scores, Order);

		for (unsigned i = 0; i < AlphaSize; ++i)
			{
			unsigned Letter = Order[i];
			float Score = Scores[Letter];
			if (Score <= 0.0f)
				break;
			fprintf(f, " %c(%.2f)", LetterToChar[Letter], Score);
			}
		fprintf(f, "\n");
		}
	fprintf(f, "//\n");
	}

void PSSM::FromFile(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);
	FromFile(f);
	CloseStdioFile(f);
	}

void PSSM::FromFile(FILE *f)
	{
	vector<string> Strings;
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		Strings.push_back(Line);
		if (Line == "//")
			break;
		}
	FromStrings("FILE*", Strings);
	}

static bool GetNextString(const vector<string> &Strings, const unsigned N,
  unsigned &i, string &s)
	{
	if (i >= N)
		return false;
	s = Strings[i++];
	return true;
	}

void PSSM::FromCStrings(const string &Name, const char **CStrings)
	{
	vector<string> Strings;
	for (unsigned i = 0; i < 32; ++i)
		{
		const string s = string(CStrings[i]);
		Strings.push_back(s);
		asserta(!s.empty());
		if (s[0] == '/')
			{
			FromStrings(Name, Strings);
			return;
			}
		}
	Die("Missing // in psm %s", Name.c_str());
	}

void PSSM::FromStrings(const string &Name, const vector<string> &Strings)
	{
	Clear();

	const unsigned N = SIZE(Strings);
	string Line;
	unsigned LineNr = 0;
	bool Ok = GetNextString(Strings, N, LineNr, Line);
	if (!Ok)
		Die("PSSM::FromString, empty .psm %s", Name.c_str());

	vector<string> Fields;
	Split(Line, Fields, ' ');
	if (SIZE(Fields) != 3 || Fields[0] != "PSSM")
		Die("Invalid header in .psm %s", Name.c_str());

	unsigned AlphaSize = 0;
	if (Fields[1] == "nt")
		{
		m_IsNucleo = true;
		AlphaSize = 4;
		}
	else if (Fields[1] == "aa")
		{
		m_IsNucleo = false;
		AlphaSize = 20;
		}
	else
		Die("Invalid alphabet in .psm %s", Name.c_str());

	m_ColCount = atou(Fields[2]);
	if (m_ColCount == 0)
		Die("Zero cols in .psm file %s", Name.c_str());

	Ok = GetNextString(Strings, N, LineNr, Line);
	if (!Ok)
		Die("Premature end-of-file in .psm %s", Name.c_str());

	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 2)
		Die("Bad consseq line %s", Line.c_str());
	m_ConsSeq = Fields[1];

	Ok = GetNextString(Strings, N, LineNr, Line);
	if (!Ok)
		Die("Premature end-of-file in .psm %s", Name.c_str());

	if (!StartsWith(Line, "Col "))
		Die("Invalid 2nd header line in .psm %s", Name.c_str());

	m_Scores.resize(m_ColCount);
	for (unsigned ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		{
		m_Scores[ColIndex].resize(AlphaSize+1, 0.0f);
		Ok = GetNextString(Strings, N, LineNr, Line);
		if (!Ok)
			Die("Premature end-of-file in .psm %s", Name.c_str());
		Split(Line, Fields, 0);
		if (SIZE(Fields) < AlphaSize + 2 || atou(Fields[0]) != ColIndex + 1)
			Die("Invalid scores for col %u in .psm  %s", ColIndex, Name.c_str());

		for (unsigned Letter = 0; Letter <= AlphaSize; ++Letter)
			{
			const string &f = Fields[Letter+1];
			float Score = (float) StrToFloat(f.c_str());
			m_Scores[ColIndex][Letter] = Score;
			}
		}
	Ok = GetNextString(Strings, N, LineNr, Line);
	if (!Ok || Line != "//")
		Die("Missing // at end of .psm %s", Name.c_str());
	}

 // joint prob, sums to 1 for each Row, not whole Mx!!
float PSSM::GetJP(unsigned Letter1, unsigned Letter2) const
	{
	extern double g_BLOSUM62_JP_Rows[20][20];
	extern double g_NtMx_JP_Rows[4][4];
	if (m_IsNucleo)
		{
		asserta(Letter1 < 4 && Letter2 < 4);
		return (float) g_NtMx_JP_Rows[Letter1][Letter2];
		}
	else
		{
		asserta(Letter1 < 20 && Letter2 < 20);
		return (float) g_BLOSUM62_JP_Rows[Letter1][Letter2];
		}
	}
