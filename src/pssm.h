#ifndef pssm_h
#define pssm_h

#include "alpha.h"

class SeqDB;

class PSSM
	{
public:
	bool m_IsNucleo;
	unsigned m_ColCount;
	vector<vector<float> > m_Scores;

private:
	float m_PseudoWeight;
	vector<vector<unsigned> > m_Counts;
	vector<vector<float> > m_RawFreqs;
	vector<vector<float> > m_PseudoFreqs;
	vector<vector<float> > m_Freqs;
	string m_ConsSeq;

public:
	PSSM()
		{
		Clear();
		}

	void Clear()
		{
		m_IsNucleo = false;
		m_ColCount = 0;
		m_PseudoWeight = 0.0f;

		m_Scores.clear();
		m_Counts.clear();
		m_RawFreqs.clear();
		m_PseudoFreqs.clear();
		}

	void FromFasta(const string &FileName);
	void FromFile(const string &FileName);
	void ToFile(FILE *f) const;
	void ToFile(const string &FileName) const;
	void FromFile(FILE *f);
	void FromStrings(const string &Name, const vector<string> &Strings);
	void FromCStrings(const string &Name, const char *CStrings[]);
	float GetJP(unsigned Letter1, unsigned Letter2) const; // joint prob, sums to 1 for each Row, not whole Mx!!

	unsigned GetColCount() const
		{
		return m_ColCount;
		}

	unsigned GetAlphaSize() const
		{
		return (m_IsNucleo ? 4 : 20);
		}

	const char *GetLetterToChar() const
		{
		if (m_IsNucleo)
			return (const char *) g_LetterToCharNucleo;
		else
			return (const char *) g_LetterToCharAmino;
		}

	const byte *GetCharToLetter() const
		{
		return (m_IsNucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
		}

	void FromSeqs(const vector<string> &Labels, const vector<string> &Seqs);
	void LogMe() const;
	const vector<float> &GetFreqs(unsigned ColIndex) const;
	const char *GetRandomSeq(string &s) const;
	byte GetRandomLetter() const;
	float GetScore(const char *Seq) const;
	float GetScore1(uint Col, char c) const;

	const string &GetLogo(const vector<float> &Freqs, string &Logo,
	  unsigned n = 16) const;
	char GetConsChar(uint ColIndex) const;
	char CalcConsChar(uint ColIndex) const;
	const char *CalcConsSeq(string &Seq) const;

private:
	void LogFreqs(const vector<vector<float> > &Freqs) const;
	void SetIsNucleo(const vector<string> &Seqs);

	void SetCounts(const vector<string> &Seqs);
	void SetCountsCol(const vector<string> &Seqs, unsigned ColIndex);

	void SetScores();
	void SetScoresCol(unsigned ColIndex);

	void SetFreqs();
	void SetFreqsCol(unsigned ColIndex);

	void SetPseudoFreqs();
	void SetPseudoFreqsCol(unsigned ColIndex);

	void AddPseudo(float w);
	void AddPseudoCol(unsigned ColIndex, float w);
	};

#endif // pssm_h
