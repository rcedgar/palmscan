#pragma once

#include "pssm.h"
#include "abcxyz.h"

const uint MPL = AL + BL + CL;

class MotifProfile
	{
public:
	string m_Name;
	vector<vector<float> > m_FreqVec;
	uint m_InputSeqIndex = UINT_MAX;
	uint m_PosA = UINT_MAX;
	uint m_PosB = UINT_MAX;
	uint m_PosC = UINT_MAX;

public:
	MotifProfile()
		{
		m_FreqVec.resize(MPL);
		for (uint i = 0; i < MPL; ++i)
			m_FreqVec[i].resize(20, 0.0);
		m_InputSeqIndex = UINT_MAX;
		m_PosA = UINT_MAX;
		m_PosB = UINT_MAX;
		m_PosC = UINT_MAX;
		}

	void ValidateFreqs() const
		{
		asserta(m_FreqVec.size() == MPL);
		for (uint i = 0; i < MPL; ++i)
			{
			const vector<float> &v = m_FreqVec[i];
			asserta(v.size() == 20);
			double Sum = 0;
			for (uint j = 0; j < 20; ++j)
				Sum += v[j];
			asserta(feq(Sum, 1.0));
			}
		}

	void Clear()
		{
		asserta(m_FreqVec.size() == MPL);
		for (uint i = 0; i < MPL; ++i)
			{
			for (uint j = 0; j < 20; ++j)
				m_FreqVec[i][j] = 0;
			}
		m_PosA = UINT_MAX;
		m_PosB = UINT_MAX;
		m_PosC = UINT_MAX;
		}

	void FromSeq(uint SeqIndex, uint PosA, uint PosB, uint PosC,
	  const string &Seq);
	void FromSeqs(const vector<string> &Seqs);
	void FromPSSMs(const PSSM &PA, const PSSM &PB, const PSSM &PC);

	void LogMe() const;

	char PosToMotif(uint i) const
		{
		if (i <= AL)
			return 'A';
		if (i <= AL + BL)
			return 'B';
		return 'C';
		}

	uint PosToOffset(uint i) const
		{
		if (i < AL)
			return i;
		if (i < AL + BL)
			return i - AL;
		return i - (AL + BL);
		}

	void GetConsSeq(string &Logo) const;
	void GetLogo(string &Logo) const;
	void GetMaxLetter(uint i, uint &Letter, float &Freq) const;

public:
	static void GetLettersFromSeq(const string &Seq,
	  vector<uint> &Letters);
	};
