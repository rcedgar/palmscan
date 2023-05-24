#pragma once

#include <map>
#include "seqdb.h"

class PAMerger
	{
public:

#define F(x)	string x;
#include "pastrs.h"

#define F(x)	uint x;
#include "paints.h"

#define F(x)	double x;
#include "padbls.h"

public:
	const string *m_Line = 0;
	string m_Label;
	map<string, string> m_NameToValue;

	string m_MotifA;
	string m_MotifB;
	string m_MotifC;

	uint m_MotifPosA;
	uint m_MotifPosB;
	uint m_MotifPosC;

	uint m_FullTrimLo;
	uint m_FullTrimHi;
	string m_FullTrimMethod;

	uint m_PPLo;
	uint m_PPHi;

public:
	void Reset()
		{
#define F(x)	x.clear();
#include "pastrs.h"

#define F(x)	x = UINT_MAX;
#include "paints.h"

#define F(x)	x = DBL_MAX;
#include "padbls.h"

		m_Line = 0;
		m_Label.clear();
		m_NameToValue.clear();

		m_MotifA.clear();
		m_MotifB.clear();
		m_MotifC.clear();

		m_MotifPosA = UINT_MAX;
		m_MotifPosB = UINT_MAX;
		m_MotifPosC = UINT_MAX;
		
		m_FullTrimLo = UINT_MAX;
		m_FullTrimHi = UINT_MAX;

		m_PPLo = UINT_MAX;
		m_PPHi = UINT_MAX;
		}

	void LogMe() const;
	void ToFev(FILE *f) const;
	void ToTsv(FILE *f) const;
	void ToFasta(FILE *f) const;
	void FromLine(const string &Line);
	void SetInt(const string &Name, uint &Value);
	void SetDbl(const string &Name, double &Value);
	void SetStr(const string &Name, string &Value);
	void SetMotifs();
	void SetMotif(const string &Seq1, const string &Seq2,
	  const string &Seq3, string &Motif);
	void SetMotif(const vector<string> &Motifs,
	  string &Motif) const;
	void SetMotifPos(const string &MotifSeq, uint &Pos) const;
	bool MotifHMMPermuted() const;
	const char *GetDGD(string &DGD) const;
	const char *GetGDD(string &GDD) const;
	void SetFullTrim();
	void SetPP();
	char GetGate() const;
	uint GetLoMotifPos() const;
	uint GetHiMotifPos() const;
	uint GetMinHitLo() const;
	uint GetMaxHitHi() const;
	double GetPSSME() const;
	double GetRdRpScore() const;
	double GetBestPlusE() const;
	double GetBestMinusE() const;
	double GetBestPlusLE() const;
	double GetBestMinusLE() const;
	bool HMMMinusWins() const;

public:
	static void LogFields();
	static bool PSSMGroupIsRT(const string &Group);
	static bool MotifHMMIsRT(const string &HMM);
	static double GetLE(double E);
	};
