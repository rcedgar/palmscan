#ifndef fastqrec_h
#define fastqrec_h

#include "fastq.h"

class SeqInfo;

class FastqRec
	{
public:
	string Label;
	string Seq;
	string Qual;

public:
	bool FromFile(FILE *f);
	void ToFile(FILE *f);
	void ToSI(SeqInfo *SI) const;
	void GetRevComp(FastqRec &RC) const;

	void TruncateLength(unsigned L)
		{
		if (L > SIZE(Seq))
			return;
		Seq = Seq.substr(0, L);
		Qual = Qual.substr(0, L);
 		asserta(SIZE(Seq) == L);
		asserta(SIZE(Qual) == L);
		}

	void StripLeft(unsigned n)
		{
		unsigned L = SIZE(Seq);
		if (L < n)
			{
			Seq.clear();
			Qual.clear();
			return;
			}
		Seq = Seq.substr(n, L-n);
		Qual = Qual.substr(n, L-n);
 		asserta(SIZE(Seq) == L-n);
		asserta(SIZE(Qual) == L-n);
		}

	unsigned GetLength() const
		{
		asserta(SIZE(Seq) == SIZE(Qual));
		return SIZE(Seq);
		}

	static char IntQualToChar(int IntQual)
		{
		return FastQ::IntQualToChar(IntQual);
		}

	static char IntQualToCharOut(int IntQual)
		{
		return FastQ::IntQualToCharOut(IntQual);
		}

	int CharToIntQual(char c) const
		{
		return FastQ::CharToIntQual(c);
		}

	double CharToProb(char c) const
		{
		return FastQ::CharToProb(c);
		}

	double GetProb(unsigned Pos) const
		{
		return FastQ::CharToProb(Qual[Pos]);
		}

	bool TruncateToQual(unsigned MinQ)
		{
		unsigned Pos = FindQualLE(MinQ);
		if (Pos == UINT_MAX)
			return false;
		Seq.resize(Pos);
		Qual.resize(Pos);
		return true;
		}

	unsigned FindQualLE(unsigned MinQ) const
		{
		const unsigned L = GetLength();
		for (unsigned Pos = 0; Pos < L; ++Pos)
			{
			unsigned Q = GetIntQual(Pos);
			if (Q <= MinQ)
				return Pos;
			}
		return UINT_MAX;
		}

	unsigned FindQualLT(unsigned MinQ) const
		{
		const unsigned L = GetLength();
		for (unsigned Pos = 0; Pos < L; ++Pos)
			{
			unsigned Q = GetIntQual(Pos);
			if (Q < MinQ)
				return Pos;
			}
		return UINT_MAX;
		}

	unsigned GetIntQual(unsigned Pos) const
		{
		assert(Pos < SIZE(Qual));
		char c = Qual[Pos];
		return CharToIntQual(c);
		}

	void UpdateQualCharCounts(double f[256])
		{
		const unsigned L = GetLength();
		for (unsigned Pos = 0; Pos < L; ++Pos)
			{
			unsigned char c = (unsigned char) Qual[Pos];
			++f[c];
			}
		}

	void UpdateIntQualCounts(double f[256])
		{
		const unsigned L = GetLength();
		for (unsigned Pos = 0; Pos < L; ++Pos)
			{
			unsigned Qual = GetIntQual(Pos);
			asserta(Qual < 256);
			++f[Qual];
			}
		}

	void UpdateIntQualCountsPos(double **f)
		{
		const unsigned L = GetLength();
		asserta(L <= MAX_FASTQ_LEN);
		for (unsigned Pos = 0; Pos < L; ++Pos)
			{
			unsigned char c = (unsigned char) Qual[Pos];
			++f[Pos][c];
			}
		}

	unsigned GetTruncLengthQ(unsigned IntQual)
		{
		unsigned L = GetLength();
		for (unsigned i = 0; i < L; ++i)
			{
			unsigned q = GetIntQual(i);
			if (q <= IntQual)
				return i;
			}
		return L;
		}

	bool HasN() const
		{
		size_t n = Seq.find('N');
		if (n == string::npos)
			return false;
		return true;
		}

	unsigned GetMinQual() const
		{
		unsigned Pos;
		return GetMinQual(Pos);
		}

	double GetAvgQual() const
		{
		unsigned L = GetLength();
		if (L == 0)
			return 0.0;

		unsigned Sum = 0;
		for (unsigned i = 0; i < L; ++i)
			{
			unsigned q = GetIntQual(i);
			Sum += q;
			}
		return double(Sum)/double(L);
		}

	void TruncateQual(int Qual)
		{
		unsigned L = GetLength();
		for (unsigned i = 0; i < L; ++i)
			{
			int q = GetIntQual(i);
			if (q <= Qual)
				{
				TruncateLength(i);
				return;
				}
			}
		}

	unsigned GetMinQual(unsigned &Pos) const
		{
		unsigned MinQual = 256;
		unsigned L = GetLength();
		Pos = UINT_MAX;
		for (unsigned i = 0; i < L; ++i)
			{
			unsigned q = GetIntQual(i);
			if (q < MinQual)
				{
				MinQual = q;
				Pos = i;
				}
			}
		return MinQual;
		}

	double GetExpectedErrCount(unsigned Length) const
		{
		unsigned L = GetLength();
		asserta(Length <= L);
		double SumP = 0.0;
		for (unsigned i = 0; i < Length; ++i)
			{
			double P = GetProb(i);
			SumP += P;
			}
		return SumP;
		}

	double GetEE() const
		{
		unsigned L = GetLength();
		double SumP = 0.0;
		for (unsigned i = 0; i < L; ++i)
			{
			double P = GetProb(i);
			SumP += P;
			}
		return SumP;
		}

	double GetExpectedErrRate() const
		{
		unsigned L = GetLength();
		if (L == 0)
			return 0.0;
		double SumP = 0.0;
		for (unsigned i = 0; i < L; ++i)
			{
			double P = GetProb(i);
			SumP += P;
			}
		return SumP/L;
		}

	double GetExpectedErrCountLoHi(unsigned Lo, unsigned Hi) const
		{
		asserta(Lo <= Hi);
		unsigned L = GetLength();
		asserta(Hi < L);
		double SumP = 0.0;
		for (unsigned i = Lo; i <= Hi; ++i)
			{
			double P = GetProb(i);
			SumP += P;
			}
		return SumP;
		}

	void UpdateErrCounts(double *PosToTotalEE, double *PosToSumP, double *PosToCount) const
		{
		unsigned L = GetLength();
		double SumP = 0.0;
		for (unsigned i = 0; i < L; ++i)
			{
			double P = GetProb(i);
			SumP += P;
			PosToSumP[i] += P;
			PosToTotalEE[i] += SumP;
			PosToCount[i] += 1;
			}
		}

	void UpdateSumLs(const vector<double> &Es, vector<vector<unsigned> > &SumLs)
		{
		unsigned N = SIZE(Es);
		vector<unsigned> Ls(N, 0);

		unsigned L = GetLength();
		double EE = 0.0;
		for (unsigned Pos = 0; Pos < L; ++Pos)
			{
			double P = GetProb(Pos);
			EE += P;
			for (unsigned k = 0; k < N; ++k)
				if (EE <= Es[k])
					Ls[k] = Pos+1;
			}

		for (unsigned k = 0; k < N; ++k)
			{
			unsigned L = Ls[k];
			asserta(L < SIZE(SumLs[k]));
			SumLs[k][L] += 1;
			}
		}

	void LogMeVerbose(bool Plus) const;
	void ToFastq(FILE *f) const;
	void ToFasta(FILE *f) const;
	};

FASTQ_FILTER FastqFilter(FastqRec &Rec);

#endif // fastqrec_h
