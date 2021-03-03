#include "myutils.h"
#include "fastq.h"
#include "seqinfo.h"
#include "fastqrec.h"

double FastQ::m_Wildcard_P_e = 1.0;
byte FastQ::m_ASCII_Offset = 33;
byte FastQ::m_IntQual_Min = 0;
byte FastQ::m_IntQual_Max = 41;
byte FastQ::m_IntQualOut_Max = 41;
double FastQ::m_CharToProb[256];
double FastQ::m_IntQualToProb[256];
byte **FastQ::m_PairMatch = 0;
byte **FastQ::m_PairMismatch = 0;

static bool GetLine(FILE *f, unsigned &LineNr, string &Line)
	{
	++LineNr;
	bool b = ReadLineStdioFile(f, Line);
	if (!b)
		return false;
//	Log("Line %u '%s'\n", LineNr, Line.c_str());
	return b;
	}

void FastqRec::ToFile(FILE *f)
	{
	fprintf(f, "@%s\n", Label.c_str());
	fprintf(f, "%s\n", Seq.c_str());
	fprintf(f, "+\n");
	fprintf(f, "%s\n", Qual.c_str());
	}

bool FastqRec::FromFile(FILE *f)
	{
	string s;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(f, s);
		if (!Ok)
			return false;

		Label.clear();
		if (s.empty())
			{
		// Hack to deal with empty line at EOF
			Ok = ReadLineStdioFile(f, Label);
			if (!Ok)
				return false;
			Die("Bad fastq, expected '@'");
			}

		unsigned n = SIZE(s);
		if (n <= 1)
			Die("Empty label in FASTQ file");
		for (unsigned i = 1; i < n; ++i)
			{
			char c = s[i];
			if (c == ' ')
				c = '_';
			Label += c;
			}

		Ok = ReadLineStdioFile(f, Seq);
		if (!Ok)
			Die("Incomplete record at end-of-file");

		string Plus;
		Ok = ReadLineStdioFile(f, Plus);
		if (!Ok)
			Die("Incomplete record at end-of-file");
		if (Plus.empty() || Plus[0] != '+')
			Die("Bad fastq, expected '+'");

		Ok = ReadLineStdioFile(f, Qual);
		if (!Ok)
			Die("Incomplete fastq record at end-of-file");

		if (Seq.size() != Qual.size())
			Die("Bad fastq, seq len != qual len");
		return true;
		}
	}

void FastqRec::ToFastq(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "@%s\n", Label.c_str());
	fprintf(f, "%s\n", Seq.c_str());
	fprintf(f, "+\n");
	fprintf(f, "%s\n", Qual.c_str());
	}

void FastqRec::ToFasta(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, ">%s\n", Label.c_str());
	fprintf(f, "%s\n", Seq.c_str());
	}

void FastqRec::LogMeVerbose(bool Plus) const
	{
	unsigned L = SIZE(Seq);
	asserta(L < 1024);
	byte RCSeq[1024];
	void RevComp(const byte *Seq, unsigned L, byte *RCSeq);
	RevComp((const byte *) Seq.c_str(), L, RCSeq);

	const byte *S;

	if (Plus)
		S = (const byte *) Seq.c_str();
	else
		S = RCSeq;

	Log("");
	Log("%c>%s\n", pom(Plus), Label.c_str());
	Log("%10u  Length\n", L);
	//unsigned MinQualPos;
	//int MinQual = GetMinQual(MinQualPos);
	//Log("%10d  Min qual (pos %u)\n", MinQual, MinQualPos);
	if (L == 0)
		return;

	Log("+ %*.*s\n", L, L, Seq.c_str());
	Log("- %*.*s\n", L, L, RCSeq);
	Log("q %s\n", Qual.c_str());

	const unsigned BLOCK = 16;
	unsigned Lo = 0;
	for (;;)
		{
		if (Lo >= L)
			break;

		unsigned End = Lo + BLOCK;
		if (End > L)
			End = L;

		Log("\n");

		for (unsigned Pos = Lo; Pos < End; ++Pos)
			Log("  %c", Seq[Pos]);
		Log("\n");

		for (unsigned Pos = Lo; Pos < End; ++Pos)
			Log("  %c", Qual[Pos]);
		Log("\n");

		//for (unsigned Pos = Lo; Pos < End; ++Pos)
		//	Log(" %2d", GetIntQual(Pos));
		//Log("\n");

		Lo += BLOCK;
		}
	}

void FastqRec::ToSI(SeqInfo *SI) const
	{
	unsigned L = GetLength();
	SI->AllocL(L);
	memcpy(SI->m_SeqBuffer, (const byte *) Seq.c_str(), L);
	SI->SetLabel(Label.c_str());
	SI->m_L = L;
	}

FASTQ_FILTER FastqFilter(FastqRec &Rec)
	{
	if (optset_fastq_truncqual)
		Rec.TruncateQual(opt_fastq_truncqual);

	if (optset_fastq_stripleft)
		Rec.StripLeft(opt_fastq_stripleft);

	unsigned L = Rec.GetLength();
	if (L == 0)
		return FF_Short;

	if (optset_fastq_minlen && L < opt_fastq_minlen)
		return FF_Short;

	if (optset_fastq_trunclen)
		{
		if (L < opt_fastq_trunclen)
			return FF_Short;

		Rec.TruncateLength(opt_fastq_trunclen);
		unsigned NewL = Rec.GetLength();
		asserta(NewL == opt_fastq_trunclen);
		L = opt_fastq_trunclen;
		}

	if (optset_fastq_maxee)
		{
		double ExErr = Rec.GetExpectedErrCount(L);
		if (ExErr > opt_fastq_maxee)
			return FF_HighErr;
		}

	return FF_Good;
	}

void FastqRec::GetRevComp(FastqRec &RC) const
	{
	RC.Label = Label + string(".revcomp");
	RC.Seq.clear();
	RC.Qual.clear();

	unsigned L = SIZE(Seq);

	for (unsigned i = 0; i < L; ++i)
		{
		byte c = Seq[L-i-1];
		byte q = Qual[L-i-1];
		RC.Seq += g_CharToCompChar[c];
		RC.Qual += q;
		}
	}
