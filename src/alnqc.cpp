#include "myutils.h"
#include "alpha.h"

static uint g_BadCount;
static uint g_QueryCount;
static FILE *g_fFev;
static FILE *g_fRep;
static FILE *g_fTsv;

static void FindLongestIdenticalSkippingGaps(const string &QRow,
  const string &TRow, uint &BestLo, uint &BestHi, uint &BestN)
	{
	BestLo = UINT_MAX;
	BestHi = UINT_MAX;
	BestN = 0;

	const uint Cols = SIZE(QRow);
	asserta(SIZE(TRow) == Cols);

	uint Lo = UINT_MAX;
	uint Hi = UINT_MAX;
	uint N = 0;
	for (uint i = 0; i < Cols; ++i)
		{
		char q = QRow[i];
		char t = TRow[i];
		if (isgap(q) || isgap(t))
			continue;
		if (q == t)
			{
			if (Lo == UINT_MAX)
				{
				Lo = i;
				Hi = i;
				N = 1;
				}
			else
				{
				Hi = i;
				++N;
				}
			}
		else
			{
			if (Lo != UINT_MAX)
				{
				if (BestLo == UINT_MAX)
					{
					BestLo = Lo;
					BestHi = Hi;
					BestN = N;
					}
				else
					{
					if (N > BestN)
						{
						BestLo = Lo;
						BestHi = Hi;
						BestN = N;
						}
					}
				}
			Lo = UINT_MAX;
			Hi = UINT_MAX;
			N = 0;
			}
		}
	if (N > BestN)
		{
		BestLo = Lo;
		BestHi = Hi;
		BestN = N;
		}
	}

static void FindLongestGap(const string &QRow,
  const string &TRow, uint &GapPos, uint &GapLen)
	{
	GapPos = UINT_MAX;
	GapLen = 0;

	const uint Cols = SIZE(QRow);
	asserta(SIZE(TRow) == Cols);

	uint Pos = UINT_MAX;
	uint Len = 0;
	for (uint i = 0; i < Cols; ++i)
		{
		char q = QRow[i];
		char t = TRow[i];
		if (isgap(q) || isgap(t))
			{
			if (Pos == UINT_MAX)
				{
				Pos = i;
				Len = 1;
				}
			else
				++Len;
			}
		else
			{
			if (Len > GapLen)
				{
				GapPos = Pos;
				GapLen = Len;
				}
			Pos = UINT_MAX;
			Len = 0;
			}
		}
	if (Len > GapLen)
		{
		GapPos = Pos;
		GapLen = Len;
		}
	}

static void FindHighestIdentityWindowIgnoringGaps(const string &QRow,
  const string &TRow, uint w, uint &BestLo, uint &BestHi, uint &BestIdCount)
	{
	BestLo = 0;
	BestHi = 0;
	BestIdCount = 0;

	uint Cols = SIZE(QRow);
	asserta(SIZE(TRow) == Cols);
	if (Cols < w)
		Cols = w;

	uint Pos = UINT_MAX;
	uint Len = 0;
	uint End = Cols - w + 1;
	for (uint Lo = 0; Lo < End; ++Lo)
		{
		uint IdCount = 0;
		uint NonGapCount = 0;
		uint Hi = 0;
		for (uint Col = Lo; Col < Lo + w; ++Col)
			{
			char q = QRow[Col];
			char t = TRow[Col];
			if (q == 'X' || t == 'X')
				continue;
			if (q == '-' || t == '-')
				continue;
			++NonGapCount;
			Hi = Col;
			if (q == t)
				++IdCount;
			if (NonGapCount >= w)
				{
				for (uint Col2 = Col + 1; Col2 < Cols; ++Col2)
					{
					char q = QRow[Col2];
					char t = TRow[Col2];
					if (q == t)
						{
						Hi = Col2;
						++IdCount;
						}
					else
						break;
					}
				break;
				}
			}
		if (IdCount > BestIdCount)
			{
			BestIdCount = IdCount;
			BestLo = Lo;
			BestHi = Hi;
			}
		}
	}

static void FindLowestIdentityWindow(const string &QRow,
  const string &TRow, uint w, 
  uint &WorstLo, uint &WorstHi, uint &WorstIdCount)
	{
	WorstLo = 0;
	WorstHi = 0;
	WorstIdCount = w+1;

	const uint Cols = SIZE(QRow);
	asserta(SIZE(TRow) == Cols);

	for (uint Lo = 0; Lo + w < Cols; ++Lo)
		{
		uint IdCount = 0;
		for (uint i = 0; i < w; ++i)
			{
			uint Col = Lo + i;
			char q = QRow[Col];
			char t = TRow[Col];
			if (q == t || q == 'X' || t == 'X')
				++IdCount;
			}
		if (IdCount < WorstIdCount)
			{
			WorstLo = Lo;
			WorstHi = Lo + w - 1;
			WorstIdCount = IdCount;
			}
		}
	}

static uint GetScore(uint BestPctId, uint BestPctCols, uint WorstPctId, uint GapLen,
  uint &IdScore, uint &GapScore)
	{
	IdScore = 0;
	GapScore = 0;
	if (BestPctId == 100)
		{
		if (BestPctCols >= 75)
			{
			if (WorstPctId <= 25)
				IdScore = 75;

			if (GapLen >= 16)
				GapScore = 100;
			else if (GapLen > 8)
				GapScore = 80;
			}
		else if (BestPctCols >= 50)
			{
			if (GapLen >= 16)
				GapScore = 100;
			}
		}
	else if (BestPctId >= 98)
		{
		if (BestPctCols >= 75 && WorstPctId <= 25)
			IdScore = 75;

		if (GapLen > 16)
			GapScore = 85;
		else if (GapLen > 8)
			GapScore = 65;
		}
	else if (BestPctId >= 90)
		{
		if (GapLen > 32)
			GapScore = 75;
		}
	uint Score = max(IdScore, GapScore);
	return Score;
	}

static void QC(const string &QLabel, const string &TLabel,
  const string &QRow, const string &TRow,
  uint BestWindowLength, uint WorstWindowLength)
	{
	const uint Cols = SIZE(QRow);
	asserta(SIZE(TRow) == Cols);

	if (Cols < 50)
		{
		if (g_fFev != 0)
			{
			fprintf(g_fFev, "score=100");
			fprintf(g_fFev, "\tlabel1=%s", QLabel.c_str());
			fprintf(g_fFev, "\tlabel2=%s", TLabel.c_str());
			fprintf(g_fFev, "\tcols=%u", Cols);
			fprintf(g_fFev, "\tcomment=short_alignment.");
			fprintf(g_fFev, "\n");
			}
		if (g_fRep != 0)
			{
			fprintf(g_fRep, ">%s\n", QLabel.c_str());
			fprintf(g_fRep, ">%s\n", TLabel.c_str());
			fprintf(g_fRep, "%s\n", QRow.c_str());
			fprintf(g_fRep, "%s\n", TRow.c_str());
			fprintf(g_fRep, "Score 100, short alignment");
			fprintf(g_fRep, "\n");
			}
		return;
		}

	uint BestLo, BestHi, BestIdCount;
	FindHighestIdentityWindowIgnoringGaps(QRow, TRow, BestWindowLength,
	  BestLo, BestHi, BestIdCount);
	uint BestLength = BestHi - BestLo + 1;
	uint BestPctId = (BestIdCount*100)/BestLength;
	uint BestPctCols = (BestLength*100)/Cols;

	uint WorstLo, WorstHi, WorstIdCount;
	FindLowestIdentityWindow(QRow, TRow, WorstWindowLength,
	  WorstLo, WorstHi, WorstIdCount);
	uint WorstLength = WorstHi - WorstLo + 1;
	asserta(WorstLength == WorstWindowLength);
	uint WorstPctId = (WorstIdCount*100)/WorstWindowLength;

	uint GapPos, GapLength;
	FindLongestGap(QRow, TRow, GapPos, GapLength);

	uint IdScore, GapScore;
	uint Score = GetScore(BestPctId, BestPctCols, WorstPctId, GapLength, IdScore, GapScore);

	if (g_fTsv != 0)
		{
		static bool HdrDone = false;
		if (!HdrDone)
			{
			fprintf(g_fTsv, "V	QCPair	QC2	QCScore\n");
			HdrDone = true;
			}
		if (Score >= 50)
			{
			static uint PairCount = 0;
			++PairCount;

			fprintf(g_fTsv, "%s", QLabel.c_str());
			fprintf(g_fTsv, "\t%u", PairCount);
			fprintf(g_fTsv, "\t%s", TLabel.c_str());
			fprintf(g_fTsv, "\t%u", Score);
			fprintf(g_fTsv, "\n");

			fprintf(g_fTsv, "%s", TLabel.c_str());
			fprintf(g_fTsv, "\t%u", PairCount);
			fprintf(g_fTsv, "\t%s", QLabel.c_str());
			fprintf(g_fTsv, "\t%u", Score);
			fprintf(g_fTsv, "\n");
			}
		}

	if (g_fFev != 0)
		{
		fprintf(g_fFev, "score=%u", Score);
		fprintf(g_fFev, "\tlabel1=%s", QLabel.c_str());
		fprintf(g_fFev, "\tlabel2=%s", TLabel.c_str());
		fprintf(g_fFev, "\tbestid=%u", BestPctId);
		fprintf(g_fFev, "\tworstid=%u", WorstPctId);
		fprintf(g_fFev, "\tgap=%u", GapLength);
		fprintf(g_fFev, "\n");
		}

	if (Score <= 50)
		return;

	++g_BadCount;

	uint Ids = 0;
	for (uint i = 0; i < Cols; ++i)
		{
		char q = QRow[i];
		char t = TRow[i];
		if (q == t)
			++Ids;
		}

	if (g_fRep != 0)
		{
		string ARow;
		for (uint Col = 0; Col < Cols; ++Col)
			ARow += ' ';
		for (uint Col = WorstLo; Col <= WorstHi; ++Col)
			ARow[Col] = '!';
		for (uint Col = BestLo; Col <= BestHi; ++Col)
			ARow[Col] = '+';

		fprintf(g_fRep, "\n");
		for (uint i = 0; i < Cols; ++i)
			fprintf(g_fRep, "_");
		fprintf(g_fRep, "\n");

		fprintf(g_fRep, ">%s\n", QLabel.c_str());
		fprintf(g_fRep, ">%s\n", TLabel.c_str());
		fprintf(g_fRep, "%s\n", QRow.c_str());
		for (uint Col = 0; Col < Cols; ++Col)
			{
			char q = QRow[Col];
			char t = TRow[Col];

			if (q == t)
				fprintf(g_fRep, "|");
			else
				fprintf(g_fRep, " ");
			}
		fprintf(g_fRep, "\n");
		fprintf(g_fRep, "%s\n", TRow.c_str());
		fprintf(g_fRep, "%s\n", ARow.c_str());

		fprintf(g_fRep, "Score %u (%u,%u)", Score, IdScore, GapScore);
		fprintf(g_fRep, ", best %u%%", BestPctId);
		fprintf(g_fRep, ", worst %u%%", WorstPctId);
		fprintf(g_fRep, ", gap %u", GapLength);
		fprintf(g_fRep, ", alnid %.1f%%", GetPct(Ids, Cols));
		fprintf(g_fRep, "\n");
		}
	}

void AlnQC()
	{
	const string &UserFileName = opt_alnqc;

	uint hiw = 45;
	uint low = 16;
	if (optset_hiw)
		hiw = opt_hiw;
	if (optset_low)
		low = opt_low;

	FILE *fIn = OpenStdioFile(UserFileName);

	g_fFev = CreateStdioFile(opt_fevout);
	g_fRep = CreateStdioFile(opt_report);
	g_fTsv = CreateStdioFile(opt_tsvout);

	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(fIn, Line))
		{
		Split(Line, Fields, '\t');
		assert(SIZE(Fields) == 5);

		const string &QLabel = Fields[0];
		const string &TLabel = Fields[1];
		double PctId = StrToFloat(Fields[2].c_str());
		const string &QRow = Fields[3];
		const string &TRow = Fields[4];

		++g_QueryCount;
		if (g_QueryCount%10000 == 0)
			Progress("%u/%u bad\r", g_BadCount, g_QueryCount);
		QC(QLabel, TLabel, QRow, TRow, hiw, low);
		}
	ProgressLog("%u/%u bad\n", g_BadCount, g_QueryCount);

	CloseStdioFile(fIn);
	CloseStdioFile(g_fFev);
	CloseStdioFile(g_fRep);
	CloseStdioFile(g_fTsv);
	}
