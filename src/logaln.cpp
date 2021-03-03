#include "myutils.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "seqdb.h"
#include "timing.h"

extern float **g_SubstMx;

void TrimTermGaps(const char *Path,
  unsigned &QLo, unsigned &QHi,
  unsigned &TLo, unsigned &THi,
  unsigned &ColLo, unsigned &ColHi)
	{
	ColLo = UINT_MAX;
	ColHi = UINT_MAX;
	unsigned i = 0;
	unsigned j = 0;
	QLo = 0;
	QHi = 0;
	TLo = 0;
	THi = 0;
	for (unsigned k = 0; ; ++k)
		{
		char c = Path[k];
		if (c == 0)
			return;
		if (c == 'M')
			{
			if (ColLo == UINT_MAX)
				{
				ColLo = k;
				QLo = i;
				TLo = j;
				}
			ColHi = k;
			QHi = i;
			THi = j;
			}
		if (c == 'M' || c == 'D')
			++i;
		if (c == 'M' || c == 'I')
			++j;
		}
	}

void TrimQueryTermGaps(const char *Path,
  unsigned &QLo, unsigned &QHi,
  unsigned &ColLo, unsigned &ColHi)
	{
	ColLo = UINT_MAX;
	ColHi = UINT_MAX;
	unsigned i = 0;
	unsigned j = 0;
	QLo = 0;
	QHi = 0;
	for (unsigned k = 0; ; ++k)
		{
		char c = Path[k];
		if (c == 0)
			return;
		if (c == 'M' || c == 'I')
			{
			if (ColLo == UINT_MAX)
				{
				ColLo = k;
				QLo = i;
				}
			ColHi = k;
			QHi = i;
			}
		if (c == 'M' || c == 'D')
			++i;
		}
	}

void LogAln(const byte *A, const byte *B, const char *Path, unsigned ColCount)
	{
	unsigned p = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M' || c == 'D')
			Log("%c", A[p++]);
		else
			Log("-");
		}
	Log("\n");

	unsigned pa = 0;
	unsigned pb = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M')
			{
			byte a = A[pa];
			byte b = B[pb];
			if (toupper(a) == toupper(b))
				Log("|");
			else if (g_SubstMx[a][b] > 0.0f)
				Log("+");
			else
				Log(" ");
			}
		else
			Log(" ");
		if (c == 'M' || c == 'D')
			++pa;
		if (c == 'M' || c == 'I')
			++pb;
		}
	Log("\n");

	p = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M' || c == 'I')
			Log("%c", B[p++]);
		else
			Log("-");
		}
	Log("\n");
	}

void LogAln(const byte *A, const byte *B, const char *Path)
	{
	LogAln(A, B, Path, UINT_MAX);
	}

static void LogARow(const byte *A, const char *Path,
  unsigned &i, unsigned ColLo, unsigned ColHi)
	{
	Log("%5u ", i+1);
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M' || c == 'D')
			Log("%c", A[i++]);
		else
			Log("-");
		}
	Log(" %u\n", i);
	}

static void LogBRow(const byte *B, const char *Path,
  unsigned &j, unsigned ColLo, unsigned ColHi)
	{
	Log("%5u ", j+1);
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M' || c == 'I')
			Log("%c", B[j++]);
		else
			Log("-");
		}
	Log(" %u\n", j);
	}

static void LogAnnotRow(const byte *A, const byte *B, const char *Path,
  unsigned i, unsigned j, unsigned ColLo, unsigned ColHi)
	{
	Log("%5.5s ", "");
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M')
			{
			byte a = A[i++];
			byte b = B[j++];
			if (toupper(a) == toupper(b))
				Log("|");
			else if (g_SubstMx[a][b] > 0.0f)
				Log("+");
			else
				Log(" ");
			}
		else
			{
			if (c == 'D')
				++i;
			else if (c == 'I')
				++j;
			else
				asserta(false);
			Log(" ");
			}
		}
	Log("\n");
	}

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps)
	{
	const unsigned BLOCK_SIZE = 80;
	unsigned ALo, BLo, ColLo, ColHi;
	if (StripTermGaps)
		{
		unsigned AHi_NotUsed, BHi_NotUsed;
		TrimTermGaps(Path, ALo, AHi_NotUsed, BLo, BHi_NotUsed, ColLo, ColHi);
		}
	else
		{
		ALo = 0;
		BLo = 0;
		ColLo = 0;
		ColHi = (unsigned) strlen(Path) - 1;
		}

	asserta(ColHi >= ColLo);

	unsigned i = ALo;
	unsigned j = BLo;
	unsigned ColFrom = ColLo;
	for (;;)
		{
		if (ColFrom > ColHi)
			break;
		unsigned ColTo = ColFrom + BLOCK_SIZE - 1;
		if (ColTo > ColHi)
			ColTo = ColHi;

		unsigned i0 = i;
		unsigned j0 = j;
		LogARow(A, Path, i, ColFrom, ColTo);
		LogAnnotRow(A, B, Path, i0, j0, ColFrom, ColTo);
		LogBRow(B, Path, j, ColFrom, ColTo);
		Log("\n");

		ColFrom += BLOCK_SIZE;
		}
	}

void SeqDB::WriteMSAPretty(FILE *f) const
	{
	const unsigned ColCount = GetColCount();
	const unsigned SeqCount = GetSeqCount();

	unsigned BLOCK_SIZE = 120;
	if (BLOCK_SIZE > ColCount)
		BLOCK_SIZE = ColCount;

	const unsigned BlockCount = (ColCount + BLOCK_SIZE - 1)/BLOCK_SIZE;

	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned ColLo = BlockIndex*BLOCK_SIZE;
		unsigned ColHi = ColLo + BLOCK_SIZE - 1;
		if (ColHi >= ColCount)
			ColHi = ColCount - 1;
		unsigned n = ColHi - ColLo + 1;

		fprintf(f, "\n");
		for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			const byte *Seq = (const byte *) GetSeq(SeqIndex).c_str();
			const char *Label = GetLabel(SeqIndex).c_str();

			fprintf(f, "%*.*s  ", n, n, Seq + ColLo);
			for (unsigned i = n; i < BLOCK_SIZE; ++i)
				fputc(' ', f);
			fprintf(f, "  >%s\n", Label);
			}
		}
	}
