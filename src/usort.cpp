#include "myutils.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "alpha.h"
#include "sort.h"

static byte *g_QueryHasWord;
static unsigned g_WordCount;

unsigned GetWord(const byte *Seq)
	{
	unsigned Word = 0;
	const byte *Front = Seq;
	for (unsigned i = 0; i < opt_usort_w; ++i)
		{
		unsigned Letter = g_CharToLetterNucleo[*Front++];
		Word = (Word*4) + Letter;
		}
	return Word;
	}

static void SetQuery(const SeqInfo &Query)
	{
	if (g_QueryHasWord == 0)
		{
		g_WordCount = 4;
		for (unsigned i = 1; i < opt_usort_w; ++i)
			g_WordCount *= 4;

		g_QueryHasWord = myalloc(byte, g_WordCount);
		}

	memset(g_QueryHasWord, 0, g_WordCount);

	if (Query.m_L <= opt_usort_w)
		return;

	const unsigned L = Query.m_L - opt_usort_w + 1;
	const byte *Seq = Query.m_Seq;
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned Word = GetWord(Seq++);
		if (Word < g_WordCount)
			g_QueryHasWord[Word] = 1;
		}
	}

static unsigned GetUniqueWordsInCommon(const byte *Target, unsigned TL)
	{
	if (TL <= opt_usort_w)
		return 0;

	unsigned Count = 0;
	const unsigned L = TL - opt_usort_w + 1;
	const byte *Seq = Target;
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned Word = GetWord(Seq++);
		if (Word < g_WordCount)
			if (g_QueryHasWord[Word])
				++Count;
		}
	return Count;
	}

void USort(const SeqInfo &Query, const SeqDB &DB, vector<unsigned> &WordCounts, 
  vector<unsigned> &Order)
	{
	WordCounts.clear();
	Order.clear();

	SetQuery(Query);

	const unsigned SeqCount = DB.GetSeqCount();
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const byte *Target = (const byte *) DB.GetSeq(SeqIndex).c_str();
		unsigned TL = DB.GetSeqLength(SeqIndex);
		unsigned WordCount = GetUniqueWordsInCommon(Target, TL);
		WordCounts.push_back(WordCount);
		}
	SortDescending(WordCounts, Order);
	}
