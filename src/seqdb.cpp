#include "myutils.h"
#include "seqdb.h"
#include "sfasta.h"
#include "alpha.h"

unsigned SeqDB::AddSeq(const string &Label, const string &Seq)
	{
	unsigned SeqIndex = SIZE(m_Seqs);
	unsigned L = SIZE(Seq);
	if (SeqIndex == 0)
		{
		m_ColCount = L;
		m_IsAligned = true;
		}
	else if (L != m_ColCount)
		m_IsAligned = false;
	m_Labels.push_back(Label);
	m_Seqs.push_back(Seq);
	return SeqIndex;
	}

const string &SeqDB::GetSeq(unsigned SeqIndex) const
	{
	assert(SeqIndex < SIZE(m_Seqs));
	return m_Seqs[SeqIndex];
	}

const string &SeqDB::GetLabel(unsigned SeqIndex) const
	{
	assert(SeqIndex < SIZE(m_Labels));
	return m_Labels[SeqIndex];
	}

unsigned SeqDB::GetSeqLength(unsigned SeqIndex) const
	{
	assert(SeqIndex < SIZE(m_Seqs));
	return SIZE(m_Seqs[SeqIndex]);
	}

bool SeqDB::IsAligned() const
	{
	return m_IsAligned;
	}

unsigned SeqDB::GetColCount() const
	{
	asserta(m_IsAligned);
	return m_ColCount;
	}

bool SeqDB::GetIsNucleo()
	{
	if (!m_IsNucleoSet)
		SetIsNucleo();
	return m_IsNucleo;
	}

void SeqDB::SetIsNucleo()
	{
	if (m_IsNucleoSet)
		return;

	const unsigned SeqCount = GetSeqCount();
	if (SeqCount == 0)
		{
		m_IsNucleo = false;
		m_IsNucleoSet = true;
		return;
		}

	unsigned N = 0;
	unsigned i = 0;
	for (;;)
		{
		unsigned SeqIndex = unsigned(rand()%SeqCount);
		const string &Seq = GetSeq(SeqIndex);
		unsigned L = GetSeqLength(SeqIndex);
		const unsigned Pos = unsigned(rand()%L);
		byte c = Seq[Pos];
		if (isgap(c))
			continue;
		++i;
		if (i >= 100)
			break;

		if (g_IsNucleoChar[c])
			++N;
		}
	m_IsNucleo = (N > 80);
	m_IsNucleoSet = true;
	}

void SeqDB::FromFasta(const string &FileName)
	{
	SFasta SF;
	SF.Open(FileName);

	m_IsAligned = false;
	for (;;)
		{
		const char *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;
		const string Label = SF.GetLabel();
		const unsigned L = SF.GetSeqLength();
		string s;
		for (unsigned i = 0; i < L; ++i)
			s.push_back(Seq[i]);
		AddSeq(Label, s);
		}
	}

void SeqDB::WritePretty(FILE *f) const
	{
	if (f == 0)
		return;

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
			const char *Seq = GetSeq(SeqIndex).c_str();
			const char *Label = GetLabel(SeqIndex).c_str();

			fprintf(f, "%*.*s  ", n, n, Seq + ColLo);
			for (unsigned i = n; i < BLOCK_SIZE; ++i)
				fputc(' ', f);
			fprintf(f, "  >%s\n", Label);
			}
		}
	}

void SeqDB::LogMe() const
	{
	WritePretty(g_fLog);
	}
