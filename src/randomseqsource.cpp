#include "myutils.h"
#include "randomseqsource.h"
#include "seqinfo.h"
#include "alpha.h"

extern double g_AAFreqs[20];

void GetRandSeqNucleo(byte *Seq, uint L)
	{
	for (uint i = 0; i < L; ++i)
		{
		uint r = randu32()%4;
		byte c = "ACGT"[r];
		Seq[i] = c;
		}
	}

static bool g_InitDone = false;
static byte *g_RandomAminoChars;
static uint g_RandomAminoCharCount = 1000*1000;
static void Init()
	{
	if (g_InitDone)
		return;

	g_RandomAminoChars = myalloc(byte, g_RandomAminoCharCount);
	memset(g_RandomAminoChars, 0, g_RandomAminoCharCount);
	uint Accum = 0;
	vector<double> Freqs;
	for (uint i = 0; i < 20; ++i)
		{
		byte Char = g_LetterToCharAmino[i];

		double Freq = g_AAFreqs[i];
		uint N = uint(Freq*g_RandomAminoCharCount + 0.5);
		if (Accum + N >= g_RandomAminoCharCount)
			N = g_RandomAminoCharCount - Accum - 1;
		asserta(Accum + N < g_RandomAminoCharCount);
		for (uint j = 0; j < N; ++j)
			g_RandomAminoChars[Accum+j] = Char;
		Accum += N;
		Freqs.push_back(double(N)/g_RandomAminoCharCount);
		}
	for (uint i = Accum; i < g_RandomAminoCharCount; ++i)
		{
		uint r = randu32()%20;
		g_RandomAminoChars[i] = g_LetterToCharAmino[r];
		}

	for (uint i = 0; i < g_RandomAminoCharCount; ++i)
		asserta(g_IsAminoChar[g_RandomAminoChars[i]]);

	for (uint i = 0; i < 20; ++i)
		Log("%c  %6.4f  %6.4f\n",
		  g_LetterToCharAmino[i], g_AAFreqs[i], Freqs[i]);

	g_InitDone = true;
	}

byte GetRandCharAmino()
	{
	Init();
	uint r = randu32()%g_RandomAminoCharCount;
	byte c = g_RandomAminoChars[r];
	return c;
	}

void GetRandSeqAmino(byte *Seq, uint L)
	{
	for (uint i = 0; i < L; ++i)
		{
		byte c = GetRandCharAmino();
		Seq[i] = c;
		}
	}

bool RandomSeqSource::GetNextLo(SeqInfo *SI)
	{
	if (m_SeqIndex >= m_SeqCount)
		return false;

	uint L = m_MinSeqLength + randu32()%(m_MaxSeqLength - m_MinSeqLength);
	asserta(L >= m_MinSeqLength && L <= m_MaxSeqLength);
	SI->AllocL(L);
	SI->m_L = L;

	if (m_Nucleo)
		GetRandSeqNucleo(SI->m_SeqBuffer, L);
	else
		GetRandSeqAmino(SI->m_SeqBuffer, L);
	SI->m_Seq = SI->m_SeqBuffer;
	SI->m_Index = m_SeqIndex;
	++m_SeqIndex;

	string Label;
	Ps(Label, "Rand%u", m_SeqIndex);
	SI->SetLabel(Label.c_str());

	return true;
	}
