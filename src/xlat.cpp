#include "myutils.h"
#include "alpha.h"

unsigned g_MaxLaa = 0;
byte *g_AASeqs[7];

// NtPos is 0-based relative to start of seq or rev-comp'd seq.
unsigned AAPosToNtPos(unsigned AAPos, unsigned Lnt, int Frame)
	{
	if (AAPos == UINT_MAX)
		return UINT_MAX;
	asserta(Frame >= -3 && Frame <= 3 && Frame != 0);
	asserta(AAPos < Lnt/3);
	if (Frame > 0)
		{
		unsigned Offset = AAPos*3 + Frame - 1;
		if (Offset >= Lnt)
			Die("AAPosToNtPos(AAPos=%u, Lnt=%u, Frame=%+d), Offset=%u",
			  AAPos, Lnt, Frame, Offset);
		return Offset;
		}
	else
		{
		unsigned Offset = AAPos*3 + (-Frame) - 1;
		if (Offset >= Lnt)
			Die("AAPosToNtPos(AAPos=%u, Lnt=%u, Frame=%+d), Offset=%u",
			  AAPos, Lnt, Frame, Offset);
		return Offset;
		}
	}

unsigned XlatX(const byte *nt, unsigned Lnt, int Frame, byte *aa)
	{
	unsigned Laa = 0;
	if (Frame > 0)
		{
		for (int Pos = Frame-1; Pos + 2 < (int) Lnt; Pos += 3)
			{
			byte c1 = nt[Pos];
			byte c2 = nt[Pos+1];
			byte c3 = nt[Pos+2];

			byte x1 = g_CharToLetterNucleo[c1];
			byte x2 = g_CharToLetterNucleo[c2];
			byte x3 = g_CharToLetterNucleo[c3];

			unsigned Word = 16*x1 + 4*x2 + x3;
			byte AminoChar = 'X';
			if (Word < 64)
				AminoChar = g_CodonWordToAminoChar[Word];
			aa[Laa++] = AminoChar;
			}
		}
	else
		{
		for (int Pos = Lnt + Frame; Pos >= 2; Pos -= 3)
			{
			byte c1 = nt[Pos];
			byte c2 = nt[Pos-1];
			byte c3 = nt[Pos-2];

			byte x1 = g_CharToCompLetter[c1];
			byte x2 = g_CharToCompLetter[c2];
			byte x3 = g_CharToCompLetter[c3];

			unsigned Word = 16*x1 + 4*x2 + x3;
			byte AminoChar = 'X';
			if (Word < 64)
				AminoChar = g_CodonWordToAminoChar[Word];
			aa[Laa++] = AminoChar;
			}
		}
	return Laa;
	}

void AllocXlat(unsigned Lnt)
	{
	if ((Lnt+2)/3 <= g_MaxLaa)
		return;

	if (g_AASeqs[0] != 0)
		{
		for (int Frame = -3; Frame <= +3; ++Frame)
			{
			if (Frame == 0)
				continue;
			myfree(g_AASeqs[Frame+3]);
			}
		}

	g_MaxLaa = Lnt;
	for (int Frame = -3; Frame <= +3; ++Frame)
		{
		if (Frame == 0)
			continue;
		g_AASeqs[Frame+3] = myalloc(byte, g_MaxLaa);
		}
	}

void PrintFastaSegNt(FILE *f, const char *nt, unsigned Lnt, unsigned NtPos, 
  unsigned CodonCount, int Frame)
	{
	assert(NtPos != UINT_MAX);

	unsigned Col = 0;
	if (Frame > 0)
		{
		for (unsigned CodonIndex = 0; CodonIndex < CodonCount; ++CodonIndex)
			{
			Col += 3;
			if (Col >= 80)
				{
				fputc('\n', f);
				Col = 3;
				}
			fputc(nt[NtPos++], f);
			fputc(nt[NtPos++], f);
			fputc(nt[NtPos++], f);
			}
		}
	else
		{
		asserta(NtPos < Lnt);
		unsigned Pos = Lnt - NtPos - 1;
		for (unsigned CodonIndex = 0; CodonIndex < CodonCount; ++CodonIndex)
			{
			Col += 3;
			if (Col >= 80)
				{
				fputc('\n', f);
				Col = 3;
				}
			asserta(Pos >= 2 && Pos < Lnt);
			fputc(g_CharToCompChar[nt[Pos--]], f);
			fputc(g_CharToCompChar[nt[Pos--]], f);
			fputc(g_CharToCompChar[nt[Pos--]], f);
			}
		}
	asserta(Col > 0);
	fputc('\n', f);
	}

void GetNtSeg(const char *nt, unsigned Lnt, unsigned NtPos, 
  unsigned CodonCount, int Frame, string &s)
	{
	assert(NtPos != UINT_MAX);
	s.clear();

	if (Frame > 0)
		{
		for (unsigned CodonIndex = 0; CodonIndex < CodonCount; ++CodonIndex)
			{
			s += nt[NtPos++];
			s += nt[NtPos++];
			s += nt[NtPos++];
			}
		}
	else
		{
		asserta(NtPos < Lnt);
		unsigned Pos = Lnt - NtPos - 1;
		for (unsigned CodonIndex = 0; CodonIndex < CodonCount; ++CodonIndex)
			{
			asserta(Pos >= 2 && Pos < Lnt);
			s += g_CharToCompChar[nt[Pos--]];
			s += g_CharToCompChar[nt[Pos--]];
			s += g_CharToCompChar[nt[Pos--]];
			}
		}
	}

const char *XlatStr(const char *nt, unsigned Lnt, unsigned NtPos, 
  unsigned CodonCount, int Frame, string &s)
	{
	s.clear();
	if (NtPos == UINT_MAX)
		{
		s = "*";
		return s.c_str();
		}

	if (Frame > 0)
		{
		for (unsigned CodonIndex = 0; CodonIndex < CodonCount; ++CodonIndex)
			{
			byte c1 = nt[NtPos++];
			byte c2 = nt[NtPos++];
			byte c3 = nt[NtPos++];

			byte x1 = g_CharToLetterNucleo[c1];
			byte x2 = g_CharToLetterNucleo[c2];
			byte x3 = g_CharToLetterNucleo[c3];

			unsigned Word = 16*x1 + 4*x2 + x3;
			byte AminoChar = 'X';
			if (Word < 64)
				AminoChar = g_CodonWordToAminoChar[Word];
			s.push_back(AminoChar);
			}
		}
	else
		{
		asserta(NtPos < Lnt);
		unsigned Pos = Lnt - NtPos - 1;
		for (unsigned CodonIndex = 0; CodonIndex < CodonCount; ++CodonIndex)
			{
			asserta(Pos >= 2);
			byte c1 = nt[Pos--];
			byte c2 = nt[Pos--];
			byte c3 = nt[Pos--];

			byte x1 = g_CharToCompLetter[c1];
			byte x2 = g_CharToCompLetter[c2];
			byte x3 = g_CharToCompLetter[c3];

			unsigned Word = 16*x1 + 4*x2 + x3;
			byte AminoChar = 'X';
			if (Word < 64)
				AminoChar = g_CodonWordToAminoChar[Word];
			s.push_back(AminoChar);
			}
		}
	return s.c_str();
	}

void Translate(const string &NtSeq, int Frame, string &AASeq)
	{
	AASeq.clear();

	int Off = abs(Frame) - 1;
	if (Off < 0 || Off > 2)
		Die("Frame=%d, Off=%d", Frame, Off);
	const byte *Q = (const byte *) NtSeq.c_str();
	const uint QL = SIZE(NtSeq);

	if (Frame > 0)
		{
		for (uint Pos = Off; Pos + 2 < QL; Pos += 3)
			{
			char a = GetAminoCharFrom3NucChars(Q[Pos], Q[Pos+1], Q[Pos+2]);
			if (a == '?' || a == '*')
				a = 'X';
			AASeq += a;
			}
		}
	else
		{
		for (int Pos = int(QL) - int(Off) - 1; Pos - 2 >= 0; Pos -= 3)
			{
			char c1 = g_CharToCompChar[Q[Pos]];
			char c2 = g_CharToCompChar[Q[Pos-1]];
			char c3 = g_CharToCompChar[Q[Pos-2]];
			char a = GetAminoCharFrom3NucChars(c1, c2, c3);
			if (a == '?')
				a = 'X';
			AASeq += a;
			}
		}
	}

void SeqToUpper(string &Seq)
	{
	const uint QL = SIZE(Seq);
	for (uint i = 0; i < QL; ++i)
		Seq[i] = toupper(Seq[i]);
	}

bool GetIsNucleo(const string &Seq)
	{
	uint NucCount = 0;
	uint OtherCount = 0;
	const uint QL = SIZE(Seq);
	for (uint i = 0; i < QL; ++i)
		{
		char c = toupper(Seq[i]);
		if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')
			++NucCount;
		else
			++OtherCount;
		}

	bool IsNuc = (NucCount > 2*OtherCount);
	return IsNuc;
	}
