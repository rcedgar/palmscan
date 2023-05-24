#pragma once

#include "seqdb.h"
#include "rdrpmodel.h"

static const uint LA = 12;
static const uint LB = 14;
static const uint LC = 8;

class MSAQC
	{
public:
	RdRpModel m_Mod;
	string m_GroupName;
	uint m_FoundCount;
	const SeqDB *m_MSA = 0;
	uint m_ColCount = 0;
	uint m_SeqCount = 0;
	vector<string> m_AlignedSeqs;
	vector<string> m_UnalignedSeqs;
	vector<vector<uint> > m_UnalignedPosToAlignedColVec;
	vector<vector<uint> > m_AlignedColToUnalignedPosVec;

	vector<uint> m_Ls;
	vector<vector<uint> > m_SeqIndexToMotifStartPosVec;
	vector<vector<uint> > m_ColsVec;
	vector<string> m_SeqIndexToGroupName;

public:
	void Init(const SeqDB &MSA, const string &GroupName);
	uint GetMotifIndex(char c) const;
	uint GetTopIx(const vector<uint> &v, uint MinCount = UINT_MAX) const;
	char mtoc(uint m) const { asserta(m < 3); return "ABC"[m]; }

	void ReportCols(FILE *f) const;
	bool Approve(string &s) const;

	bool CheckIncreasingOrder(const vector<uint> &v) const;
	uint GetCorrectlyAlignedLetterCount(uint SeqIndex) const;
	uint GetCorrectlyAlignedLetterCount1(uint SeqIndex, uint m) const;
	bool MotifCorrectlyAligned(uint SeqIndex, uint m, string &s) const;
	bool MotifsCorrectlyAligned(uint SeqIndex, string &s) const;
	uint GetCorrectSeqCount() const;
	uint GetCorrectLetterCount() const;
	uint GetTotalLetterCount() const;
	uint WriteAln(FILE *f);
	void WriteAlnHdr(FILE *f, 
	  const vector<uint> &LoCols, const vector<uint> &HiCols,
	  const vector<bool> &UpperCols) const;
	void WriteAln1(FILE *f, uint SeqIndex, const string &s,
	  const vector<uint> &LoCols, const vector<uint> &HiCols,
	  const vector<bool> &UpperCols) const;
	void WritePalmscanReport(FILE *f, uint SeqIndex, const string &s);

public:
	static void GetColMaps(const string &UnalignedSeq, const string &AlignedSeq,
	  vector<uint> &UnalignedPosToAlignedCol,
	  vector<uint> &AlignedColToUnalginedPos);
	};
