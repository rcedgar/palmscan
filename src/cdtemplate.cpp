#include "myutils.h"
#include "cdtemplate.h"

void CDTemplate::Init(
  const CDInfo &Info,
  const vector<char> &AnchorAAs,
  const vector<uint> &AnchorAAOffsets,
  const vector<uint> &MinAAsNext,
  const vector<uint> &MaxAAsNext,
  const vector<double> &MinScores)
	{
	Clear();
	m_Info = &Info;
	m_AnchorAAs = AnchorAAs;
	m_AnchorAAOffsets = AnchorAAOffsets;
	m_MinAAsNext = MinAAsNext;
	m_MaxAAsNext = MaxAAsNext;
	m_MinScores = MinScores;
	}

// Minimum number of amino acids in interval
//    [MotifPos=first aa in motif, first aa in next motif].
uint CDTemplate::GetMinAAsNext(uint MotifIndex) const
	{
	asserta(MotifIndex < SIZE(m_MinAAsNext));
	uint n = m_MinAAsNext[MotifIndex];
	return n;
	}

// Maximum number of amino acids in interval
//    [MotifPos=first aa in motif, first aa in next motif].
uint CDTemplate::GetMaxAAsNext(uint MotifIndex) const
	{
	asserta(MotifIndex < SIZE(m_MinAAsNext));
	uint n = m_MaxAAsNext[MotifIndex];
	return n;
	}

// Minimum number of amino acids in interval
//   [first aa in motif, last aa in sequence].
uint CDTemplate::GetMinAAsSeqEnd(uint MotifIndex) const
	{
	const uint MotifCount = m_Info->GetMotifCount();
	uint Sum = 0;
	for (uint i = MotifIndex; i < MotifCount; ++i)
		{
		uint n = m_MinAAsNext[i];
		Sum += n + 1;
		}
	return Sum;
	}

// Minimum number of amino acids in interval
//   [first aa in sequence, first aa in motif].
uint CDTemplate::GetMinAAsSeqStart(uint MotifIndex) const
	{
	const uint MotifCount = m_Info->GetMotifCount();
	uint Sum = 0;
	for (uint i = 0; i < MotifIndex; ++i)
		{
		uint n = m_MinAAsNext[i];
		Sum += n + 1;
		}
	return Sum;
	}
