#pragma once

#include "cdinfo.h"

class CDTemplate
	{
public:
	const CDInfo *m_Info = 0;

// Required AA in motif, or 'x' if none
	vector<char> m_AnchorAAs;

// 0-based offset of required aa relative to
//   first AA in motif.
	vector<uint> m_AnchorAAOffsets;

// Minimum number of AAs in interval
//    (first aa in motif, first aa in next motif]
	vector<uint> m_MinAAsNext;

// Maximum number of AAs in interval
//    (first aa in motif, first aa in next motif]
	vector<uint> m_MaxAAsNext;

// Minimum score to consider motif match
	vector<double> m_MinScores;

public:
	void Clear()
		{
		m_Info = 0;
		m_AnchorAAs.clear();
		m_AnchorAAOffsets.clear();
		m_MinAAsNext.clear();
		m_MaxAAsNext.clear();
		m_MinScores.clear();
		}

	void Init(
	  const CDInfo &Info,
	  const vector<char> &AnchorAAs,
	  const vector<uint> &AnchorAAOffsets,
	  const vector<uint> &MinAAsNext,
	  const vector<uint> &MaxAAsNext,
	  const vector<double> &MinScores);

	uint GetMinAAsNext(uint MotifIndex) const;
	uint GetMaxAAsNext(uint MotifIndex) const;
	uint GetMinAAsSeqEnd(uint MotifIndex) const;
	uint GetMinAAsSeqStart(uint MotifIndex) const;
	};
