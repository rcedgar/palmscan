#pragma once

#include "rphit.h"
#include "rdrpsearcher.h"

class PalmHit
	{
public:
	RPHit m_A;
	RPHit m_B;
	RPHit m_C;
	float m_Score = -1;
	uint m_GroupIndex = UINT_MAX;
	bool m_Permuted = false;

public:
	void Clear()
		{
		m_A.Clear();
		m_B.Clear();
		m_C.Clear();
		m_Score = -1;
		m_GroupIndex = UINT_MAX;
		m_Permuted = false;
		}

	const RPHit &GetHit(uint MotifIndex) const
		{
		switch (MotifIndex)
			{
		case MOTIF_A: return m_A;
		case MOTIF_B: return m_B;
		case MOTIF_C: return m_C;
			}
		asserta(false);
		return *(RPHit *) 0;
		}

	void GetSpan(uint &QLo, uint &QHi) const
		{
		QLo = UINT_MAX;
		QHi = UINT_MAX;
		if (m_Score <= 0)
			return;
		if (m_Permuted)
			{
			QLo = m_C.m_QPos;
			QHi = m_B.GetQHi();
			}
		else
			{
			QLo = m_A.m_QPos;
			QHi = m_C.GetQHi();
			}
		}
	};
