#pragma once

class RPHit
	{
public:
	const PSSM *m_PSSM = 0;
	uint m_QPos = 0;
	string m_Path;
	float m_Score = 0;

public:
	void SetUngapped(const PSSM &P, uint QPos, float Score)
		{
		m_PSSM = &P;
		m_QPos = QPos;
		m_Score = Score;
		uint m = P.m_ColCount;
		m_Path.clear();
		for (uint i = 0; i < m; ++i)
			m_Path += 'M';
		}

	void SetGapped(const PSSM &P, uint QPos, float Score, PathInfo *PI,
	  uint ColLo, uint ColHi)
		{
		m_PSSM = &P;
		m_QPos = QPos;
		m_Score = Score;
		uint m = P.m_ColCount;
		m_Path.clear();
		const char *PIPath = PI->GetPath();
		for (uint Col = ColLo; Col <= ColHi; ++Col)
			m_Path.push_back(PIPath[Col]);
		}

	uint GetGapCount() const
		{
		uint n = 0;
		for (uint i = 0; i < SIZE(m_Path); ++i)
			if (m_Path[i] != 'M')
				++n;
		return n;
		}

	uint GetQuerySegLength() const
		{
		uint n = 0;
		for (uint i = 0; i < SIZE(m_Path); ++i)
			if (m_Path[i] == 'M' || m_Path[i] == 'D')
				++n;
		return n;
		}
	};
