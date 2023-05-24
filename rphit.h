#pragma once

const uint MOTIF_A = 0;
const uint MOTIF_B = 1;
const uint MOTIF_C = 2;

class RPHit
	{
public:
	const PSSM *m_PSSM = 0;
	uint m_QPos = 0;
	float m_Score = 0;

public:
	void Init(const PSSM &P, uint QPos, float Score)
		{
		m_PSSM = &P;
		m_QPos = QPos;
		m_Score = Score;
		uint m = P.m_ColCount;
		}

	void Clear()
		{
		m_PSSM = 0;
		m_QPos = 0;
		m_Score = 0;
		}

	void SetNoHit()
		{
		m_QPos = UINT_MAX;
		m_Score = -1.0f;
		}

	uint GetQHi() const
		{
		asserta(m_PSSM != 0);
		return m_QPos + m_PSSM->GetColCount() - 1;
		}

	void Copy(const RPHit &rhs)
		{
		m_PSSM = rhs.m_PSSM;
		m_QPos = rhs.m_QPos;
		m_Score = rhs.m_Score;
		}
	};
