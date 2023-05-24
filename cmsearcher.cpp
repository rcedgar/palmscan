#include "myutils.h"
#include "cmsearcher.h"

void CMSearcher::Search(PDBChain &Query)
	{
	ClearSearch();
	m_Query = &Query;
	const uint RefCount = SIZE(m_ProfDB.m_DistMxVec);
#if 1
	Log("\n>%s\n", Query.m_Label.c_str());
#endif
	for (uint RefIndex = 0; RefIndex < RefCount; ++RefIndex)
		{
		m_PS.m_DistMx = &m_ProfDB.m_DistMxVec[RefIndex];
		m_PS.m_StdDevs = &m_ProfDB.m_StdDevs;
		m_PS.Search(Query);

		uint PosA, PosB, PosC;
		double Score = m_PS.GetPSSMStarts(PosA, PosB, PosC);
#if 1
		const char *RefLabel = m_ProfDB.m_RefLabels[RefIndex].c_str();
		Log("Score = %8.3f  %4u  %4u  %4u  >%s\n",
		  Score, PosA, PosB, PosC, RefLabel);
#endif
		if (PosA != UINT_MAX && Score < m_Score)
			{
			m_RefIndex = RefIndex;
			m_Score = Score;
			m_PosA = PosA;
			m_PosB = PosB;
			m_PosC = PosC;
			}
		}
	}

uint CMSearcher::GetPSSMStarts(uint &PosA, uint &PosB, uint &PosC,
  double &Score) const
	{
	PosA = m_PosA;
	PosB = m_PosB;
	PosC = m_PosC;
	Score = m_Score;
	return m_RefIndex;
	}
