#ifndef randomseqsource_h
#define randomseqsource_h

#include "fileseqsource.h"

class RandomSeqSource : public SeqSource
	{
public:
	bool m_Nucleo = false;
	uint m_MinSeqLength = 100;
	uint m_MaxSeqLength = 10000;
	uint m_SeqCount = 0;
	uint m_SeqIndex = 0;

public:
	virtual bool GetIsNucleo() { return m_Nucleo; };

protected:
	virtual bool GetNextLo(SeqInfo *SI);

public:
	virtual unsigned GetPctDoneX10()
		{
		if (m_SeqCount == 0)
			return 0;
		uint Pct10 = (m_SeqIndex*1000)/m_SeqCount;
		if (Pct10 >= 999)
			Pct10 = 998;
		return Pct10;
		}

	virtual const char *GetFileNameC() const
		{
		return ".random.";
		}

	virtual void Rewind()
		{
		m_SeqIndex = 0;
		}

public:
	RandomSeqSource() {}
	virtual ~RandomSeqSource() {}
	};

#endif // randomseqsource_h
