#ifndef alnheuristics_h
#define alnheuristics_h

struct AlnParams;

struct AlnHeuristics
	{
	unsigned BandRadius;
	unsigned HSPFinderWordLength;

	float XDropG;			//  GappedBlast default
	float XDropU;			//  UngappedBlast default
	unsigned MinGlobalHSPLength;

// For HSPs in global alignment
	float XDropGlobalHSP;
	float MinGlobalHSPScore;
	float MinGlobalHSPFractId;

	AlnHeuristics();
	void LogMe() const;
	void InitFromCmdLine(const AlnParams &AP);
	void InitGlobalFullDP();

	bool IsGlobalFull() const
		{
		return MinGlobalHSPLength == 0 && BandRadius == 0;
		}
	};

#endif // alnheuristics_h
