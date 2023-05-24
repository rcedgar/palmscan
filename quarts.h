#ifndef quarts_h
#define quarts_h

struct Quarts
	{
	unsigned Min;
	unsigned LoQ;
	unsigned Med;
	unsigned HiQ;
	unsigned Max;
	unsigned Total;
	double Avg;

	Quarts()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		}

	void Clear()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		}

	void LogMe() const
		{
		Log("Min=%u", Min);
		Log(", LoQ=%u", LoQ);
		Log(", Med=%u", Med);
		Log(", HiQ=%u", HiQ);
		Log(", Max=%u", Max);
		Log(", Avg=%.1f", Avg);
		Log("\n");
		}
	};

struct QuartsDouble
	{
	double Min;
	double LoQ;
	double Med;
	double HiQ;
	double Max;
	double Total;
	double Avg;
	double StdDev;

	QuartsDouble()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		StdDev = 0;
		}

	void Clear()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		StdDev = 0;
		}

	void LogMe() const
		{
		Log("Min=%.3g", Min);
		Log(", LoQ=%.3g", LoQ);
		Log(", Med=%.3g", Med);
		Log(", HiQ=%.3g", HiQ);
		Log(", Max=%.3g", Max);
		Log(", Avg=%.3g", Avg);
		Log(", StdDev=%.3g", StdDev);
		Log("\n");
		}
	};

void GetQuarts(const vector<unsigned> &v, Quarts &Q);
void GetQuartsDouble(const vector<double> &v, QuartsDouble &Q);

#endif // quarts_h
