#include "myutils.h"

static double Gumbel(double Mu, double Beta, double x)
	{
	return 0;
	}

void cmd_gumbel()
	{
	const string &InputFN = opt_gumbel;

	FILE *f = OpenStdioFile(InputFN);

	string Line;
	vector<string> Fields;
	vector<double> Freqs;
	double SumScores = 0.0;
	double SumFreqs = 0.0;
	double Score = 0;
	double Mode = 0;
	double MaxFreq = 0;
	const double BIN_SIZE = 5.0;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		double Freq = StrToFloat(Fields[3]);
		Freqs.push_back(Freq);
		if (Freq > MaxFreq)
			{
			MaxFreq = Freq;
			Mode = Score;
			}
		SumFreqs += Freq;
		SumScores += Freq*Score;
		Score += BIN_SIZE;
		}
	CloseStdioFile(f);
	uint N = SIZE(Freqs);

	double Sum2 = 0;
	double Median = 0;
	Score = 0;
	for (uint i = 0; i < N; ++i)
		{
		Sum2 += Freqs[i];
		if (Sum2 >= SumFreqs/2)
			{
			Median = Score;
			break;
			}
		Score += BIN_SIZE;
		}
	double MeanScore = SumScores/SumFreqs;

	ProgressLog("Mean %.4g, median %.4g, mode %.4g\n",
	  MeanScore, Median, Mode);

	double Mu = double(Mode);
// Median = mu - beta ln(ln 2)
// beta = (mu - Median)/lnln2
// ln(2) = 0.6931471805599453
// ln(ln(2)) = -0.36651292058166435
	const double LNLN2 = -0.36651292058166435;
	double Beta = (Mu - Median)/LNLN2;

// gamma = 0.5772156649
// Mean = mu + beta gamma 
	const double GAMMA = 0.5772156649;
	double Beta2 = (MeanScore - Mu)/GAMMA;

	ProgressLog("Mu %.4g, beta %.4g %.4g\n", Mu, Beta, Beta2);
	}
