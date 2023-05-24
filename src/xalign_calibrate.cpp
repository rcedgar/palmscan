#include "myutils.h"
#include "mx.h"
#include "xprof.h"
#include "xprofdata.h"
#include "xtrainer.h"
#include "outputfiles.h"
#include "omplock.h"

static uint g_QueryCount;
static uint g_HitCount;

double XAlign(Mx<float> &SMx, 
  Mx<float> &a_FwdM,
  Mx<float> &a_FwdD,
  Mx<float> &a_FwdI,
  Mx<char> &a_TBM,
  Mx<char> &a_TBD,
  Mx<char> &a_TBI,
  const XProfData &Prof1, const XProfData &Prof2,
  const vector<vector<vector<double> > > &FeatureIndexToLogOddsMx,
  double MinScore);

static const double BIN_SIZE = 5.0;
static const uint MAX_BIN = 50;
static const double MAX_SCORE = MAX_BIN*BIN_SIZE;

static uint ScoreToBin(double Score)
	{
	double r = Score/MAX_SCORE;
	if (r > 1.0)
		r = 1.0;
	uint Bin = uint(r*MAX_BIN);
	asserta(Bin <= MAX_BIN);
	return Bin;
	}

static void Calibrate(const string &Label, 
  const vector<double> &Scores)
	{
	FILE *f = g_fcal;
	if (f == 0)
		return;

	const uint N = SIZE(Scores);
	vector<uint> Counts(MAX_BIN+1);
	for (uint i = 0; i < N; ++i)
		{
		double Score = Scores[i];
		uint Bin = ScoreToBin(Score);
		if (Bin > MAX_BIN)
			Bin = MAX_BIN;
		Counts[Bin] += 1;
		}

	for (uint Bin = 0; Bin <= MAX_BIN; ++Bin)
		{
		uint n = Counts[Bin];
		double Freq = double(n)/N;
		fprintf(f, "%s\t%u\t%u\t%.4g\n",
		  Label.c_str(), Bin, n, Freq);
		}
	}

static void Thread(FILE *fQ, const vector<XProfData *> &DBXDs,
  const vector<vector<vector<double> > > &FeatureIndexToLogOddsMx)
	{
	const uint DBN = SIZE(DBXDs);
	Mx<float> SMx;
	Mx<float> a_FwdM;
	Mx<float> a_FwdD;
	Mx<float> a_FwdI;
	Mx<char> a_TBM;
	Mx<char> a_TBD;
	Mx<char> a_TBI;

	XProfData QP;
	vector<double> Scores;
	vector<double> Scores2;
	for (;;)
		{
		bool Ok;
#pragma omp critical
		{
		Ok = QP.FromCfv(fQ);
		}

		if (!Ok)
			return;

#pragma omp critical
		{
		++g_QueryCount;
		if (g_QueryCount%1 == 0)
			Progress("%u queries, %u hits\r", g_QueryCount, g_HitCount);
		}

		Scores.clear();
		Scores2.clear();
		Scores.reserve(DBN);
		Scores.reserve(DBN);
		for (uint i = 0; i < DBN; ++i)
			{
			const XProfData &DP = *DBXDs[i];
			double Score = XAlign(SMx, a_FwdM, a_FwdD, a_FwdI,
			  a_TBM, a_TBD, a_TBI,
				QP, DP, FeatureIndexToLogOddsMx, 0);
			Scores.push_back(Score);
			}

#if 0
		QP.Shuffle();
		for (uint i = 0; i < DBN; ++i)
			{
			const XProfData &DP = *DBXDs[i];
			double Score = XAlign(SMx, a_FwdM, a_FwdD, a_FwdI,
			  a_TBM, a_TBD, a_TBI,
				QP, DP, FeatureIndexToLogOddsMx, 0);
			Scores2.push_back(Score);
			}
#endif
		Calibrate(QP.m_Label, Scores);
		}
	}

void cmd_xalign_calibrate()
	{
	const string &QueryFileName = opt_xalign_calibrate;
	const string &DBFileName = opt_db;

	vector<XProfData *> DBXDs;
	Progress("Reading db... ");
	FILE *fDB = OpenStdioFile(DBFileName);
	for (;;)
		{
		XProfData *XD = new XProfData;
		bool Ok = XD->FromCfv(fDB);
		if (!Ok)
			break;
		DBXDs.push_back(XD);
		}
	CloseStdioFile(fDB);
	fDB = 0;
	const uint DBN = SIZE(DBXDs);
	Progress(" %u profiles\n", DBN);

	vector<vector<vector<double> > > FeatureIndexToLogOddsMx;
	XTrainer::LogOddsFromTsv(opt_logodds, FeatureIndexToLogOddsMx);

	XProfData QP;
	FILE *fQ = OpenStdioFile(QueryFileName);

	uint ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount)
	Thread(fQ, DBXDs, FeatureIndexToLogOddsMx);

	ProgressLog("%u queries, %u hits\n", g_QueryCount, g_HitCount);
	}
