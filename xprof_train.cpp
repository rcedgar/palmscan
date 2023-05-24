#include "myutils.h"
#include "pdbchain.h"
#include "xtrainer.h"
#include "xprof.h"
#include "outputfiles.h"

float SW(const Mx<float> &SMx, 
  Mx<float> &a_FwdM,
  Mx<float> &a_FwdD,
  Mx<float> &a_FwdI,
  Mx<char> &a_TBM,
  Mx<char> &a_TBD,
  Mx<char> &a_TBI,
  uint &Starti, uint &Startj, string &Path);

static char ScoreChar(double Score, double MaxScore)
	{
	asserta(Score <= MaxScore);
	if (Score <= 0)
		return ' ';

	if (Score >= MaxScore*0.75)
		return '*';
	if (Score >= MaxScore*0.5)
		return '+';
	if (Score >= MaxScore*0.25)
		return '.';
	if (Score >= MaxScore*0.1)
		return '_';
	return ' ';
	}

void cmd_xprof_train()
	{
	vector<PDBChain *> Chains;
	ReadChains(opt_xprof_train, Chains);

	XTrainer XT;
	XT.SetChains(Chains);

	const double MINTM = 0.6;
	const double MAXTM = 0.9;

	FILE *f = OpenStdioFile(opt_ref);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 7);
		double TM = StrToFloat(Fields[0]);
		if (TM < MINTM || TM > MAXTM)
			continue;
		const string &Label1 = Fields[1];
		const string &Label2 = Fields[2];
		const string &Row1 = Fields[5];
		const string &Row2 = Fields[6];

		XT.SetAln(Label1, Label2, Row1, Row2);
		}

	const uint FeatureCount = XProf::GetFeatureCount();
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		ProgressStep(FeatureIndex, FeatureCount, "Feature value dist");
		XT.SetValues(FeatureIndex);

		Log("\n");
		Log("Feature %s\n", XProf::GetFeatureName(FeatureIndex));
		XT.LogValueHist();
		}

	XT.Train();
	XTrainer::LogOddsToTsv(opt_logodds,
	  XT.m_FeatureIndexToLogOddsMx);

	Log("\n\n");
	for (uint i = 0; i < FeatureCount; ++i)
		{
		ProgressStep(i, FeatureCount, "Correl");
		const char *Namei = XProf::GetFeatureName(i);
		for (uint j = i+1; j < FeatureCount; ++j)
			{
			const char *Namej = XProf::GetFeatureName(j);
			double r = XT.CorrelFeatures(i, j);
			Log("Correl  %10.10s  %10.10s  %6.3f\n", Namei, Namej, r);
			}
		}

	XT.LogFinal();

	const uint L0 = XT.m_Chains[0]->GetSeqLength();
	const uint L1 = XT.m_Chains[1]->GetSeqLength();

	vector<vector<double> > StdMx;
	XT.GetScoreMx(0, 1, StdMx);

	if (g_ftsv_dot != 0)
		{
		FILE *f = g_ftsv_dot;

		double MaxScore = 0;
		for (uint i = 0; i < L0; ++i)
			for (uint j = 0; j < L1; ++j)
				MaxScore = max(MaxScore, StdMx[i][j]);

		fprintf(f, "%u\t%u\n", L0, L1);
		for (uint i = 0; i < L0; ++i)
			{
			for (uint j = 0; j < L1; ++j)
				{
				double Score = StdMx[i][j];
				fprintf(f, "%.3f\n", Score/MaxScore);
				}
			}
		}

	Mx<float> SMx;
	XT.GetScoreMx2(0, 1, SMx);

	uint Start0, Start1;
	string Path;
	Mx<float> a_FwdM;
	Mx<float> a_FwdD;
	Mx<float> a_FwdI;
	Mx<char> a_TBM;
	Mx<char> a_TBD;
	Mx<char> a_TBI;
	float Score = SW(SMx, a_FwdM, a_FwdD, a_FwdI,
	  a_TBM, a_TBD, a_TBI, Start0, Start1, Path);
	uint Cols = SIZE(Path);
	Log("Score = %.1f, Cols=%u, Path = %s\n", Score, Cols, Path.c_str());

	const string &A = Chains[0]->m_Seq;
	const string &B = Chains[1]->m_Seq;
	const uint ColCount = SIZE(Path);
	uint i = Start0;
	uint j = Start1;
	string ARow;
	string BRow;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'D')
			{
			ARow += A[i];
			++i;
			}
		else
			ARow += '-';

		if (c == 'M' || c == 'I')
			{
			BRow += B[j];
			++j;
			}
		else
			BRow += '-';
		}
	Log("A  >%s\n", Chains[0]->m_Label.c_str());
	Log("B  >%s\n", Chains[1]->m_Label.c_str());
	Log("A  %s\n", ARow.c_str());
	Log("B  %s\n", BRow.c_str());
	}
