#include "myutils.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "abcxyz.h"
#include "outputfiles.h"

static const uint MAX_DP = 1024;

#define TRACE	0

// static XDPMem g_XDPMem;
// static float **DPScoreMx;

float ViterbiFastMem(XDPMem &Mem, float **ScoreMx,
  uint LA, uint LB, string &Path);

float **AllocDPScoreMx()
	{
	float **DPScoreMx = myalloc(float *, MAX_DP + 16);
	for (uint i = 0; i < MAX_DP + 16; ++i)
		DPScoreMx[i] = myalloc(float, MAX_DP + 16);
	return DPScoreMx;
	}

// TM-align paper eq(2).
static double Lmin = 100.0;
static double d0 = 1.24*pow(Lmin - 15.0, 0.33333) - 1.8;
static double d02 = d0*d0;

static float GetTM_Sim(double d)
	{
	double Bottom = 1.0 + d*d/d02;
	double TM_sim = 1.0/Bottom;
	return float(TM_sim);
	}

double AlignTmPpc(XDPMem &Mem, float **DPScoreMx,
  const PDBChain &Query, const PDBChain &Ref, string &Path)
	{
	const uint QLen = Query.GetSeqLength();
	const uint RLen = Ref.GetSeqLength();
	asserta(QLen > 0 && QLen < MAX_DP);
	asserta(RLen > 0 && RLen < MAX_DP);

	Mem.Alloc(MAX_DP + 16, MAX_DP + 16);
//	AllocDPScoreMx();

	vector<double> QPt(3);
	vector<double> RPt(3);
	for (uint QPos = 0; QPos < QLen; ++QPos)
		{
		Query.GetPt(QPos, QPt);
		for (uint RPos = 0; RPos < RLen; ++RPos)
			{
			Ref.GetPt(RPos, RPt);

			float d = (float) GetDist(QPt, RPt);
			float TM_Sim = GetTM_Sim(d);
			DPScoreMx[QPos][RPos] = TM_Sim;
			}
		}

#if TRACE
	{
	Log("\n");
	Log("ScoreMx\n");
	Log("      ");
	for (uint j = 0; j < RLen; ++j)
		Log("   %8u", j);
	Log("\n");

	for (uint i = 0; i < QLen; ++i)
		{
		Log("[%4u]", i);
		for (uint j = 0; j < RLen; ++j)
			{
			float d = DPScoreMx[i][j];
			Log("  %9.3g", d);
			}
		Log("\n");
		}
	Log("\n");
	}
#endif

	ViterbiFastMem(Mem, DPScoreMx, QLen, RLen, Path);
	const uint ColCount = SIZE(Path);
	uint i = 0;
	uint j = 0;
	const string &QSeq = Query.m_Seq;
	const string &RSeq = Ref.m_Seq;
	string QRow;
	string RRow;
	double SumScore = 0.0;
	uint MatchCount = 0;
	double MinScore = GetTM_Sim(12.0);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			{
			double Score = DPScoreMx[i][j];
			if (Score >= MinScore)
				{
				SumScore += Score;
				++MatchCount;
				}

			QRow += QSeq[i++];
			RRow += RSeq[j++];
			}
		else if (c == 'D')
			{
			QRow += QSeq[i++];
			RRow += '-';
			}
		else if (c == 'I')
			{
			QRow += '-';
			RRow += RSeq[j++];
			}
		}
	double TM = SumScore/MatchCount;

#if TRACE
	{
	Log("Path=%s\n", Path.c_str());
	const uint ColCount = SIZE(Path);
	uint i = 0;
	uint j = 0;
	const string &QSeq = Query.m_Seq;
	const string &RSeq = Ref.m_Seq;
	string QRow;
	string RRow;
	double SumScore = 0.0;
	uint MatchCount = 0;
	double MinScore = GetTM_Sim(4.0);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			{
			double Score = DPScoreMx[i][j];
			if (Score >= MinScore)
				{
				++MatchCount;
				SumScore += Score;
				}
			QRow += QSeq[i++];
			RRow += RSeq[j++];
			}
		else if (c == 'D')
			{
			QRow += QSeq[i++];
			RRow += '-';
			}
		else if (c == 'I')
			{
			QRow += '-';
			RRow += RSeq[j++];
			}
		}
	{
	Log("QRow %s\n", QRow.c_str());
	Log("TRow %s\n", RRow.c_str());
	asserta(i == QLen);
	asserta(j == RLen);
	double TM = SumScore/MatchCount;
	Log("TM = %.4f\n", TM);
	}
	}
#endif
	return TM;
	}

void cmd_tm_ppc()
	{
	const string &QueryFileName = opt_tm_ppc;
	const string &RefFileName = opt_ref;

	PDBChain Query;
	PDBChain Ref;
	Query.FromCal(QueryFileName);
	Ref.FromCal(RefFileName);
	Query.CheckPPCMotifCoords();
	Ref.CheckPPCMotifCoords();

	XDPMem Mem;
	float **DPScoreMx = AllocDPScoreMx();

	string Path;
	double TM = AlignTmPpc(Mem, DPScoreMx, Query, Ref, Path);
	ProgressLog("TM %6.4f  %s  %s\n",
	  TM, Query.m_Label.c_str(), Ref.m_Label.c_str());
	}

void cmd_tm_ppc_all_vs_all()
	{
	const string &InputFileName = opt_tm_ppc_all_vs_all;
	vector<PDBChain *> Chains;
	ReadChains(InputFileName, Chains);
	const uint N = SIZE(Chains);

	XDPMem Mem;
	float **DPScoreMx = AllocDPScoreMx();

	string Path;
	string QAccStr;
	string RAccStr;
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Q = *Chains[i];
		const char *QAcc = Q.GetAcc(QAccStr);
		ProgressStep(i, N, "Aligning %s", QAcc);
		for (uint j = 0; j < N; ++j)
			{
			const PDBChain &R = *Chains[j];
			const char *RAcc = R.GetAcc(RAccStr);
			double TM = AlignTmPpc(Mem, DPScoreMx, Q, R, Path);
			// Log("TM %6.4f  %s  %s\n", TM, QAcc, RAcc);
			if (g_ftsv != 0)
				{
				fprintf(g_ftsv, "%s", QAcc);
				fprintf(g_ftsv, "\t%s", RAcc);
				fprintf(g_ftsv, "\t%6.4f", TM);
				fprintf(g_ftsv, "\t%6.4f", TM);
				fprintf(g_ftsv, "\n");
				}
			}
		}
	}

void cmd_tm_search_ppc()
	{
	const string &QueryFileName = opt_tm_search_ppc;
	const string &RefFileName = opt_ref;

	vector<PDBChain *> Qs;
	vector<PDBChain *> Rs;
	ReadChains(QueryFileName, Qs);
	ReadChains(RefFileName, Rs);

	const uint NQ = SIZE(Qs);
	const uint NR = SIZE(Rs);

	XDPMem Mem;
	float **DPScoreMx = AllocDPScoreMx();

	string Path;
	string QAccStr;
	string RAccStr;
	for (uint i = 0; i < NQ; ++i)
		{
		const PDBChain &Q = *Qs[i];
		const char *QAcc = Q.GetAcc(QAccStr);
		ProgressStep(i, NQ, "Aligning %s", QAcc);
		string TopRAccStr;
		double TopTM = 0;
		for (uint j = 0; j < NR; ++j)
			{
			const PDBChain &R = *Rs[j];
			const char *RAcc = R.GetAcc(RAccStr);
			double TM = AlignTmPpc(Mem, DPScoreMx, Q, R, Path);
			if (opt_top_hit_only)
				{
				if (TM > TopTM)
					{
					TopTM = TM;
					TopRAccStr = RAccStr;
					}
				}
			else
				{
				if (g_ftsv != 0)
					{
					fprintf(g_ftsv, "%s", QAcc);
					fprintf(g_ftsv, "\t%s", RAcc);
					fprintf(g_ftsv, "\t%6.4f", TM);
					fprintf(g_ftsv, "\n");
					}
				}
			}

		if (opt_top_hit_only)
			{
			if (g_ftsv != 0)
				{
				fprintf(g_ftsv, "%s", QAcc);
				fprintf(g_ftsv, "\t%s", TopRAccStr.c_str());
				fprintf(g_ftsv, "\t%6.4f", TopTM);
				fprintf(g_ftsv, "\n");
				}
			}
		}
	}
