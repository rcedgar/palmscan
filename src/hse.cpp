#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include <map>

void StaticLogVec(const string &Name, const vector<double> &v)
	{
	asserta(SIZE(v) == 3);
	Log("  x=%8.3f,  y=%8.3f,  z=%8.3f  %s\n", v[X], v[Y], v[Z], Name.c_str());
	}

void GetHSE(const PDBChain &Chain, uint Pos, double Radius,
  uint MinPos, uint MaxPos, uint &NU, uint &ND)
	{
	NU = 0;
	ND = 0;
	if (Pos == 0 || Pos+1 >= SIZE(Chain.m_Seq))
		return;

	vector<double> PtPrevCA;
	vector<double> PtNextCA;
	vector<double> PtCA;
	vector<double> PtCB;
	Chain.GetPt(Pos-1, PtPrevCA);
	Chain.GetPt(Pos, PtCA);
	Chain.GetPt(Pos+1, PtNextCA);

	vector<double> d1;
	vector<double> d2;
	Sub_Vecs(PtCA, PtPrevCA, d1);
	Sub_Vecs(PtCA, PtNextCA, d2);

	vector<double> VecPAB;
	Add_Vecs(d1, d2, VecPAB);
	NormalizeVec(VecPAB);

	vector<uint> SpherePosVec;
	Chain.GetSphere(Pos, Radius, MinPos, MaxPos, SpherePosVec);

	const uint N = SIZE(SpherePosVec);
	vector<double> Pt2;
	vector<double> Vec12;
	for (uint i = 0; i < N; ++i)
		{
		uint Pos2 = SpherePosVec[i];
		Chain.GetPt(Pos2, Pt2);
		Sub_Vecs(Pt2, PtCA, Vec12);
		double Theta = GetTheta_Vecs(VecPAB, Vec12);
		double Deg = degrees(Theta);
		if (Deg < 90)
			++NU;
		else
			++ND;
		}
	}

static void HSE(const PDBChain &Chain, double Radius)
	{
	const uint L = Chain.GetSeqLength();
	for (uint Pos = 1; Pos + 1 < L; ++Pos)
		{
		char c = Chain.m_Seq[Pos];
		uint L = Chain.GetSeqLength();

		uint NU, ND;
		GetHSE(Chain, Pos, Radius, 0, L-1, NU, ND);
		if (g_ftsv != 0)
			{
			fprintf(g_ftsv, "%s", Chain.m_Label.c_str());
			fprintf(g_ftsv, "\t%u", Pos+1);
			fprintf(g_ftsv, "\t%c", c);
			fprintf(g_ftsv, "\t%u", NU);
			fprintf(g_ftsv, "\t%u", ND);
			fprintf(g_ftsv, "\n");
			}
		}
	}

void cmd_hse()
	{
	const string &FN = opt_hse;

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);
	uint FlankSize = 0;

	const uint ChainCount = SIZE(Chains);
	vector<double> Angles;
	for (uint i = 0; i < ChainCount; ++i)
		{
		ProgressStep(i, ChainCount, "Processing");
		PDBChain &Chain = *Chains[i];
		HSE(Chain, 12.0);
		}
	}
