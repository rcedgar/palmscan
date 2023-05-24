#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "outputfiles.h"

static void StaticLogVec(const string &Name, const vector<double> &v)
	{
	asserta(SIZE(v) == 3);
	Log("  x=%8.3f,  y=%8.3f,  z=%8.3f  %s\n", v[X], v[Y], v[Z], Name.c_str());
	}

static void HSETest(const PDBChain &Chain, vector<double> &Angles)
	{
//	Angles.clear();
	const uint L = Chain.GetSeqLength();
	for (uint Pos = 1; Pos + 1 < L; ++Pos)
		{
		char c = Chain.m_Seq[Pos];

		vector<double> PtPrevCA;
		vector<double> PtNextCA;
		vector<double> PtCA;
		vector<double> PtCB;
		Chain.GetPt(Pos-1, PtPrevCA);
		Chain.GetPt(Pos, PtCA);
		Chain.GetPt(Pos+1, PtNextCA);

		bool Ok = Chain.Get_CB_Pt(Pos, PtCB);
		if (!Ok)
			continue;

		vector<double> VecAB;
		Sub_Vecs(PtCB, PtCA, VecAB);
		NormalizeVec(VecAB);

		vector<double> d1;
		vector<double> d2;
		Sub_Vecs(PtCA, PtPrevCA, d1);
		Sub_Vecs(PtCA, PtNextCA, d2);

		vector<double> VecPAB;
		Add_Vecs(d1, d2, VecPAB);

		NormalizeVec(VecPAB);

		double Theta = GetTheta_Vecs(VecPAB, VecAB);
		double Deg = degrees(Theta);
		Angles.push_back(Deg);

#if TRACE
		Log("Pos [%4u] %c\n", Pos, c);
		StaticLogVec("CA", PtCA);
		StaticLogVec("CB", PtCB);
		StaticLogVec("PrevCA", PtPrevCA);
		StaticLogVec("NextCA", PtNextCA);
		StaticLogVec("norm(d1+d2)", VecPAB);
		StaticLogVec("norm(CA->CB)", VecAB);
		Log(" theta = %.1f\n", degrees(Theta));
		Log("\n");
		Log("\n");
#endif
		}
	}

void cmd_pseudo_cb_angles()
	{
	const string &FN = opt_pseudo_cb_angles;

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);
	uint FlankSize = 0;

	const uint ChainCount = SIZE(Chains);
	vector<double> Angles;
	for (uint i = 0; i < ChainCount; ++i)
		{
		ProgressStep(i, ChainCount, "Processing");
		PDBChain &Chain = *Chains[i];
		HSETest(Chain, Angles);
		}

	uint W = 5;
	uint BINS = 360/W;
	vector<uint> BinCounts(BINS+1, 0);
	uint MaxBin = 0;
	for (uint i = 0; i < SIZE(Angles); ++i)
		{
		double Angle = Angles[i];
		uint Bin = uint(Angle/W + 0.5);
		asserta(Bin <= BINS);
		BinCounts[Bin] += 1;
		MaxBin = max(Bin, MaxBin);
		}

	for (uint Bin = 0; Bin <= MaxBin; ++Bin)
		{
		uint n = BinCounts[Bin];
		double Center = double(Bin*W) + W/2.0;
		Log("%6.1f  %u\n", Center, n);
		if (g_ftsv != 0)
			fprintf(g_ftsv, "%.1f\t%u\n", Center, n);
		}
	}
