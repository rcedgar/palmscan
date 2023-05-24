#include "myutils.h"
#include "pdbchain.h"
#include "spher.h"
#include "abcxyz.h"
#include "sort.h"
#include "cavity.h"
#include "outputfiles.h"

void Palmcator(FILE *f, const PDBChain &Chain);
void RemoveHidden(const PDBChain &Chain, 
  const vector<uint> &PosVec, 
  vector<uint> &PosVec2);

void GetCavityLatLong(double x, double y, double z,
  double &r, double &lat_deg, double &long_deg)
	{
	double theta, phi;
	CartToSpher(
	  x - g_CavityCenterX,
	  y - g_CavityCenterY,
	  z - g_CavityCenterZ,
	  r, theta, phi);

	lat_deg = degrees_0_to_360(theta);
	asserta(lat_deg >= 0 && lat_deg < 180);

	long_deg = degrees_0_to_360(phi);
	asserta(long_deg >= 0 && long_deg < 360);
	}

static void SmoothDistVec(const vector<double> DistVec,
  vector<uint> &Smooth)
	{
	Smooth.clear();
	uint BinCount = uint(g_CavityMaxDist+0.5);
	Smooth.resize(BinCount+1, 0);
	const uint N = SIZE(DistVec);
	for (uint i = 0; i < N; ++i)
		{
		double r = DistVec[i];
		uint Bin = uint(r);
		asserta(Bin <= BinCount);
		++(Smooth[Bin]);
		}
	}

void FindCavity(const PDBChain &Chain, const string &PDBFileName)
	{
//	Palmcator(g_fpalmcator, Chain);

	const uint L = Chain.GetSeqLength();
	const uint PPStart = Chain.GetMotifPos(A);
	const uint PPEnd = Chain.GetMotifPos(C) + CL - 1;

	const uint MinPos = (PPStart < g_Flank_F ? 0 : PPStart - g_Flank_F);
	uint MaxPos = PPEnd + g_Flank_E;
	if (MaxPos >= L)
		MaxPos = L - 1;

	string Label;
	GetLabelFromFileName(PDBFileName, Label);
	const char *Lab = Label.c_str();

	vector<uint> PosVec;
	for (uint Pos = MinPos; Pos <= MaxPos; ++Pos)
		{
		vector<double> Pt;
		Chain.GetPt(Pos, Pt);

		double r;
		double Lat, Long;
		GetCavityLatLong(Pt[X], Pt[Y], Pt[Z], r, Lat, Long);
		if (r <= g_CavityMaxDist)
			PosVec.push_back(Pos);
		}

	if (g_fpml != 0)
		{
		bool First = true;
		fprintf(g_fpml, "cmd.select(\"cavity\", \"name CA and (");
		for (uint i = 0; i < SIZE(PosVec); ++i)
			{
			uint Pos = PosVec[i];
			uint ResNr = Chain.GetResidueNr(Pos);
			if (First)
				{
				fprintf(g_fpml, "resi %u", ResNr);
				First = false;
				}
			else
				fprintf(g_fpml, " or resi %u", ResNr);
			}
		fprintf(g_fpml, ")\")\n");
		fprintf(g_fpml, "cmd.color(\"density\", \"cavity\")\n");
		fprintf(g_fpml, "deselect\n");
		fprintf(g_fpml, "color tv_blue, %s_mA\n", Lab);
		fprintf(g_fpml, "color tv_green, %s_mB\n", Lab);
		fprintf(g_fpml, "color tv_red, %s_mC\n", Lab);
		}

	vector<uint> PosVec2;
	RemoveHidden(Chain, PosVec, PosVec2);
	if (g_fpml != 0 && !PosVec2.empty())
		{
		bool First = true;
		fprintf(g_fpml, "cmd.select(\"grip\", \"name CA and (");
		for (uint i = 0; i < SIZE(PosVec2); ++i)
			{
			uint Pos = PosVec2[i];
			uint ResNr = Chain.GetResidueNr(Pos);
			if (First)
				{
				fprintf(g_fpml, "resi %u", ResNr);
				First = false;
				}
			else
				fprintf(g_fpml, " or resi %u", ResNr);
			}
		fprintf(g_fpml, ")\")\n");
		fprintf(g_fpml, "cmd.color(\"lightorange\", \"cavity\")\n");
		fprintf(g_fpml, "deselect\n");
		fprintf(g_fpml, "color tv_blue, %s_mA\n", Lab);
		fprintf(g_fpml, "color tv_green, %s_mB\n", Lab);
		fprintf(g_fpml, "color tv_red, %s_mC\n", Lab);
		fprintf(g_fpml, "hide everything\n");
		fprintf(g_fpml, "show cartoon, grip\n");
		fprintf(g_fpml, "show cartoon, %s_mA\n", Lab);
		fprintf(g_fpml, "show cartoon, %s_mB\n", Lab);
		fprintf(g_fpml, "show cartoon, %s_mC\n", Lab);
		fprintf(g_fpml, "show spheres, %s_DGD\n", Lab);
		fprintf(g_fpml, "zoom vis\n");
		}
	}
