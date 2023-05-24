#include "myutils.h"
#include "structprof.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "cmpsearcher.h"
#include "gsprof.h"
#include "abcxyz.h"

static uint g_APos;
static uint g_BPos;
static uint g_CPos;
static uint g_DPos;
static uint g_EPos;
static uint g_F1Pos;
static uint g_F2Pos;

void StructProf::SetChain(const PDBChain &Chain)
	{
	Clear();
	m_Chain = &Chain;
	uint L = Chain.GetSeqLength();
	asserta(L > 0);
	}

void StructProf::SetCavityCenterPt()
	{
	vector<vector<double> > MotifCoords;
	vector<vector<double> > Basis;
	vector<double> CentroidCoords;
	m_Chain->GetMotifCoords(MotifCoords);
	GetTriangleBasis(MotifCoords, CentroidCoords, Basis);

	m_CavityCenterPt.resize(3);
	m_CavityCenterPt[X] = CentroidCoords[X];
	m_CavityCenterPt[Y] = CentroidCoords[Y];
	m_CavityCenterPt[Z] = CentroidCoords[Z];

	m_CavityCenterPt[X] += 3.0*Basis[X][X];
	m_CavityCenterPt[Y] += 3.0*Basis[X][Y];
	m_CavityCenterPt[Z] += 3.0*Basis[X][Z];

	m_CavityCenterPt[X] += -10.0*Basis[Z][X];
	m_CavityCenterPt[Y] += -10.0*Basis[Z][Y];
	m_CavityCenterPt[Z] += -10.0*Basis[Z][Z];
	}

uint StructProf::SearchDist(uint Pos, uint Lo, uint Hi,
  bool Maximize, double X, double &BestDist, bool Trace) const
	{
	if (Trace)
		Log("SearchDist(Pos=%u, Lo=%u, Hi=%u, %s, X=%.1f)\n",
		  Pos, Lo, Hi, Maximize ? "MAX" : "min", X);

	const PDBChain &Chain = *m_Chain;
	const uint L = Chain.GetSeqLength();
//	asserta(Lo < Hi && Hi < L);
	if (Lo >= Hi || Hi >= L)
		return UINT_MAX;
	BestDist = (Maximize ? 0 : DBL_MAX);
	uint BestPos = UINT_MAX;
	for (uint Pos2 = Lo; Pos2 <= Hi; ++Pos2)
		{
		double Dist = Chain.GetDist(Pos2, Pos);
		if (Trace)
			Log("  Pos2 [%4u]  %5.1f", Pos2, Dist);
		if (Maximize)
			{
			if (Dist > BestDist)
				{
				BestDist = Dist;
				BestPos = Pos2;
				if (Trace)
					Log(" << best");
				}
			else if (Dist < BestDist - X)
				{
				if (Trace)
					Log(" << X-drop\n");
				break;
				}
			if (Trace)
				Log("\n");
			}
		else
			{
			if (Dist < BestDist)
				{
				BestDist = Dist;
				BestPos = Pos2;
				if (Trace)
					Log(" << best");
				}
			else if (Dist > BestDist + X)
				{
				if (Trace)
					Log(" << X-drop\n");
				break;
				}
			if (Trace)
				Log("\n");
			}
		}
	return BestPos;
	}

void StructProf::GetSphere(const vector<double> &CenterPt,
  double Radius, vector<uint> &PosVec) const
	{
	PosVec.clear();
	const uint L = m_Chain->GetSeqLength();
	//asserta(m_MinPos <= m_MaxPos && m_MaxPos < L);
	vector<double> Pt;
	//for (uint Pos2 = m_MinPos; Pos2 <= m_MaxPos; ++Pos2)
	for (uint Pos2 = 0; Pos2 < L; ++Pos2)
		{
		m_Chain->GetPt(Pos2, Pt);
		double d = GetDist(CenterPt, Pt);
		if (d <= Radius)
			PosVec.push_back(Pos2);
		}
	}

void StructProf::GetHSE(uint Pos, double Radius,
  uint &NU, uint &ND) const
	{
	NU = 0;
	ND = 0;
	const uint L = m_Chain->GetSeqLength();
	if (Pos == 0 || Pos+1 >= L)
		return;

	vector<double> PtPrevCA;
	vector<double> PtNextCA;
	vector<double> PtCA;
	vector<double> PtCB;
	m_Chain->GetPt(Pos-1, PtPrevCA);
	m_Chain->GetPt(Pos, PtCA);
	m_Chain->GetPt(Pos+1, PtNextCA);

	vector<double> d1;
	vector<double> d2;
	Sub_Vecs(PtCA, PtPrevCA, d1);
	Sub_Vecs(PtCA, PtNextCA, d2);

	vector<double> VecPAB;
	Add_Vecs(d1, d2, VecPAB);
	NormalizeVec(VecPAB);

	vector<uint> SpherePosVec;
	GetSphere(PtCA, Radius, SpherePosVec);

	const uint N = SIZE(SpherePosVec);
	vector<double> Pt2;
	vector<double> Vec12;
	for (uint i = 0; i < N; ++i)
		{
		uint Pos2 = SpherePosVec[i];
		if (Pos2 == Pos)
			continue;
		m_Chain->GetPt(Pos2, Pt2);
		Sub_Vecs(Pt2, PtCA, Vec12);
		double Theta = GetTheta_Vecs(VecPAB, Vec12);
		double Deg = degrees(Theta);
		if (Deg < 90)
			++NU;
		else
			++ND;
		}
	}

uint StructProf::GetCavityNumber(uint Pos) const
	{
	vector<double> Pt;
	m_Chain->GetPt(Pos, Pt);

	vector<double> CenterPt;
	Add_Vecs(Pt, m_CavityCenterPt, CenterPt);

	for (uint i = 0; i < 3; ++i)
		CenterPt[i] /= 2.0;

	double Dist = GetDist3D(Pt[X], Pt[Y], Pt[Z],
	  CenterPt[X], CenterPt[Y], CenterPt[Z]);

	vector<uint> SpherePosVec;
	GetSphere(CenterPt, Dist/2.0, SpherePosVec);
	uint CN = SIZE(SpherePosVec);
	return CN;
	}

uint StructProf::GetTSB(uint Pos, double Radius) const
	{
	const uint L = m_Chain->GetSeqLength();
	if (Pos == 0 || Pos+1 >= L)
		return 0;

	vector<double> PtPrevCA;
	vector<double> PtNextCA;
	vector<double> PtCA;
	vector<double> PtCB;
	m_Chain->GetPt(Pos-1, PtPrevCA);
	m_Chain->GetPt(Pos, PtCA);
	m_Chain->GetPt(Pos+1, PtNextCA);

	vector<double> d1;
	vector<double> d2;
	Sub_Vecs(PtCA, PtPrevCA, d1);
	Sub_Vecs(PtCA, PtNextCA, d2);

	vector<double> VecPAB;
	Add_Vecs(d1, d2, VecPAB);
	NormalizeVec(VecPAB);

	vector<double> CenterPt;
	for (uint i = 0; i < 3; ++i)
		{
		double Coord = PtCA[i] + Radius*VecPAB[i];
		CenterPt.push_back(Coord);
		}

	vector<uint> SpherePosVec;
	GetSphere(CenterPt, Radius, SpherePosVec);

	const uint N = SIZE(SpherePosVec);
	return N;
	}

static void DoStructProfPos(FILE *f, const StructProf &SP, uint Pos)
	{
	if (f == 0)
		return;

	const PDBChain &Chain = *SP.m_Chain;
	const char *Label = Chain.m_Label.c_str();
	const char *Seq = Chain.m_Seq.c_str();
	char aa = Seq[Pos];
	char ss = Chain.m_SS[Pos];
	uint NU, ND;
	SP.GetHSE(Pos, 12.0, NU, ND);
	uint TSB = SP.GetTSB(Pos, 10.0);
	uint CN = SP.GetCavityNumber(Pos);

	uint Pos_aD = Chain.GetMotifPos(A) + 3;
	uint Pos_bG = Chain.GetMotifPos(B) + 1;
	uint Pos_cD = Chain.GetMotifPos(C) + 3;
	asserta(Seq[Pos_aD] == 'D');
	asserta(Seq[Pos_bG] == 'G');
	asserta(Seq[Pos_cD] == 'D');

	double Dist_aD = Chain.GetDist(Pos, Pos_aD);
	double Dist_bG = Chain.GetDist(Pos, Pos_bG);
	double Dist_cD = Chain.GetDist(Pos, Pos_cD);

	static bool HdrDone = false;
	if (!HdrDone)
		{
		fprintf(f, "Label");
		fprintf(f, "\tPos");
		fprintf(f, "\taa");
		fprintf(f, "\tss");
		fprintf(f, "\tNU");
		fprintf(f, "\tND");
		fprintf(f, "\tTSB");
		fprintf(f, "\tDist_aD");
		fprintf(f, "\tDist_bG");
		fprintf(f, "\tDist_cD");
		fprintf(f, "\tCN");
		fprintf(f, "\tMotif");
		fprintf(f, "\n");
		HdrDone = true;
		}

	int ResNr = SP.m_Chain->GetResidueNr(Pos);
	char Motif = '.';
	if (Pos == g_APos)
		Motif = 'A';
	if (Pos == g_BPos)
		Motif = 'B';
	if (Pos == g_CPos)
		Motif = 'C';
	if (Pos == g_DPos)
		Motif = 'D';
	if (Pos == g_EPos)
		Motif = 'E';
	if (Pos == g_F2Pos)
		Motif = 'F';

	fprintf(f, "%s", Label);
	fprintf(f, "\t%u", Pos+1);
	fprintf(f, "\t%c", aa);
	fprintf(f, "\t%c", ss);
	fprintf(f, "\t%u", NU);
	fprintf(f, "\t%u", ND);
	fprintf(f, "\t%u", TSB);
	fprintf(f, "\t%.1f", Dist_aD);
	fprintf(f, "\t%.1f", Dist_bG);
	fprintf(f, "\t%.1f", Dist_cD);
	fprintf(f, "\t%u", CN);
	fprintf(f, "\t%c", Motif);
	if (ResNr == INT_MAX)
		fprintf(f, "\t.");
	else
		fprintf(f, "\t%d", ResNr);
	fprintf(f, "\n");
	}

static bool DoStructProf(FILE *f, CMPSearcher &CS,
 PDBChain &Chain)
	{
	//CS.Search(Chain);

	//g_APos = UINT_MAX;
	//g_BPos = UINT_MAX;
	//g_CPos = UINT_MAX;
	//double PalmScore = CS.GetPSSMStarts(g_APos, g_BPos, g_CPos);
	//if (PalmScore <= 0)
	//	return false;

	//Chain.SetMotifPosVec(g_APos, g_BPos, g_CPos);
	//Chain.PrintSeqCoords(g_fLog);

	const uint L = Chain.GetSeqLength();
	uint MinPos = (g_APos < 150 ? 0 : g_APos - 150);
	uint MaxPos = g_CPos + CL + 150 - 1;
	if (MaxPos >= L)
		MaxPos = L - 1;

	StructProf SP;
	SP.SetChain(Chain);
	//SP.SetMinMaxPos(MinPos, MaxPos);
	SP.SetCavityCenterPt();
	SP.WriteGSProf(g_fgsprof_out);
	SP.WriteTsv(g_fmotifs);

	g_DPos = SP.FindMofifD_Hueuristics();
	g_EPos = UINT_MAX;
	if (g_DPos != UINT_MAX)
		g_EPos = SP.FindMofifE_Hueuristics(g_DPos);
	g_F2Pos = UINT_MAX;
	if (g_APos > 25)
		g_F2Pos = SP.FindMofifF2_Hueuristics(g_APos);

	for (uint Pos = MinPos; Pos <= MaxPos; ++Pos)
		DoStructProfPos(f, SP, Pos);

	if (g_fpml != 0)
		{
		string Label;
		GetLabelFromFileName(opt_struct_prof, Label);
		fprintf(g_fpml, "window hide\n");
		fprintf(g_fpml, "cmd.load(\"%s\")\n", opt_struct_prof);
		fprintf(g_fpml, "select %s\n", Label.c_str());
		fprintf(g_fpml, "color gray40, %s\n", Label.c_str());
		if (optset_png)
			fprintf(g_fpml, "bg_color white\n");

		if (g_APos != UINT_MAX)
			{
			uint Res1 = Chain.GetResidueNr(g_APos);
			uint Res2 = Chain.GetResidueNr(g_APos + AL - 1);
			fprintf(g_fpml, "cmd.select(\"motifA\", \"%s and resi %u-%u\")\n",
			  Label.c_str(), Res1, Res2);
			fprintf(g_fpml, "color tv_blue, motifA\n");

			fprintf(g_fpml, "cmd.select(\"motifA_D\", \"%s and resi %d and name CA\")\n",
			  Label.c_str(), Res1+3);
			fprintf(g_fpml, "show sphere, motifA_D\n");
			fprintf(g_fpml, "set sphere_transparency, 0.3\n");
			}
		if (g_BPos != UINT_MAX)
			{
			int Res1 = Chain.GetResidueNr(g_BPos);
			int Res2 = Chain.GetResidueNr(g_BPos + BL - 1);
			fprintf(g_fpml, "cmd.select(\"motifB\", \"%s and resi %d-%d\")\n",
			  Label.c_str(), Res1, Res2);
			fprintf(g_fpml, "color tv_green, motifB\n");

			fprintf(g_fpml, "cmd.select(\"motifB_G\", \"%s and resi %d and name CA\")\n",
			  Label.c_str(), Res1+1);
			fprintf(g_fpml, "show sphere, motifB_G\n");
			fprintf(g_fpml, "set sphere_transparency, 0.3\n");
			}
		if (g_CPos != UINT_MAX)
			{
			int Res1 = Chain.GetResidueNr(g_CPos);
			int Res2 = Chain.GetResidueNr(g_CPos + CL - 1);
			fprintf(g_fpml, "cmd.select(\"motifC\", \"%s and resi %d-%d\")\n",
			  Label.c_str(), Res1, Res2);
			fprintf(g_fpml, "color tv_red, motifC\n");

			fprintf(g_fpml, "cmd.select(\"motifC_D\", \"%s and resi %d and name CA\")\n",
			  Label.c_str(), Res1+3);
			fprintf(g_fpml, "show sphere, motifC_D\n");
			fprintf(g_fpml, "set sphere_transparency, 0.3\n");
			}
		if (g_DPos != UINT_MAX)
			{
			int Res1 = Chain.GetResidueNr(g_DPos);
			int Res2 = Chain.GetResidueNr(g_DPos + 6);
			fprintf(g_fpml, "cmd.select(\"motifD\", \"%s and resi %u-%u\")\n",
			  Label.c_str(), Res1, Res2);
			fprintf(g_fpml, "color yellow, motifD\n");

			fprintf(g_fpml, "cmd.select(\"motifD_X\", \"%s and resi %d and name CA\")\n",
			  Label.c_str(), Res1+3);
			fprintf(g_fpml, "show sphere, motifD_X\n");
			fprintf(g_fpml, "set sphere_transparency, 0.3\n");
			}
		if (g_EPos != UINT_MAX)
			{
			int Res1 = Chain.GetResidueNr(g_EPos);
			int Res2 = Chain.GetResidueNr(g_EPos + 6);
			fprintf(g_fpml, "cmd.select(\"motifE\", \"%s and resi %u-%u\")\n",
			  Label.c_str(), Res1, Res2);
			fprintf(g_fpml, "color orange, motifE\n");

			fprintf(g_fpml, "cmd.select(\"motifE_X\", \"%s and resi %d and name CA\")\n",
			  Label.c_str(), Res1+3);
			fprintf(g_fpml, "show sphere, motifE_X\n");
			fprintf(g_fpml, "set sphere_transparency, 0.3\n");
			}
		if (g_F2Pos != UINT_MAX)
			{
			int Res1 = Chain.GetResidueNr(g_F2Pos);
			int Res2 = Chain.GetResidueNr(g_F2Pos + 6);
			fprintf(g_fpml, "cmd.select(\"motifF\", \"%s and resi %u-%u\")\n",
			  Label.c_str(), Res1, Res2);
			fprintf(g_fpml, "color cyan, motifF\n");

			fprintf(g_fpml, "cmd.select(\"motifF_R\", \"%s and resi %d and name CA\")\n",
			  Label.c_str(), Res1+2);
			fprintf(g_fpml, "show sphere, motifF_R\n");
			fprintf(g_fpml, "set sphere_transparency, 0.3\n");
			}
		fprintf(g_fpml, "deselect\n");
		if (optset_png)
			{
			fprintf(g_fpml, "ray\n");
			fprintf(g_fpml, "png %s\n", opt_png);
			fprintf(g_fpml, "quit\n");
			}
		}

	return true;
	}

void StructProf::WriteMotifTsv(FILE *f, uint Pos, uint n) const
	{
	if (f == 0)
		return;
	if (Pos == UINT_MAX)
		{
		fprintf(f, "\t.\t.");
		return;
		}

	const uint L = m_Chain->GetSeqLength();
	asserta(Pos + n <= L);
	asserta(m_Chain != 0);
	const PDBChain &Chain = *m_Chain;
	const string &Seq = Chain.m_Seq;

	fprintf(f, "\t%u\t", Pos+1);
	for (uint i = 0; i < n; ++i)
		fputc(Seq[Pos+i], f);
	return;
	}

void StructProf::WriteTsv(FILE *f) const
	{
	if (f == 0)
		return;
	const uint L = m_Chain->GetSeqLength();
	asserta(m_Chain != 0);
	const PDBChain &Chain = *m_Chain;

	uint PosA = Chain.GetMotifPos(A);
	uint PosB = Chain.GetMotifPos(B);
	uint PosC = Chain.GetMotifPos(C);
	uint PosD = FindMofifD_Hueuristics();
	uint PosE = FindMofifE_Hueuristics(PosD);
	uint PosF = FindMofifF2_Hueuristics(PosA);

	if (PosC + CL + 2 > L)
		PosC = UINT_MAX;
	if (PosD + 7 > L)
		PosD = UINT_MAX;
	if (PosE + 7 > L)
		PosE = UINT_MAX;

	fprintf(f, "%s", Chain.m_Label.c_str());
	WriteMotifTsv(f, PosA, AL);
	WriteMotifTsv(f, PosB, BL);
	WriteMotifTsv(f, PosC, CL+2);
	WriteMotifTsv(f, PosD, 7);
	WriteMotifTsv(f, PosE, 7);
	WriteMotifTsv(f, PosF, 7);
	fprintf(f, "\n");
	}

void StructProf::WriteGSProf(FILE *f) const
	{
	if (f == 0)
		return;

	const uint L = m_Chain->GetSeqLength();
	asserta(m_Chain != 0);
	const PDBChain &Chain = *m_Chain;

	uint Pos_aD = Chain.GetMotifPos(A) + 3;
	uint Pos_bG = Chain.GetMotifPos(B) + 1;
	uint Pos_cD = Chain.GetMotifPos(C) + 3;

	string Annot(L, '.');
	Annot[Pos_aD] = 'A';
	Annot[Pos_bG] = 'C';
	Annot[Pos_cD] = 'D';

	GSProf GP;
	GP.Init(Chain.m_Label, Chain.m_Seq, Annot);

	vector<double> DistAVec;
	vector<double> DistBVec;
	vector<double> DistCVec;
	vector<double> Dist5Vec;
	vector<double> CNVec;

	//for (uint Pos = m_MinPos; Pos <= m_MaxPos; ++Pos)
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double DistA = Chain.GetDist(Pos, Pos_aD);
		double DistB = Chain.GetDist(Pos, Pos_bG);
		double DistC = Chain.GetDist(Pos, Pos_cD);

	// 0..12
		double CN = GetCavityNumber(Pos);
		if (CN > 12)
			CN = 12;
		CN /= 12;
		asserta(CN >= 0 && CN <= 1);

	// 4..16
		double Dist5 = 0;
		if (Pos + 5 < L)
			Dist5 = Chain.GetDist(Pos, Pos+5);
		if (Dist5 < 4)
			Dist5 = 4;
		else if (Dist5 > 16)
			Dist5 = 16;

		Dist5 = (Dist5 - 4)/12;
		asserta(Dist5 >= 0 && Dist5 <= 1);

		DistA /= 30;
		DistB /= 30;
		DistC /= 30;
		DistA = min(DistA, 1.0);
		DistB = min(DistB, 1.0);
		DistC = min(DistC, 1.0);

		DistAVec.push_back(DistA);
		DistBVec.push_back(DistB);
		DistCVec.push_back(DistC);
		Dist5Vec.push_back(Dist5);
		CNVec.push_back(CN);
		}

	GP.AppendFeature("DA", DistAVec);
	GP.AppendFeature("DB", DistBVec);
	GP.AppendFeature("DC", DistCVec);
	GP.AppendFeature("D5", Dist5Vec);
	GP.AppendFeature("CN", CNVec);
	GP.ToTsv(f);
	}

void cmd_struct_prof()
	{
	const string &InputFN = opt_struct_prof;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;

	CMP Prof;
	Prof.FromFile(ModelFileName);

	CMPSearcher CS;
	CS.SetProf(Prof);

	ChainReader CR;
	CR.Open(InputFN);

	PDBChain Chain;
	PDBChain XChain;
	ProgressStep(0, 1001, "Processing");
	uint LastMil = 0;
	uint DoneCount = 0;
	while (CR.GetNext(Chain))
		{
		uint Mil = CR.GetMilDone();
		if (Mil > 0)
			ProgressStep(Mil, 1001, "Processing %u", DoneCount);

		CS.Search(Chain);

		g_APos = UINT_MAX;
		g_BPos = UINT_MAX;
		g_CPos = UINT_MAX;
		double PalmScore = CS.GetPSSMStarts(g_APos, g_BPos, g_CPos);
		if (PalmScore <= 0)
			continue;

		Chain.SetMotifPosVec(g_APos, g_BPos, g_CPos);
		Chain.PrintSeqCoords(g_fLog);
		Chain.SetSS();
		bool Ok = DoStructProf(g_ftsv, CS, Chain);
		++DoneCount;
		if (opt_first_only && Ok)
			break;
		}
	ProgressStep(1000, 1001, "Processing %u", DoneCount);
	}
