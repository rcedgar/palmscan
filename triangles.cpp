#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "outputfiles.h"

static double g_MinX;
static double g_MaxX;
static double g_MinY;
static double g_MaxY;

static void WriteSvgHeader()
	{
	if (g_fsvg == 0)
		return;

	double Width = g_MaxX - g_MinX;
	double Height = g_MaxY - g_MinY;

	//Width = 15;
	//Height = 5;

	fprintf(g_fsvg, "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"");
	fprintf(g_fsvg, " width=\"%.6g\" height=\"%.6g\">\n", Width, Height);
	}

static void WriteSvgFooter()
	{
	if (g_fsvg == 0)
		return;
	fprintf(g_fsvg, "</svg>\n");
	}

static void WriteSvgTriangle(const vector<double> &Coords)
	{
	if (g_fsvg == 0)
		return;
	asserta(SIZE(Coords) == 6);

	double A_x = Coords[0] - g_MinX;
	double A_y = Coords[1] - g_MinY;

	double B_x = Coords[2] - g_MinX;
	double B_y = Coords[3] - g_MinY;

	double C_x = Coords[4] - g_MinX;
	double C_y = Coords[5] - g_MinY;

	fprintf(g_fsvg, "<path");
	fprintf(g_fsvg, " stroke-width=\"0.1\"");
	fprintf(g_fsvg, " stroke=\"black\"");
	fprintf(g_fsvg, " fill=\"none\"");
	fprintf(g_fsvg, " d=\"");
	fprintf(g_fsvg, "M %.3g %.3g", A_x, A_y);
	fprintf(g_fsvg, " L %.3g %.3g", B_x, B_y);
	fprintf(g_fsvg, " L %.3g %.3g", C_x, C_y);
	fprintf(g_fsvg, " Z");
	fprintf(g_fsvg, " \"");
	fprintf(g_fsvg, " />");
	fprintf(g_fsvg, "\n");
	}

void cmd_triangles()
	{
	const string &InputFN = opt_triangles;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);

	vector<vector<double> > CoordsVec;
	const uint N = SIZE(Chains);
	double MinX = 0;
	double MinY = 0;
	double MaxX = 0;
	double MaxY = 0;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Writing triangles");
		const PDBChain &Chain = *Chains[i];
		asserta(Chain.m_MotifPosVec.size() == 3);

		vector<vector<double> > MotifCoords;
		Chain.GetMotifCoords(MotifCoords);

		vector<double> Coords;
		if (g_ftsv != 0)
			fprintf(g_ftsv, "%s", Chain.m_Label.c_str());
		for (uint MotifIndex = 0; MotifIndex < 3; ++MotifIndex)
			{
			double PosX = MotifCoords[MotifIndex][X];
			double PosY = MotifCoords[MotifIndex][Y];
			double PosZ = MotifCoords[MotifIndex][Z];
			asserta(PosZ > -0.01 && PosZ < 0.01);
			if (MotifIndex == B)
				asserta(PosY > -0.01 && PosY < 0.01);
			Coords.push_back(PosX);
			Coords.push_back(PosY);

			if (i == 0 && MotifIndex == 0)
				{
				MinX = PosX;
				MaxX = PosX;
				MinY = PosY;
				MaxY = PosY;
				}
			else
				{
				MinX = min(PosX, MinX);
				MaxX = max(PosX, MaxX);
				MinY = min(PosY, MinY);
				MaxY = max(PosY, MaxY);
				}

			if (g_ftsv != 0)
				{
				fprintf(g_ftsv, "%s", Chain.m_Label.c_str());
				fprintf(g_ftsv, "\t%.2f", PosX);
				fprintf(g_ftsv, "\t%.2f", PosY);
				}
			}
		if (g_ftsv != 0)
			fprintf(g_ftsv, "\n");
		CoordsVec.push_back(Coords);
		}

	g_MinX = MinX;
	g_MaxX = MaxX;
	g_MinY = MinY;
	g_MaxY = MaxY;

	WriteSvgHeader();
	const uint M = SIZE(CoordsVec);
	for (uint i = 0; i < M; ++i)
		WriteSvgTriangle(CoordsVec[i]);
	WriteSvgFooter();
	}
