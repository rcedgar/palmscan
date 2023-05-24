#include "myutils.h"
#include "trifinder.h"
#include "ppcprofileparams.h"
#include "sort.h"
#include "abcxyz.h"

#define TRACE	0

uint TriFinder::L_ABmin;
uint TriFinder::L_ACmin;
uint TriFinder::L_BCmin;

uint TriFinder::L_ABmax;
uint TriFinder::L_ACmax;
uint TriFinder::L_BCmax;

double TriFinder::d_ABmin;
double TriFinder::d_ACmin;
double TriFinder::d_BCmin;

double TriFinder::d_ABmax;
double TriFinder::d_ACmax;
double TriFinder::d_BCmax;

double GetPPCProfileScore(
  uint LAB, uint LAC, uint LBC,
  double dAB, double dAC, double dBC);

void TriFinder::LogMe(FILE *f) const
	{
	if (f == 0)
		return;

	if (m_Score > 0)
		LogTri(m_PosA, m_PosB, m_PosC);
	else
		fprintf(f, ">%s  (no tri)\n", m_Query->m_Label.c_str());
	}

void TriFinder::SetSigmas(double Sigmas)
	{
	double rN;
	double rd;
#define X(x)	\
    rN = Sigmas*StdDev_L_##x;	\
    rd = Sigmas*StdDev_d_##x;	\
    L_##x##min = uint(floor(Mean_L_##x - rN));	\
    L_##x##max = uint(ceil(Mean_L_##x + rN));	\
    d_##x##min = Mean_d_##x - rd;	\
    d_##x##max = Mean_d_##x + rd;	\
	Log("L_" #x "  %u .. %u\n", L_##x##min, L_##x##max); \
	Log("d_" #x "  %.3g .. %.3g\n", d_##x##min, d_##x##max);	

	X(AB)
	X(AC)
	X(BC)
#undef X
	}

void TriFinder::Find(const PDBChain &Query, bool DGD)
	{
	Clear();
	asserta(L_ABmin != 0);

	m_Query = &Query;

	double MinScore = 0.3;
	const char *Q = m_Query->m_Seq.c_str();
	const uint QL = m_Query->GetSeqLength();
#if TRACE
	Log("\n");
	Log("\n");
	Log("__________________________________________________________________\n");
	Log("TriSearcher::Search() QL=%u\n", QL);
	if (Query.m_MotifPosVec.size() == 3)
		{
		uint PosA = Query.m_MotifPosVec[0];
		uint PosB = Query.m_MotifPosVec[1];
		uint PosC = Query.m_MotifPosVec[2];
		Log("PPC:");
		LogTri(PosA, PosB, PosC);
		}
#endif
	for (uint PosA = 0; PosA < QL; ++PosA)
		{
		if (DGD && Q[PosA+3] != 'D')
			continue;
		if (PosA + L_ACmin > QL)
			{
#if TRACE
			Log("PosA=%u + NACmin > QL\n", PosA);
#endif
			break;
			}

		for (uint PosC = PosA + L_ACmin; PosC <= PosA + L_ACmax; ++PosC)
			{
			if (DGD && Q[PosC+3] != 'D')
				continue;
			if (PosC + CL > QL)
				{
#if TRACE
				Log("PosC+CL=%u > QL\n", PosC+CL);
#endif
				break;
				}

			double dAC = m_Query->GetDist(PosA, PosC);
#if TRACE
			Log("PosA=%u, PosC=%u, dAC=%.1f\n", PosA, PosC, dAC);
#endif
			if (dAC < d_ACmin || dAC > d_ACmax)
				{
#if TRACE
				Log("  dAC=%.3g outside radius\n", dAC);
#endif
				continue;
				}

			for (uint PosB = PosA + L_ABmin; PosA + L_ABmax; ++PosB)
				{
				if (DGD && Q[PosB+1] != 'G')
					continue;
				if (PosB >= PosC)
					break;
				uint L_BC = PosC - PosB + 1;
				if (L_BC < L_BCmin || L_BC > L_BCmax)
					continue;

				double dAB = m_Query->GetDist(PosA, PosB);
#if TRACE
				Log("PosA=%u, PosB=%u, dAB=%.1f\n", PosA, PosB, dAB);
#endif
				if (dAB < d_ABmin || dAB > d_ABmax)
					continue;

				double dBC = m_Query->GetDist(PosB, PosC);
#if TRACE
				Log("PosB=%u, PosC=%u, dBC=%.1f\n", PosB, PosC, dBC);
#endif
#if TRACE
				{
				uint L_AB = PosB - PosA + 1;
				uint L_AC = PosC - PosA + 1;
				LogTri(PosA, PosB, PosC);
				}
#endif
				if (dBC < d_BCmin || dBC > d_BCmax)
					continue;

				if (PosA + AL >= PosB)
					continue;
				if (PosB + BL >= PosC)
					continue;

				uint L_AB = PosB - PosA + 1;
				uint L_AC = PosC - PosA + 1;
				double Score =
				  GetPPCProfileScore(L_AB, L_AC, L_BC, dAB, dAC, dBC);
				if (Score > MinScore)
					{
#if TRACE
					LogTri(PosA, PosB, PosC);
#endif
					++m_TriCount;
					if (Score > m_Score)
						{
						m_Score = Score;
						m_PosA = PosA;
						m_PosB = PosB;
						m_PosC = PosC;
						}
					}
				}
			}
		}
	}

void TriFinder::LogTri(uint PosA, uint PosB, uint PosC) const
	{
	uint L_AB = PosB - PosA + 1;
	uint L_AC = PosC - PosA + 1;
	uint L_BC = PosC - PosB + 1;

	double dAB = m_Query->GetDist(PosA, PosB);
	double dAC = m_Query->GetDist(PosA, PosC);
	double dBC = m_Query->GetDist(PosB, PosC);

	bool dAB_ok = (dAB >= d_ABmin || dAB <= d_ABmax);
	bool dAC_ok = (dAC >= d_ACmin || dAC <= d_ACmax);
	bool dBC_ok = (dBC >= d_BCmin || dBC <= d_BCmax);

	bool L_AB_ok = (L_AB >= L_ABmin || L_AB <= L_ABmax);
	bool L_AC_ok = (L_AC >= L_ACmin || L_AC <= L_ACmax);
	bool L_BC_ok = (L_BC >= L_BCmin || L_BC <= L_BCmax);

	double Score = GetPPCProfileScore(L_AB, L_AC, L_BC, dAB, dAC, dBC);

	Log("\n>%s  Score(%u %c, %u %c, %u %c, %.3g %c, %.3g %c, %.3g %c) = %.4f [%u tris]\n",
	  m_Query->m_Label.c_str(),
	  L_AB, yon(L_AB_ok),
	  L_AC, yon(L_AC_ok),
	  L_BC, yon(L_BC_ok),
	  dAB, yon(dAB_ok),
	  dAC, yon(dAC_ok),
	  dBC, yon(dBC_ok),
	  Score,
	  m_TriCount);
	}
