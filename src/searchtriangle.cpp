#include "myutils.h"
#include "pdbchain.h"
#include "trisearcher.h"
#include "searchparams.h"
#include "abcxyz.h"

#define TRACE	0

void TriSearcher::Search(PDBChain &Query, const PDBChain &Ref)
	{
	m_Query = &Query;
	m_Ref = &Ref;

	double RefAB = 0;
	double RefBC = 0;
	double RefAC = 0;
	m_Ref->GetMotifDists(RefAB, RefBC, RefAC);
//	m_Ref->GetMotifDists2(RefAB, RefBC, RefAC);

	vector<vector<double> > MotifCoords;
	//m_Ref->GetMotifCoords(MotifCoords);
	//LogMx("MotifCoords", MotifCoords);

	asserta(NABmin <= NABmax);
	asserta(NBCmin <= NBCmax);
	asserta(NACmin <= NACmax);

	m_PosAs.clear();
	m_PosBs.clear();
	m_PosCs.clear();
	m_TriRMSD2s.clear();
	m_MotifRMSD2s.clear();

	const uint QL = m_Query->GetSeqLength();
#if TRACE
	Log("\n");
	Log("TriSearcher::Search() QL=%u\n", QL);
	LogMe();
	Log("RefAB %.1f, RefBC %.1f, RefAC %.1f\n", RefAB, RefBC, RefAC);
#endif
	for (uint PosA = 0; PosA < QL; ++PosA)
		{
		if (PosA + NACmin > QL)
			{
#if TRACE
			Log("PosA=%u + NACmin > QL\n", PosA);
#endif
			break;
			}

		for (uint PosC = PosA + NACmin; PosC <= PosA + NACmax; ++PosC)
			{
			if (PosC + CL > QL)
				{
#if TRACE
				Log("PosC+CL=%u > QL\n", PosC+CL);
#endif
				break;
				}

			double QueryAC = m_Query->GetDist(PosA, PosC);
#if TRACE
			Log("PosA=%u, PosC=%u, QueryAC=%.1f RefAC=%.1f\n",
			  PosA, PosC, QueryAC, RefAC);
#endif
			if (QueryAC < RefAC - Radius || QueryAC > RefAC + Radius)
				{
#if TRACE
				Log("  outside radius\n");
#endif
				continue;
				}

			for (uint PosB = PosA + NABmin; PosA + NABmax; ++PosB)
				{
				if (PosB >= PosC)
					break;
				uint NBC = PosC - PosB;
				if (NBC < NBCmin || NBC > NBCmax)
					continue;

				double QueryAB = m_Query->GetDist(PosA, PosB);
#if TRACE
				Log("PosA=%u, PosB=%u, QueryAB=%.1f RefAB=%.1f\n",
				  PosA, PosB, QueryAB, RefAB);
#endif
				if (QueryAB < RefAB - Radius || QueryAB > RefAB + Radius)
					continue;

				double QueryBC = m_Query->GetDist(PosB, PosC);
#if TRACE
				Log("PosB=%u, PosC=%u, QueryBC=%.1f RefBC=%.1f\n",
				  PosB, PosC, QueryBC, RefBC);
#endif
				if (QueryBC < RefBC - Radius || QueryBC > RefBC + Radius)
					continue;

				if (PosA + AL >= PosB)
					continue;
				if (PosB + BL >= PosC)
					continue;
				double DistAB = RefAB - QueryAB;
				double DistBC = RefBC - QueryBC;
				double DistAC = RefAC - QueryAC;
				double TriRMSD2 = DistAB*DistAB + DistBC*DistBC + DistAC*DistAC;
				if (TriRMSD2 < MaxTriRMSD2)
					{
					double MotifRMSD2 = GetRMSDMotifs(PosA, PosB, PosC);
					if (MotifRMSD2 < MaxMotifRMSD)
						{
						m_TriRMSD2s.push_back(TriRMSD2);
						m_MotifRMSD2s.push_back(MotifRMSD2);
						m_PosAs.push_back(PosA);
						m_PosBs.push_back(PosB);
						m_PosCs.push_back(PosC);
						}
					}
				}
			}
		}
	SetHitOrder();
	}
