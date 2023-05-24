#include "myutils.h"
#include "pdbchain.h"
#include "spher.h"
#include "abcxyz.h"
#include "cavity.h"

void Palmcator(FILE *f, const PDBChain &Chain)
	{
	if (f == 0)
		return;

	const char *Label = Chain.m_Label.c_str();

	const uint L = Chain.GetSeqLength();

	for (uint Pos = 0; Pos < L; ++Pos)
		{
		vector<double> Pt;
		Chain.GetPt(Pos, Pt);

		double theta, phi, r;
		CartToSpher(Pt, r, theta, phi);
		if (r > g_CavityMaxDist)
			continue;

		double theta_deg = degrees_0_to_360(theta);
		asserta(theta_deg >= 0 && theta_deg < 180);

		double phi_deg = degrees_0_to_360(phi);
		asserta(phi_deg >= 0 && phi_deg < 360);

		char aa = Chain.m_Seq[Pos];
		char MotifChar = Chain.GetMotifChar_Pos(Pos);

		uint ResidueNr = Chain.GetResidueNr(Pos);

		fprintf(f, "%s", Label);
		fprintf(f, "\t%u", ResidueNr);
		fprintf(f, "\t%c", aa);
		fprintf(f, "\t%c", MotifChar);
		fprintf(f, "\t%.3f", r);
		fprintf(f, "\t%.1f", theta_deg);
		fprintf(f, "\t%.1f", phi_deg);
		fprintf(f, "\n");
		}
	}
