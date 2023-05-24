#include "myutils.h"
#include "outputfiles.h"
#include "pdbchain.h"
#include "abcxyz.h"

void PDBChain::ToPML(FILE *f, const string &PDBFileName) const
	{
	if (f == 0)
		return;
	if (PDBFileName == "")
		Die("ToPML() PDBFileName empty");

	uint PosA = GetMotifPos(A);
	uint PosB = GetMotifPos(B);
	uint PosC = GetMotifPos(C);

	int ResLoA;
	int ResLoB;
	int ResLoC;

	int ResHiA;
	int ResHiB;
	int ResHiC;

	GetResidueRange(PosA, AL, ResLoA, ResHiA);
	GetResidueRange(PosB, BL, ResLoB, ResHiB);
	GetResidueRange(PosC, CL, ResLoC, ResHiC);

	fprintf(f, "#!/usr/bin/pymol\n");
	fprintf(f, "\n");
	fprintf(f, "cmd.load(\"%s\")\n", PDBFileName.c_str());

	string Label;
	GetLabelFromFileName(PDBFileName, Label);
	const char *Lab = Label.c_str();

	fprintf(f, "\n");
	fprintf(f, "select %s\n", Lab);
	fprintf(f, "color gray20, %s\n", Lab);

	fprintf(f, "\n");
	fprintf(f, "cmd.select(\"%s_PP\", \"%s and resi %u-%d\")\n", Lab, Lab, ResLoA, ResHiC);

	fprintf(f, "cmd.select(\"%s_mA\", \"%s and resi %u-%d\")\n", Lab, Lab, ResLoA, ResHiA);
	fprintf(f, "cmd.select(\"%s_mB\", \"%s and resi %u-%d\")\n", Lab, Lab, ResLoB, ResHiB);
	fprintf(f, "cmd.select(\"%s_mC\", \"%s and resi %u-%d\")\n", Lab, Lab, ResLoC, ResHiC);

	fprintf(f, "cmd.select(\"%s_aD\", \"%s and name CA and resi %u\")\n", Lab, Lab, ResLoA+3);
	fprintf(f, "cmd.select(\"%s_bG\", \"%s and name CA and resi %u\")\n", Lab, Lab, ResLoB+1);
	fprintf(f, "cmd.select(\"%s_cD\", \"%s and name CA and resi %u\")\n", Lab, Lab, ResLoC+3);

	fprintf(f, "\n");
	fprintf(f, "cmd.select(\"%s_DGD\", \"%s_aD + %s_bG + %s_cD\")\n", Lab, Lab, Lab, Lab);
	fprintf(f, "show sphere, %s_DGD\n", Lab);
	fprintf(f, "set sphere_transparency, 0.3\n");

	fprintf(f, "\n");
	fprintf(f, "color gray60, %s_PP\n", Lab);
	fprintf(f, "color tv_blue, %s_mA\n", Lab);
	fprintf(f, "color tv_green, %s_mB\n", Lab);
	fprintf(f, "color tv_red, %s_mC\n", Lab);
	fprintf(f, "\n");
	fprintf(f, "deselect\n");
	fprintf(f, "delete sele\n");

	fprintf(f, "\n");
	}
