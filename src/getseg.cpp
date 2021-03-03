#include "myutils.h"
#include "pssm.h"
#include "sfasta.h"
#include "pssmsearch.h"
#include "xlat.h"

const char *GetSeg(const PSSMHitNt &Hit1, const PSSMHitNt &Hit2,
  int &Frame, unsigned &Lo, unsigned &Hi)
	{
	Frame = 0;
	Lo = UINT_MAX;
	Hi = UINT_MAX;

//	unsigned ColCount1 = Hit1.P->GetColCount();
	unsigned ColCount2 = Hit2.P->GetColCount();

	if (Hit1.NtPos == UINT_MAX || Hit2.NtPos == UINT_MAX ||
	  Hit1.Score < opt_minscore || Hit2.Score < opt_minscore)
		return "(notfound)";

	if (Hit1.Frame != Hit2.Frame)
		return "(frameshift)";

	unsigned Lo2 = Hit1.NtPos;
	unsigned Hi2 = Hit2.NtPos + ColCount2*3 - 1;
	if (Lo2 >= Hi2)
		return "(order)";

	Frame = Hit1.Frame;
	Lo = Lo2;
	Hi = Hi2;
	return Hit1.NtSeq + Lo;
	}
